"""OpenMC transport operator

This module implements a transport operator for OpenMC so that it can be used by
depletion integrators. The implementation makes use of the Python bindings to
OpenMC's C API so that reading tally results and updating material number
densities is all done in-memory instead of through the filesystem.

"""

import sys
import copy
from collections import OrderedDict
from itertools import chain
import os
import time
import xml.etree.ElementTree as ET
from warnings import warn

import h5py
import numpy as np
import math 
from uncertainties import ufloat

import openmc
from openmc import zernike as zer 
import openmc.lib
from . import comm
from .abc import TransportOperator, OperatorResult
from .atom_number import AtomNumber
from .reaction_rates import ReactionRates
from .results_list import ResultsList
from .helpers import (
    DirectReactionRateHelper, ChainFissionHelper, ConstantFissionYieldHelper,
    FissionYieldCutoffHelper, AveragedFissionYieldHelper, EnergyScoreHelper)


__all__ = ["Operator", "OperatorResult"]


def _distribute(items):
    """Distribute items across MPI communicator

    Parameters
    ----------
    items : list
        List of items of distribute

    Returns
    -------
    list
        Items assigned to process that called

    """
    min_size, extra = divmod(len(items), comm.size)
    j = 0
    for i in range(comm.size):
        chunk_size = min_size + int(i < extra)
        if comm.rank == i:
            return items[j:j + chunk_size]
        j += chunk_size


class Operator(TransportOperator):
    """OpenMC transport operator for depletion.

    Instances of this class can be used to perform depletion using OpenMC as the
    transport operator. Normally, a user needn't call methods of this class
    directly. Instead, an instance of this class is passed to an integrator
    class, such as :class:`openmc.deplete.CECMIntegrator`.

    Parameters
    ----------
    geometry : openmc.Geometry
        OpenMC geometry object
    settings : openmc.Settings
        OpenMC Settings object
    chain_file : str, optional
        Path to the depletion chain XML file.  Defaults to the file
        listed under ``depletion_chain`` in
        :envvar:`OPENMC_CROSS_SECTIONS` environment variable.
    prev_results : ResultsList, optional
        Results from a previous depletion calculation. If this argument is
        specified, the depletion calculation will start from the latest state
        in the previous results.
    diff_burnable_mats : bool, optional
        Whether to differentiate burnable materials with multiple instances.
        Volumes are divided equally from the original material volume.
        Default: False.
    energy_mode : {"energy-deposition", "fission-q"}
        Indicator for computing system energy. ``"energy-deposition"`` will
        compute with a single energy deposition tally, taking fission energy
        release data and heating into consideration. ``"fission-q"`` will
        use the fission Q values from the depletion chain
    fission_q : dict, optional
        Dictionary of nuclides and their fission Q values [eV]. If not given,
        values will be pulled from the ``chain_file``. Only applicable
        if ``"energy_mode" == "fission-q"``
    dilute_initial : float, optional
        Initial atom density [atoms/cm^3] to add for nuclides that are zero
        in initial condition to ensure they exist in the decay chain.
        Only done for nuclides with reaction rates.
        Defaults to 1.0e3.
    fission_yield_mode : {"constant", "cutoff", "average"}
        Key indicating what fission product yield scheme to use. The
        key determines what fission energy helper is used:

        * "constant": :class:`~openmc.deplete.helpers.ConstantFissionYieldHelper`
        * "cutoff": :class:`~openmc.deplete.helpers.FissionYieldCutoffHelper`
        * "average": :class:`~openmc.deplete.helpers.AveragedFissionYieldHelper`

        The documentation on these classes describe their methodology
        and differences. Default: ``"constant"``
    fission_yield_opts : dict of str to option, optional
        Optional arguments to pass to the helper determined by
        ``fission_yield_mode``. Will be passed directly on to the
        helper. Passing a value of None will use the defaults for
        the associated helper.
    fet_order : integer, optional
        The order of functional expansion tallied (fet) for nuclide 
        number density. 
        Default: 0        

    Attributes
    ----------
    geometry : openmc.Geometry
        OpenMC geometry object
    settings : openmc.Settings
        OpenMC settings object
    dilute_initial : float
        Initial atom density [atoms/cm^3] to add for nuclides that
        are zero in initial condition to ensure they exist in the decay
        chain. Only done for nuclides with reaction rates.
    output_dir : pathlib.Path
        Path to output directory to save results.
    round_number : bool
        Whether or not to round output to OpenMC to 8 digits.
        Useful in testing, as OpenMC is incredibly sensitive to exact values.
    number : openmc.deplete.AtomNumber
        Total number of atoms in simulation.
    nuclides_with_data : set of str
        A set listing all unique nuclides available from cross_sections.xml.
    chain : openmc.deplete.Chain
        The depletion chain information necessary to form matrices and tallies.
    reaction_rates : openmc.deplete.ReactionRates
        Reaction rates from the last operator step.
    burnable_mats : list of str
        All burnable material IDs
    heavy_metal : float
        Initial heavy metal inventory [g]
    local_mats : list of str
        All burnable material IDs being managed by a single process
    prev_res : ResultsList or None
        Results from a previous depletion calculation. ``None`` if no
        results are to be used.
    diff_burnable_mats : bool
        Whether to differentiate burnable materials with multiple instances
    """
    _fission_helpers = {
        "average": AveragedFissionYieldHelper,
        "constant": ConstantFissionYieldHelper,
        "cutoff": FissionYieldCutoffHelper,
    }

    def __init__(self, geometry, settings, chain_file=None, prev_results=None,
                 diff_burnable_mats=False, energy_mode="fission-q",
                 fission_q=None, dilute_initial=1.0e3,
                 fission_yield_mode="constant", fission_yield_opts=None):
        if fission_yield_mode not in self._fission_helpers:
            raise KeyError(
                "fission_yield_mode must be one of {}, not {}".format(
                    ", ".join(self._fission_helpers), fission_yield_mode))
        if energy_mode == "energy-deposition":
            if fission_q is not None:
                warn("Fission Q dictionary not used if energy deposition "
                     "is used")
                fission_q = None
        elif energy_mode != "fission-q":
            raise ValueError(
                "energy_mode {} not supported. Must be energy-deposition "
                "or fission-q".format(energy_mode))
        super().__init__(chain_file, fission_q, dilute_initial, prev_results)
        self.round_number = False
        self.prev_res = None
        self.settings = settings
        self.geometry = geometry
        self.diff_burnable_mats = diff_burnable_mats 
        # FETs 
        self.fet_deplete = None
        if settings.fet_deplete is not None and settings.fet_deplete['enable'] == True:
            self.fet_deplete = settings.fet_deplete

        # Differentiate burnable materials with multiple instances
        if self.diff_burnable_mats:
            self._differentiate_burnable_mats()

        # Clear out OpenMC, create task lists, distribute
        openmc.reset_auto_ids()
        self.burnable_mats, volume, nuclides = self._get_burnable_mats()
        self.local_mats = _distribute(self.burnable_mats)
        
        
        # Generate map from local materials => material index
        self._mat_index_map = {
            lm: self.burnable_mats.index(lm) for lm in self.local_mats}

        if self.prev_res is not None:
            # Reload volumes into geometry
            prev_results[-1].transfer_volumes(geometry)

            # Store previous results in operator
            # Distribute reaction rates according to those tracked
            # on this process
            if comm.size == 1:
                self.prev_res = prev_results
            else:
                self.prev_res = ResultsList()
                mat_indexes = _distribute(range(len(self.burnable_mats)))
                for res_obj in prev_results:
                    new_res = res_obj.distribute(self.local_mats, mat_indexes)
                    self.prev_res.append(new_res)

        # Determine which nuclides have incident neutron data
        self.nuclides_with_data = self._get_nuclides_with_data()

        # Select nuclides with data that are also in the chain
        self._burnable_nucs = [nuc.name for nuc in self.chain.nuclides
                               if nuc.name in self.nuclides_with_data]

        # Extract number densities from the geometry / previous depletion run
        self._extract_number(self.local_mats, volume, nuclides, self.prev_res, # FETs 
                             fet_deplete=self.fet_deplete)

        # Create reaction rates array
        self.reaction_rates = ReactionRates(
            self.local_mats, self._burnable_nucs, self.chain.reactions,
            fet_deplete=self.fet_deplete) # FETs 
        
        # Get classes to assist working with tallies
        self._rate_helper = DirectReactionRateHelper(
            self.reaction_rates.n_nuc, self.reaction_rates.n_react, 
            fet_deplete=self.fet_deplete) #FETs    
        if energy_mode == "fission-q":
            self._energy_helper = ChainFissionHelper()
        else:
            score = "heating" if settings.photon_transport else "heating-local"
            self._energy_helper = EnergyScoreHelper(score)
        
        # Select and create fission yield helper
        fission_helper = self._fission_helpers[fission_yield_mode]
        fission_yield_opts = (
            {} if fission_yield_opts is None else fission_yield_opts)
        self._yield_helper = fission_helper.from_operator(
            self, **fission_yield_opts)
        

    def __call__(self, vec, power):
        """Runs a simulation.

        Simulation will abort under the following circumstances:

            1) No energy is computed using OpenMC tallies.

        Parameters
        ----------
        vec : list of numpy.ndarray
            Total atoms to be used in function.
        power : float
            Power of the reactor in [W]

        Returns
        -------
        openmc.deplete.OperatorResult
            Eigenvalue and reaction rates resulting from transport operator

        """
        # Prevent OpenMC from complaining about re-creating tallies
        openmc.reset_auto_ids()

        # Update status
        # print(vec) # FETs testing 
        # print(self.number[0, :])
        # print('O16', self.number.get_atom_density(0, 'O16', fet_deplete=self.fet_deplete))
        self.number.set_density(vec, self.fet_deplete) # FETs 
        # print(self.number[0, :])
        # print('O16', self.number.get_atom_density(0, 'O16', fet_deplete=self.fet_deplete))
        # Update material compositions and tally nuclides
        self._update_materials() # FETs 
        nuclides = self._get_tally_nuclides(fet_deplete=self.fet_deplete) # FETs 
        self._rate_helper.nuclides = nuclides # FETs 
        self._energy_helper.nuclides = nuclides # FETs no change
        self._yield_helper.update_tally_nuclides(nuclides) #FETs no change

        # Run OpenMC
        openmc.lib.reset()
        openmc.lib.run()

        # Extract results
        op_result = self._unpack_tallies_and_normalize(power, fet_deplete=self.fet_deplete)

        return copy.deepcopy(op_result)

    @staticmethod
    def write_bos_data(step):
        """Write a state-point file with beginning of step data

        Parameters
        ----------
        step : int
            Current depletion step including restarts
        """
        openmc.lib.statepoint_write(
            "openmc_simulation_n{}.h5".format(step),
            write_source=False)

    def _differentiate_burnable_mats(self):
        """Assign distribmats for each burnable material

        """

        # Count the number of instances for each cell and material
        self.geometry.determine_paths(instances_only=True)

        # Extract all burnable materials which have multiple instances
        distribmats = set(
            [mat for mat in self.geometry.get_all_materials().values()
             if mat.depletable and mat.num_instances > 1])

        if distribmats:
            # Assign distribmats to cells
            for cell in self.geometry.get_all_material_cells().values():
                if cell.fill in distribmats and cell.num_instances > 1:
                    mat = cell.fill
                    if mat.volume is None:
                        raise RuntimeError("Volume not specified for depletable "
                                           "material with ID={}.".format(mat.id))
                    mat.volume /= mat.num_instances
                    cell.fill = [mat.clone()
                                 for i in range(cell.num_instances)]

    def _get_burnable_mats(self):
        """Determine depletable materials, volumes, and nuclides

        Returns
        -------
        burnable_mats : list of str
            List of burnable material IDs
        volume : OrderedDict of str to float
            Volume of each material in [cm^3]
        nuclides : list of str
            Nuclides in order of how they'll appear in the simulation.

        """

        burnable_mats = set()
        model_nuclides = set()
        volume = OrderedDict()

        self.heavy_metal = 0.0

        # Iterate once through the geometry to get dictionaries
        for mat in self.geometry.get_all_materials().values():
            for nuclide in mat.get_nuclides():
                model_nuclides.add(nuclide)
            if mat.depletable:
                burnable_mats.add(str(mat.id))
                if mat.volume is None:
                    raise RuntimeError("Volume not specified for depletable "
                                       "material with ID={}.".format(mat.id))
                volume[str(mat.id)] = mat.volume
                self.heavy_metal += mat.fissionable_mass

        # Make sure there are burnable materials
        if not burnable_mats:
            raise RuntimeError(
                "No depletable materials were found in the model.")

        # Sort the sets
        burnable_mats = sorted(burnable_mats, key=int)
        model_nuclides = sorted(model_nuclides)

        # Construct a global nuclide dictionary, burned first
        nuclides = list(self.chain.nuclide_dict)
        for nuc in model_nuclides:
            if nuc not in nuclides:
                nuclides.append(nuc)

        return burnable_mats, volume, nuclides

    def _extract_number(self, local_mats, volume, nuclides, prev_res=None, fet_deplete=None):
        """Construct AtomNumber using geometry

        Parameters
        ----------
        local_mats : list of str
            Material IDs to be managed by this process
        volume : OrderedDict of str to float
            Volumes for the above materials in [cm^3]
        nuclides : list of str
            Nuclides to be used in the simulation.
        prev_res : ResultsList, optional
            Results from a previous depletion calculation

        """
        self.number = AtomNumber(local_mats, nuclides, volume, len(self.chain), fet_deplete=fet_deplete)

        if self.dilute_initial != 0.0:
            #print(self.dilute_initial) #FETs 
            #print(self._burnable_nucs)
            for nuc in self._burnable_nucs:
                self.number.set_atom_density(np.s_[:], nuc, self.dilute_initial, fet_deplete=fet_deplete)

        # Now extract and store the number densities
        # From the geometry if no previous depletion results
        if prev_res is None:
            for mat in self.geometry.get_all_materials().values():
                #print(local_mats)
                if str(mat.id) in local_mats:
                    self._set_number_from_mat(mat, fet_deplete=fet_deplete) #FETs 

        # Else from previous depletion results
        else:
            for mat in self.geometry.get_all_materials().values():
                if str(mat.id) in local_mats:
                    self._set_number_from_results(mat, prev_res, fet_deplete=fet_deplete) #FETs 
        #print(self.number.number) # FETs testing 
        
    def _set_number_from_mat(self, mat, fet_deplete=None):
        """Extracts material and number densities from openmc.Material

        Parameters
        ----------
        mat : openmc.Material
            The material to read from

        """
        import itertools
        
        mat_id = str(mat.id)
        # FETs 
        for nuclide, density, *temp_args in mat.get_nuclide_atom_densities().values():
            if len(temp_args) > 0:
                coeff = []
                temp_args = list(itertools.chain(*temp_args))[0]
                temp_args = temp_args[1:]
                #temp_args[0] = density * 1.0e24
                coeff = np.array(temp_args) 
                coeff[0] *= 1.0e24
                #print(mat_id, nuclide, coeff)
                self.number.set_atom_density(mat_id, nuclide, coeff, fet_deplete=fet_deplete)
            else:
                number = density * 1.0e24
                self.number.set_atom_density(mat_id, nuclide, number)

    def _set_number_from_results(self, mat, prev_res, fet_deplete=None):
        """Extracts material nuclides and number densities.

        If the nuclide concentration's evolution is tracked, the densities come
        from depletion results. Else, densities are extracted from the geometry
        in the summary.

        Parameters
        ----------
        mat : openmc.Material
            The material to read from
        prev_res : ResultsList
            Results from a previous depletion calculation

        """
        mat_id = str(mat.id)

        # Get nuclide lists from geometry and depletion results
        depl_nuc = prev_res[-1].nuc_to_ind
        geom_nuc_densities = mat.get_nuclide_atom_densities()

        # Merge lists of nuclides, with the same order for every calculation
        geom_nuc_densities.update(depl_nuc)

        for nuclide in geom_nuc_densities.keys():
            if nuclide in depl_nuc:
                concentration = prev_res.get_atoms(mat_id, nuclide, fet_deplete=fet_deplete)[1][-1] #FETs 
                #print(concentration)
                volume = prev_res[-1].volume[mat_id]
                number = concentration / volume
            else:
                density = geom_nuc_densities[nuclide][1]
                number = density * 1.0e24

            self.number.set_atom_density(mat_id, nuclide, number, fet_deplete=fet_deplete) #FETs 

    def initial_condition(self):
        """Performs final setup and returns initial condition.

        Returns
        -------
        list of numpy.ndarray
            Total density for initial conditions.
        """

        # Create XML files
        if comm.rank == 0:
            self.geometry.export_to_xml()
            self.settings.export_to_xml()
            self._generate_materials_xml()
            

        # Initialize OpenMC library
        comm.barrier()
        openmc.lib.init(intracomm=comm)

        # Generate tallies in memory
        materials = [openmc.lib.materials[int(i)]
                     for i in self.burnable_mats]
        # FETs 
        #nuclides = self._burnable_nucs
        self._rate_helper.generate_tallies(materials, self.chain.reactions, 
                                           fet_deplete=self.fet_deplete)
        self._energy_helper.prepare(
            self.chain.nuclides, self.reaction_rates.index_nuc, materials)
        # Tell fission yield helper what materials this process is
        # responsible for
        self._yield_helper.generate_tallies(
            materials, tuple(sorted(self._mat_index_map.values())))

        # Return number density vector
        return list(self.number.get_mat_slice(np.s_[:], fet_deplete=self.fet_deplete)) # FETs 

    def finalize(self):
        """Finalize a depletion simulation and release resources."""
        openmc.lib.finalize()

    def _update_materials(self):
        """Updates material compositions in OpenMC on all processes."""
        
        from collections import Iterable # FETs 
        
        for rank in range(comm.size):
            number_i = comm.bcast(self.number, root=rank)

            for mat in number_i.materials:
                nuclides = []
                densities = []
                for nuc in number_i.nuclides:
                    if nuc in self.nuclides_with_data:
                        val = number_i.get_atom_density(mat, nuc, fet_deplete=self.fet_deplete) #FETs 
                        #if self.fet_deplete is not None:
                        #    val = self._check_negative(val, fet_deplete=self.fet_deplete) 
                        #    number_i.set_atom_density(mat, nuc, val, fet_deplete=self.fet_deplete) 
                        #    print(mat, nuc, val)
                        if not isinstance(val, Iterable):
                            val *= 1.0e-24  
                            # If nuclide is zero, do not add to the problem.
                            if val > 1.0e-50:
                                if self.round_number:
                                    val_magnitude = np.floor(np.log10(val))
                                    val_scaled = val / 10**val_magnitude
                                    val_round = round(val_scaled, 8)
                            
                                    val = val_round * 10**val_magnitude
                            
                                nuclides.append(nuc)
                                densities.append(val)
                            else:
                                # Only output warnings if values are significantly
                                # negative. CRAM does not guarantee positive values.
                                if val <= 1.0e-50:
                                    print("WARNING: nuclide ", nuc, " in material ", mat,
                                          " is negative (density = ", val, " at/barn-cm)")
                                number_i[mat, nuc] = 0.0
                        else: #FETs 
                            val[0] *= 1.0e-24
                            if val[0] > 1.0e-50:
                                if self.round_number:
                                    val_magnitude = np.floor(np.log10(val[0]))
                                    val_scaled = val[0] / 10**val_magnitude
                                    val_round = round(val_scaled, 8)
                            
                                    val[0] = val_round * 10**val_magnitude
                            
                                nuclides.append(nuc)
                                densities.append(val[0])
                            else:
                                # Only output warnings if values are significantly
                                # negative. CRAM does not guarantee positive values.
                                if val[0] <= 1.0e-50:
                                    print("WARNING: nuclide ", nuc, " in material ", mat,
                                          " is negative (density = ", val[0], " at/barn-cm)")
                                number_i[mat, nuc] = 0.0

                # Update densities on C API side
                mat_internal = openmc.lib.materials[int(mat)]
                mat_internal.set_densities(nuclides, densities)
                #TODO Update densities on the Python side, otherwise the
                # summary.h5 file contains densities at the first time step    
    
    def _check_negative(self, val,  fet_deplete=None):
        """
        Check the nagative of val and force it to positive 
        """
        mp = 1
        if fet_deplete is not None:
            order = fet_deplete['order']
            radius = fet_deplete['radius']
            if fet_deplete['name']== 'zernike':
                mp = zer.num_poly(fet_deplete['order'])
            elif fet['name']=='zernike1d':
                mp = zer.num_poly1d(fet_deplete['order'])
        #
        fet = zer.ZernikePolynomial(order, val, radius, sqrt_normed=False)
        fet.force_positive()
        return fet.coeffs
    
    def _export_materials_xml(self):
        """
        """
        if comm.rank == 0:
            materials = openmc.Materials(self.geometry.get_all_materials()
                                         .values())
            number = self.number
            nuclides = list(self.number.nuclides)        
            for i in range(len(materials)):
                materials[i]._nuclides.sort(key=lambda x: nuclides.index(x[0]))
                mat = materials[i]
                for mat_i in number.materials:
                    if str(mat.id) == mat_i:
                        # change density unit to "sum"
                        materials[i]._density_units = "sum"
                        #for j in range(len(mat.nuclides)):
                        #    nuc = mat.nuclides[j]
                        #    nuc_name = mat.nuclides[j][0]
                        #    val = number.get_atom_density(str(mat.id), nuc_name, fet_deplete=self.fet_deplete) 
                        #    val /= 1.0e24 # Unit conversion from atom/cm3 to atom/b-cm
                        #    materials[i].update_nuclide(nuc_name, val, fet_deplete=self.fet_deplete)
                        for nuc_name in number.nuclides:
                            val = number.get_atom_density(str(mat.id), nuc_name, fet_deplete=self.fet_deplete) 
                            val /= 1.0e24 # Unit conversion from atom/cm3 to atom/b-cm
                            loc_nuc = materials[i].find_nuclide(nuc_name)
                            #print(i, nuc_name, loc_nuc, val)
                            if loc_nuc == -1:
                                if self.fet_deplete is None:
                                    if val > 1.0e-50: # val[0], val.any() or val.all()
                                        materials[i].add_nuclide(nuc_name, val)
                                else:
                                    if val.any() > 1.0e-50:
                                        materials[i].add_nuclide_fet(nuc_name, val, fet_deplete=self.fet_deplete)
                            else:
                                #print(nuc_name, val)
                                materials[i].update_nuclide(nuc_name, val, fet_deplete=self.fet_deplete)
                        #break                      
            materials.export_to_xml()
        
        
    def _generate_materials_xml(self):
        """Creates materials.xml from self.number.

        Due to uncertainty with how MPI interacts with OpenMC API, this
        constructs the XML manually.  The long term goal is to do this
        through direct memory writing.

        """
        materials = openmc.Materials(self.geometry.get_all_materials()
                                     .values())
        #print(materials)
        # Sort nuclides according to order in AtomNumber object
        nuclides = list(self.number.nuclides)
        for mat in materials:
            mat._nuclides.sort(key=lambda x: nuclides.index(x[0]))

        materials.export_to_xml()

    def _get_tally_nuclides(self, fet_deplete=None):
        """Determine nuclides that should be tallied for reaction rates.

        This method returns a list of all nuclides that have neutron data and
        are listed in the depletion chain. Technically, we should tally nuclides
        that may not appear in the depletion chain because we still need to get
        the fission reaction rate for these nuclides in order to normalize
        power, but that is left as a future exercise.

        Returns
        -------
        list of str
            Tally nuclides

        """
        #FETs
        mp = 1
        if fet_deplete is not None:
            if fet_deplete['name']== 'zernike':
                mp = zer.num_poly(fet_deplete['order'])
            elif fet['name']=='zernike1d':
                mp = zer.num_poly1d(fet_deplete['order'])
        #
        nuc_set = set()

        # Create the set of all nuclides in the decay chain in materials marked
        # for burning in which the number density is greater than zero.
        for nuc in self.number.nuclides:
            if nuc in self.nuclides_with_data:
                #FETs                 
                if fet_deplete is None:
                    if np.sum(self.number[:, nuc]) > 0.0:
                        nuc_set.add(nuc)
                else: #FETs
                    i_nuc = self.number.index_nuc[nuc] 
                    if np.sum(self.number[:, i_nuc * mp]) > 0.0:
                        nuc_set.add(nuc)                

        # Communicate which nuclides have nonzeros to rank 0
        if comm.rank == 0:
            for i in range(1, comm.size):
                nuc_newset = comm.recv(source=i, tag=i)
                nuc_set |= nuc_newset
        else:
            comm.send(nuc_set, dest=0, tag=comm.rank)

        if comm.rank == 0:
            # Sort nuclides in the same order as self.number
            nuc_list = [nuc for nuc in self.number.nuclides
                        if nuc in nuc_set]
        else:
            nuc_list = None

        # Store list of tally nuclides on each process
        nuc_list = comm.bcast(nuc_list)
        return [nuc for nuc in nuc_list if nuc in self.chain]

    def _unpack_tallies_and_normalize(self, power, fet_deplete=None):
        """Unpack tallies from OpenMC and return an operator result

        This method uses OpenMC's C API bindings to determine the k-effective
        value and reaction rates from the simulation. The reaction rates are
        normalized by the user-specified power, summing the product of the
        fission reaction rate times the fission Q value for each material.

        Parameters
        ----------
        power : float
            Power of the reactor in [W]

        Returns
        -------
        openmc.deplete.OperatorResult
            Eigenvalue and reaction rates resulting from transport operator

        """
        # FETs 
        mp = 1
        if fet_deplete is not None:
            if fet_deplete['name']== 'zernike':
                mp = zer.num_poly(fet_deplete['order'])
            elif fet['name']=='zernike1d':
                mp = zer.num_poly1d(fet_deplete['order'])
        #
        rates = self.reaction_rates
        rates.fill(0.0)

        # Get k and uncertainty
        k_combined = ufloat(*openmc.lib.keff())

        # Extract tally bins
        nuclides = self._rate_helper.nuclides

        # Form fast map
        nuc_ind = [rates.index_nuc[nuc] for nuc in nuclides]
        react_ind = [rates.index_rx[react] for react in self.chain.reactions]

        # Compute fission power

        # Keep track of energy produced from all reactions in eV per source
        # particle
        self._energy_helper.reset()
        self._yield_helper.unpack()

        # Store fission yield dictionaries
        fission_yields = []

        # Create arrays to store fission Q values, reaction rates, and nuclide
        # numbers, zeroed out in material iteration
        number = np.empty(rates.n_nuc * mp) #FETs 

        fission_ind = rates.index_rx["fission"]

        # Extract results
        for i, mat in enumerate(self.local_mats):
            # Get tally index
            mat_index = self._mat_index_map[mat]

            # Zero out reaction rates and nuclide numbers
            number.fill(0.0)

            # Get new number densities
            for nuc, i_nuc_results in zip(nuclides, nuc_ind):
                #FETs
                i_nuc = self.number.index_nuc[nuc]
                #print(nuclides, nuc_ind, i_nuc)
                for j in range(mp):
                    number[i_nuc_results * mp + j] = self.number[mat, i_nuc * mp + j] 
                # number[i_nuc_results] = self.number[mat, nuc] #FETs 

            tally_rates = self._rate_helper.get_material_rates(
                mat_index, nuc_ind, react_ind, fet_deplete=fet_deplete) #FETs 

            # Compute fission yields for this material
            fission_yields.append(self._yield_helper.weighted_yields(i)) #FETs no change

            # Accumulate energy from fission
            self._energy_helper.update(tally_rates[:, fission_ind * mp], mat_index) #FETs updated index

            # Divide by total number and store
            rates[i] = self._rate_helper.divide_by_adens(number, fet_deplete=fet_deplete) #FETs 

        # Reduce energy produced from all processes
        # J / s / source neutron
        energy = comm.allreduce(self._energy_helper.energy)
        print(energy, "in operator.py")

        # Guard against divide by zero
        if energy == 0:
            if comm.rank == 0:
                sys.stderr.flush()
                print(" No energy reported from OpenMC tallies. Do your HDF5 "
                      "files have heating data?\n", file=sys.stderr, flush=True)
            comm.barrier()
            comm.Abort(1)

        # Scale reaction rates to obtain units of reactions/sec
        rates *= power / energy
        
        # Jiankai testing
        # for i, mat in enumerate(self.local_mats):
        #     mat_index = self._mat_index_map[mat]
        #     for nuc, i_nuc_results in zip(nuclides, nuc_ind):
        #         i_nuc = self.number.index_nuc[nuc] 
        #         print(mat, nuc, rates[i][i_nuc][:])
        # end
        
        # Store new fission yields on the chain
        self.chain.fission_yields = fission_yields
        
        return OperatorResult(k_combined, rates)
        

    def _get_nuclides_with_data(self):
        """Loads a cross_sections.xml file to find participating nuclides.

        This allows for nuclides that are important in the decay chain but not
        important neutronically, or have no cross section data.
        """

        # Reads cross_sections.xml to create a dictionary containing
        # participating (burning and not just decaying) nuclides.

        try:
            filename = os.environ["OPENMC_CROSS_SECTIONS"]
        except KeyError:
            filename = None

        nuclides = set()

        try:
            tree = ET.parse(filename)
        except Exception:
            if filename is None:
                msg = "No cross_sections.xml specified in materials."
            else:
                msg = 'Cross section file "{}" is invalid.'.format(filename)
            raise IOError(msg)

        root = tree.getroot()
        for nuclide_node in root.findall('library'):
            mats = nuclide_node.get('materials')
            if not mats:
                continue
            for name in mats.split():
                # Make a burn list of the union of nuclides in cross_sections.xml
                # and nuclides in depletion chain.
                if name not in nuclides:
                    nuclides.add(name)

        return nuclides

    def get_results_info(self):
        """Returns volume list, material lists, and nuc lists.

        Returns
        -------
        volume : dict of str float
            Volumes corresponding to materials in full_burn_dict
        nuc_list : list of str
            A list of all nuclide names. Used for sorting the simulation.
        burn_list : list of int
            A list of all material IDs to be burned.  Used for sorting the simulation.
        full_burn_list : list
            List of all burnable material IDs

        """
        nuc_list = self.number.burnable_nuclides
        burn_list = self.local_mats

        volume = {}
        for i, mat in enumerate(burn_list):
            volume[mat] = self.number.volume[i]

        # Combine volume dictionaries across processes
        volume_list = comm.allgather(volume)
        volume = {k: v for d in volume_list for k, v in d.items()}

        return volume, nuc_list, burn_list, self.burnable_mats
