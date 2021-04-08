import h5py
import numpy as np

from .results import Results, _VERSION_RESULTS
from openmc.checkvalue import check_filetype_version


__all__ = ["ResultsList"]


class ResultsList(list):
    """A list of openmc.deplete.Results objects

    It is recommended to use :meth:`from_hdf5` over
    direct creation.
    """

    @classmethod
    def from_hdf5(cls, filename):
        """Load in depletion results from a previous file

        Parameters
        ----------
        filename : str
            Path to depletion result file

        Returns
        -------
        new : ResultsList
            New instance of depletion results
        """
        with h5py.File(str(filename), "r") as fh:
            check_filetype_version(fh, 'depletion results', _VERSION_RESULTS[0])
            new = cls()

            # Get number of results stored
            n = fh["number"][...].shape[0]

            for i in range(n):
                new.append(Results.from_hdf5(fh, i))
        return new

    def get_atoms(self, mat, nuc, fet_deplete=None):
        """Get number of nuclides over time from a single material

        .. note::

            Initial values for some isotopes that do not appear in
            initial concentrations may be non-zero, depending on the
            value of :class:`openmc.deplete.Operator` ``dilute_initial``.
            The :class:`openmc.deplete.Operator` adds isotopes according
            to this setting, which can be set to zero.

        Parameters
        ----------
        mat : str
            Material name to evaluate
        nuc : str
            Nuclide name to evaluate

        Returns
        -------
        time : numpy.ndarray
            Array of times in [s]
        concentration : numpy.ndarray
            Total number of atoms for specified nuclide

        """
        #FETs 
        mp = 1
        if fet_deplete is not None:
            if fet_deplete['name']== 'zernike':
                mp = zer.num_poly(fet_deplete['order'])
            elif fet['name']=='zernike1d':
                mp = zer.num_poly1d(fet_deplete['order']) 
        #
        
        time = np.empty_like(self, dtype=float)
        concentration = np.empty_like(self, dtype=float)

        # Evaluate value in each region
        for i, result in enumerate(self):
            time[i] = result.time[0]
            #FETs
            if fet_deplete is not None:
                for j in range(mp):
                    concentration[i * mp + j] = result[0, mat, nuc * mp + j] 
            else:
                concentration[i] = result[0, mat, nuc] 

        return time, concentration

    def get_reaction_rate(self, mat, nuc, rx, fet_deplete=None):
        """Get reaction rate in a single material/nuclide over time

        .. note::

            Initial values for some isotopes that do not appear in
            initial concentrations may be non-zero, depending on the
            value of :class:`openmc.deplete.Operator` ``dilute_initial``
            The :class:`openmc.deplete.Operator` adds isotopes according
            to this setting, which can be set to zero.

        Parameters
        ----------
        mat : str
            Material name to evaluate
        nuc : str
            Nuclide name to evaluate
        rx : str
            Reaction rate to evaluate

        Returns
        -------
        time : numpy.ndarray
            Array of times in [s]
        rate : numpy.ndarray
            Array of reaction rates

        """
        mp = 1
        if fet_deplete is not None:
            if fet_deplete['name']== 'zernike':
                mp = zer.num_poly(fet_deplete['order'])
            elif fet['name']=='zernike1d':
                mp = zer.num_poly1d(fet_deplete['order']) 
        #
        time = np.empty_like(self, dtype=float)
        rate = np.empty_like(self, dtype=float)

        # Evaluate value in each region
        if fet_deplete is None:
            for i, result in enumerate(self):
                time[i] = result.time[0]
                rate[i] = result.rates[0].get(mat, nuc, rx) * result[0, mat, nuc]
        else:
            for i, result in enumerate(self):
                time[i] = result.time[0]
                for j in range(mp):
                    rate[i * mp + j] = result.rates[0].get(mat, nuc, rx * mp + j) \
                    * result[0, mat, nuc * mp] #FETs To be determined: result[0, mat, nuc * mp + j] 
        return time, rate

    def get_eigenvalue(self):
        """Evaluates the eigenvalue from a results list.

        Returns
        -------
        time : numpy.ndarray
            Array of times in [s]
        eigenvalue : numpy.ndarray
            k-eigenvalue at each time. Column 0
            contains the eigenvalue, while column
            1 contains the associated uncertainty

        """
        time = np.empty_like(self, dtype=float)
        eigenvalue = np.empty((len(self), 2), dtype=float)

        # Get time/eigenvalue at each point
        for i, result in enumerate(self):
            time[i] = result.time[0]
            eigenvalue[i] = result.k[0]

        return time, eigenvalue

    def get_depletion_time(self):
        """Return an array of the average time to deplete a material

        .. note::

            Will have one fewer row than number of other methods,
            like :meth:`get_eigenvalues`, because no depletion
            is performed at the final transport stage

        Returns
        -------
        times : :class:`numpy.ndarray`
            Vector of average time to deplete a single material
            across all processes and materials.

        """
        times = np.empty(len(self) - 1)
        # Need special logic because the predictor
        # writes EOS values for step i as BOS values
        # for step i+1
        # The first proc_time may be zero
        if self[0].proc_time > 0.0:
            items = self[:-1]
        else:
            items = self[1:]
        for ix, res in enumerate(items):
            times[ix] = res.proc_time
        return times
    
     def export_to_materials_xml(self, mats_list, burnup_index, nuc_with_data=None, fet_deplete=None):
        """
        """
        result = self[burnup_index]
        mat_file = openmc.Materials()
        for i, mat in enumerate(mats_list):
            id = mat.id
            new_mat = openmc.Material(id)
            mat_id = str(mat.id)
            if mat_id in result.mat_to_ind.keys():
                new_mat.volume = result.volume[str(mat_id)]
                for nuc in result.nuc_to_ind.keys():
                    atoms = result[0, str(mat_id), nuc]
                    if atoms > 0.0:
                        atoms_per_barn_cm = 1e-24 * atoms / new_mat.volume
                        if nuc_with_data is None:
                            if fet_deplete is None:
                                new_mat.add_nuclide(nuc, atoms_per_barn_cm)
                            else: 
                                new_mat.add_nuclide_fet(nuc, atoms_per_barn_cm, fet_deplete=fet_deplete)
                        elif nuc in nuc_with_data:
                            if fet_deplete is None:
                                new_mat.add_nuclide(nuc, atoms_per_barn_cm) 
                            else:
                                new_mat.add_nuclide_fet(nuc, atoms_per_barn_cm, fet_deplete=fet_deplet)
                mat_file.append(new_mat)
            else:
                mat_file.append(mat)
        mat_file.export_to_xml()
        
        