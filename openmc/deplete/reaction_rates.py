"""ReactionRates module.

An ndarray to store reaction rates with string, integer, or slice indexing.
"""
from collections import OrderedDict

import numpy as np

from openmc import zernike as zer

__all__ = ["ReactionRates"]


class ReactionRates(np.ndarray):
    """Reaction rates resulting from a transport operator call

    This class is a subclass of :class:`numpy.ndarray` with a few custom
    attributes that make it easy to determine what index corresponds to a given
    material, nuclide, and reaction rate.

    Parameters
    ----------
    local_mats : list of str
        Material IDs
    nuclides : list of str
        Depletable nuclides
    reactions : list of str
        Transmutation reactions being tracked
    from_results : boolean
        If the reaction rates are loaded from results, indexing dictionnaries
        need to be kept the same.

    Attributes
    ----------
    index_mat : OrderedDict of str to int
        A dictionary mapping material ID as string to index.
    index_nuc : OrderedDict of str to int
        A dictionary mapping nuclide name as string to index.
    index_rx : OrderedDict of str to int
        A dictionary mapping reaction name as string to index.
    n_mat : int
        Number of materials.
    n_nuc : int
        Number of nucs.
    n_react : int
        Number of reactions.

    """

    # NumPy arrays can be created 1) explicitly 2) using view casting, and 3) by
    # slicing an existing array. Because of these possibilities, it's necessary
    # to put initialization logic in __new__ rather than __init__. Additionally,
    # subclasses need to handle the multiple ways of creating arrays by using
    # the __array_finalize__ method (discussed here:
    # https://docs.scipy.org/doc/numpy/user/basics.subclassing.html)

    def __new__(cls, local_mats, nuclides, reactions, from_results=False, max_poly_order=None):
        # Create appropriately-sized zeroed-out ndarray
        n_poly = 1
        if max_poly_order is not None:
            n_poly = zer.num_poly(max_poly_order)
            shape = (len(local_mats), len(nuclides), len(reactions), n_poly)  # FETs
        else:
            shape = (len(local_mats), len(nuclides), len(reactions))
        obj = super().__new__(cls, shape)
        obj[:] = 0.0

        # Add mapping attributes, keep same indexing if from depletion_results
        if from_results:
            obj.index_mat = local_mats
            obj.index_nuc = nuclides
            obj.index_rx = reactions
            obj.n_poly = n_poly # FETs
        # Else, assumes that reaction rates are ordered the same way as
        # the lists of local_mats, nuclides and reactions (or keys if these
        # are dictionnaries)
        else:
            obj.index_mat = {mat: i for i, mat in enumerate(local_mats)}
            obj.index_nuc = {nuc: i for i, nuc in enumerate(nuclides)}
            obj.index_rx = {rx: i for i, rx in enumerate(reactions)}
            obj.n_poly = n_poly # FETs 
            
        return obj

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self.index_mat = getattr(obj, 'index_mat', None)
        self.index_nuc = getattr(obj, 'index_nuc', None)
        self.index_rx = getattr(obj, 'index_rx', None)
        self.n_poly = getattr(obj, 'n_poly', 1) # FETs 
        
    # Reaction rates are distributed to other processes via multiprocessing,
    # which entails pickling the objects. In order to preserve the custom
    # attributes, we have to modify how the ndarray is pickled as described
    # here: https://stackoverflow.com/a/26599346/1572453

    def __reduce__(self):
        state = super().__reduce__()
        new_state = state[2] + (self.index_mat, self.index_nuc, self.index_rx, self.n_poly) # FETs 
        return (state[0], state[1], new_state)

    def __setstate__(self, state):
        self.index_mat = state[-4] # FETs
        self.index_nuc = state[-3]
        self.index_rx = state[-2]
        self.n_poly = state[-1]
        super().__setstate__(state[0:-4])

    @property
    def n_mat(self):
        return len(self.index_mat)

    @property
    def n_nuc(self):
        return len(self.index_nuc)

    @property
    def n_react(self):
        return len(self.index_rx)

    def get(self, mat, nuc, rx):
        """Get reaction rate by material/nuclide/reaction

        Parameters
        ----------
        mat : str
            Material ID as a string
        nuc : str
            Nuclide name
        rx : str
            Name of the reaction

        Returns
        -------
        float array
            Reaction rate corresponding to given material, nuclide, and reaction

        """
        mat = self.index_mat[mat]
        nuc = self.index_nuc[nuc]
        rx = self.index_rx[rx]
        return self[mat, nuc, rx, :]  # FETs 

    def set(self, mat, nuc, rx, fet):
        """Set reaction rate by material/nuclide/reaction

        Parameters
        ----------
        mat : str
            Material ID as a string
        nuc : str
            Nuclide name
        rx : str
            Name of the reaction
        fet : array
            A ZernikePolynomial object

        """
        mat = self.index_mat[mat]
        nuc = self.index_nuc[nuc]
        rx = self.index_rx[rx]
        self[mat, nuc, rx, :] = fet[:] # FETs 
    
    
    