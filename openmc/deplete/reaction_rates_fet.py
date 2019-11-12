"""ReactionRates module.

Just contains a dictionary of np.arrays to store reaction rates.
"""

import numpy as np
from openmc import zernike


class ReactionRatesFet:
    """ The Nuclide class.

    Contains everything in a depletion chain relating to a single nuclide.

    Parameters
    ----------
    mats_to_ind : OrderedDict[int]
        A dictionary mapping mats ID as string to index.
    nuc_to_ind : OrderedDict[int]
        A dictionary mapping nuclide name as string to index.
    react_to_ind : OrderedDict[int]
        A dictionary mapping reaction name as string to index.

    Attributes
    ----------
    mats_to_ind : OrderedDict[int]
        A dictionary mapping mats ID as string to index.
    nuc_to_ind : OrderedDict[int]
        A dictionary mapping nuclide name as string to index.
    react_to_ind : OrderedDict[int]
        A dictionary mapping reaction name as string to index.
    n_mats : int
        Number of mats.
    n_nuc : int
        Number of nucs.
    n_react : int
        Number of reactions.
    rates : np.array
        Array storing rates indexed by the above dictionaries.
    fets : np.array
        Array storing fet polynomial object indexed by the above dictionaries
    """

    def __init__(self, mats_to_ind, nuc_to_ind, react_to_ind, max_poly_order=10):

        self.mats_to_ind = mats_to_ind
        self.nuc_to_ind = nuc_to_ind
        self.react_to_ind = react_to_ind

        self.n_mats = len(mats_to_ind)
        self.n_nuc = len(nuc_to_ind)
        self.n_react = len(react_to_ind)
        self.max_poly_order = max_poly_order
        if max_poly_order != None:
            self.n_poly = zernike.num_poly(max_poly_order)
        else:
            self.n_poly = 1

        self.rates = np.zeros((self.n_mats, self.n_nuc, self.n_react, self.n_poly))
        self.fets = [[[None for x in range(self.n_react)] for y in range(self.n_nuc)] for z in range(self.n_mats)]
        #self.fets = np.zeros((self.n_mats, self.n_nuc, self.n_react))
    def __getitem__(self, pos):
        """ Retrieves an item from reaction_rates.

        Parameters
        ----------
        pos : Tuple
            A four-length tuple containing a mats index, a nuc index, a
            reaction index, and a polynomial index.  These indexes can be
            strings (which get converted to integers via the dictionaries),
            integers used directly, or slices.

        Returns
        -------
        np.array
            The value indexed from self.rates.
        """

        mats, nuc, react, poly = pos
        if isinstance(mats, str):
            mats_id = self.mats_to_ind[mats]
        else:
            mats_id = mats
        if isinstance(nuc, str):
            nuc_id = self.nuc_to_ind[nuc]
        else:
            nuc_id = nuc
        if isinstance(react, str):
            react_id = self.react_to_ind[react]
        else:
            react_id = react

        return self.rates[mats_id, nuc_id, react_id, poly]

    def __setitem__(self, pos, val):
        """ Sets an item from reaction_rates.

        Parameters
        ----------
        pos : Tuple
            A four-length tuple containing a mats index, a nuc index, a
            reaction index, and a polynomial index.  These indexes can be
            strings (which get converted to integers via the dictionaries),
            integers used directly, or slices.
        val : float
            The value to set the array to.
        """

        mats, nuc, react, poly = pos
        if isinstance(mats, str):
            mats_id = self.mats_to_ind[mats]
        else:
            mats_id = mats
        if isinstance(nuc, str):
            nuc_id = self.nuc_to_ind[nuc]
        else:
            nuc_id = nuc
        if isinstance(react, str):
            react_id = self.react_to_ind[react]
        else:
            react_id = react

        self.rates[mats_id, nuc_id, react_id, poly] = val
        
    # TODO make these properties with decorators 
    def set_fet(self, pos, fet):
        """ Sets the fet in the reaction rates object and tries
        to automatically set the reaction rate values as well
        
        Parameters
        ----------
        pos : Tuple
            A four-length tuple containing a mats index, a nuc index, a
            reaction index.  These indexes can be strings (which get 
            converted to integers via the dictionaries), integers used 
            directly, or slices.
        fet : A ZernikePolynomial object
        """
        mats, nuc, react = pos
        if isinstance(mats, str):
            mats_id = self.mats_to_ind[mats]
        else:
            mats_id = mats
        if isinstance(nuc, str):
            nuc_id = self.nuc_to_ind[nuc]
        else:
            nuc_id = nuc
        if isinstance(react, str):
            react_id = self.react_to_ind[react]
        else:
            react_id = react
            
        self.fets[mats_id][nuc_id][react_id] = fet

        for i in range(0, fet.n_coeffs):
            self[mats, nuc, react, i] = fet.coeffs[i]
            
    def get_fet(self, pos):
        """ Gets the fet in the reaction rates object
        
        Parameters
        ----------
        pos : Tuple
            A four-length tuple containing a mats index, a nuc index, a
            reaction index.  These indexes can be strings (which get 
            converted to integers via the dictionaries), integers used 
            directly, or slices.
        fet : A ZernikePolynomial object
        """
        mats, nuc, react = pos
        if isinstance(mats, str):
            mats_id = self.mats_to_ind[mats]
        else:
            mats_id = mats
        if isinstance(nuc, str):
            nuc_id = self.nuc_to_ind[nuc]
        else:
            nuc_id = nuc
        if isinstance(react, str):
            react_id = self.react_to_ind[react]
        else:
            react_id = react
        
        return self.fets[mats_id][nuc_id][react_id]
