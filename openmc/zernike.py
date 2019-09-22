"""Zernike module.

Utility class for zernike polynomials
"""

import numpy as np
from math import pi
import numpy as np
import copy
from scipy.optimize import minimize
from openmc.zernike_data import b_matrix, c_matrix

def num_poly(n):
    return int(1/2 * (n+1) * (n+2))

def zern_to_ind(n,m):
    ind = num_poly(n-1)
    ind += (m + n) // 2
    return int(ind)

def form_b_matrix(p, pp, rate):
    # Yields the sum
    # ret = sum_r Int[P_p * P_pp * P_r * rate_r] / Int[P_pp^2]

    order = len(rate)

    v1 = b_matrix[p,pp,0:order]
    return np.dot(v1, rate)/c_matrix[pp]

def flat_to_zern(val, order, radius):
    # Get number of unknowns
    n = num_poly(order)
    coeffs = np.zeros(n)
    coeffs[0] = val

    return ZernikePolynomial(order, coeffs, radial_norm=radius, sqrt_normed=False)

class ZernikePolynomial:
    ''' ZernikePolynomial class
    
    This class contains the data that fully describes a Zernike polynomial
    and supporting functions.

    Attributes
    ----------
    order : int
        The maximum order of the polynomial
    radial_norm : flaot
        The radial normalization factor
    name : str
        The name of the polynomial
    n_coeffs : int
        The number of coefficients for a given order
    coeffs : List[flaot]
        The coefficients of the polynomial ordered as (m,n) 
        (0,0) (-1,1) (1,1) (-2,2) (0,2) (2,2) ...
    p_coeffs : List[float]
        Precomputed polynomial coefficients so that factorials
        do not need to be evaluated for funciton evaluation.
    '''
    
    import openmc.zernike_int as zni
    
    def __init__(self, order, coeffs, radial_norm=1.0, sqrt_normed=False):
        self._order = order
        self._radial_norm = radial_norm
        self._name = ''
        self._coeffs = coeffs
        self._n_coeffs = int(1/2 * (order+1) * (order+2))
        self._sqrt_normed = sqrt_normed
        self._p_coeffs = self.precompute_zn_coeffs()

    def __mul__(self, other):
        new = copy.deepcopy(self)
        new.coeffs *= other
        return new

    def __rmul__(self, other):
        new = copy.deepcopy(self)
        new.coeffs *= other
        return new

    def __div__(self, other):
        new = copy.deepcopy(self)
        new.coeffs /= other
        return new
    
    def __truediv__(self, other):
        new = copy.deepcopy(self)
        new.coeffs /= other
        return new
        
    @property
    def order(self):
        return self._order

    @property
    def radial_norm(self):
        return self._radial_norm
    
    @property
    def coeffs(self):
        return self._coeffs

    @property
    def n_coeffs(self):
        return self._n_coeffs
    
    @property
    def name(self):
        return self._name

    @property
    def sqrt_normed(self):
        return self._sqrt_normed

    @order.setter
    def order(self, order):
        self._order = order

    @coeffs.setter
    def coeffs(self, coeffs):
        self._coeffs = coeffs

    @n_coeffs.setter
    def n_coeffs(self, n_coeffs):
        self._n_coeffs = n_coeffs
            
    @radial_norm.setter
    def radial_norm(self, radial_norm):
        self._radial_norm = radial_norm
        
    @name.setter
    def name(self, name):
        self._name = name

    @sqrt_normed.setter
    def sqrt_normed(self, sqrt_normed):
        self._sqrt_normed = sqrt_normed


    def order_to_index(self, n, m):
        ''' Returns the index for accessing coefficients based 
            on the (n,m) values.  
        '''
        return int(1/2 * (n) * (n+1))  + (m + n) // 2

    def precompute_zn_coeffs(self):
        ''' Precompute the zernike polynomial leading coefficeints

            Note that all FETs in OpenMC and reconstruction in
            MOOSE assume that the square root of the normalization
            constant is included.  These are embedded in poly_coeffs.
        '''

        poly_coeffs = []

        for n in range(0,(self.order+1)):
            for m in range(-n,(n+1),2):
                for s in range(0,(n-abs(m))//2+1):
                    poly_coeffs.append(self.R_m_n_s(n,m,s,self.radial_norm))

        return poly_coeffs

    def get_poly_value(self, r, theta):
        ''' Compute the value of a polynomial at a point.

            Note that the precomputed leading coefficients,
            p_coeffs[], assuems that the square root of the
            normalization constant is included in the coefficients
            coeffs[]

            Parameters
            ----------
            r : float
                 The radial point.  Not normalizated.
            theta : float
                 The theta value in radians.

        '''
        import math

        r = r / self.radial_norm

        # Determine the vector of sin(n*theta) and cos(n*theta)
        sin_phi = math.sin(theta)
        cos_phi = math.cos(theta)

        sin_phi_vec = np.empty(self.order+1)
        cos_phi_vec = np.empty(self.order+1)

        sin_phi_vec[0] = 1.0
        cos_phi_vec[0] = 1.0

        sin_phi_vec[1] = 2.0 * cos_phi
        cos_phi_vec[1] = cos_phi

        for i in range(2,self.order+1):
            sin_phi_vec[i] = 2.0 * cos_phi * sin_phi_vec[i-1] - sin_phi_vec[i-2]
            cos_phi_vec[i] = 2.0 * cos_phi * cos_phi_vec[i-1] - cos_phi_vec[i-2]

        sin_phi_vec = sin_phi_vec * sin_phi

        # Calculate R_m_n(rho)

        zn_mat = np.empty([self.order+1, self.order+1])

        # Fill out the main diagonal first
        for p in range(0,self.order+1):
            zn_mat[p][p] = r**p
        # Fill in the second diagonal
        for q in range(0,(self.order-2+1)):
            zn_mat[q+2][q] = (q+2) * zn_mat[q+2][q+2] - (q+1) * zn_mat[q][q]
        # Fill in the rest of the values using the original results
        for p in range(4,(self.order+1)):
            k2 = float(2 * p * (p - 1) * (p - 2))
            for q in range(p-4,-1,-2):
                k1 = float((p + q) * (p - q) * (p - 2.0) / 2.0)
                k3 = float(-q**2*(p - 1) - p * (p - 1) * (p - 2))
                k4 = float(-p * (p + q - 2) * (p - q - 2) / 2.0)
                zn_mat[p][q] = ((k2 * r**2 + k3) * zn_mat[p-2][q] + k4 * zn_mat[p-4][q]) / k1

        val = 0.0
        i = 0
        for p in range(0,(self.order+1)):
            for q in range(-p,p+1,2):
                norm = self.get_norm_factor(p,q)
                if(q < 0):
                    val += zn_mat[p][abs(q)] * sin_phi_vec[p-1] * norm * self.coeffs[i]
                elif (q == 0):
                    val += zn_mat[p][q] * norm * self.coeffs[i]
                else:
                    val += zn_mat[p][q] * cos_phi_vec[p] * norm * self.coeffs[i]
                i = i + 1

        return val

    def get_poly_value_quick(self, r, theta):
        ''' Compute the value of a polynomial at a point.

            Parameters
            ----------
            r : float
                 The radial point.  Not normalizated.
            theta : float
                 The theta value in radians.

        '''
        import math

        r = r / self.radial_norm
        
        # Determine the vector of sin(n*theta) and cos(n*theta)
        sin_phi = math.sin(theta)
        cos_phi = math.cos(theta)

        sin_phi_vec = np.empty(self.order+1)
        cos_phi_vec = np.empty(self.order+1)

        sin_phi_vec[0] = 1.0
        cos_phi_vec[0] = 1.0

        sin_phi_vec[1] = 2.0 * cos_phi
        cos_phi_vec[1] = cos_phi

        for i in range(2,self.order+1):
            sin_phi_vec[i] = 2.0 * cos_phi * sin_phi_vec[i-1] - sin_phi_vec[i-2]
            cos_phi_vec[i] = 2.0 * cos_phi * cos_phi_vec[i-1] - cos_phi_vec[i-2]

        sin_phi_vec = sin_phi_vec * sin_phi
        
        # Calculate R_m_n(rho)

        zn_mat = np.empty([self.order+1, self.order+1])
        
        # Fill out the main diagonal first
        for p in range(0,self.order+1):
            zn_mat[p][p] = r**p
        # Fill in the second diagonal
        for q in range(0,(self.order-2+1)):
            zn_mat[q+2][q] = (q+2) * zn_mat[q+2][q+2] - (q+1) * zn_mat[q][q]
        # Fill in the rest of the values using the original results
        for p in range(4,(self.order+1)):
            k2 = float(2 * p * (p - 1) * (p - 2))
            for q in range(p-4,-1,-2):
                k1 = float((p + q) * (p - q) * (p - 2.0) / 2.0)
                k3 = float(-q**2*(p - 1) - p * (p - 1) * (p - 2))
                k4 = float(-p * (p + q - 2) * (p - q - 2) / 2.0)
                zn_mat[p][q] = ((k2 * r**2 + k3) * zn_mat[p-2][q] + k4 * zn_mat[p-4][q]) / k1

        val = 0.0
        i = 0
        for p in range(0,(self.order+1)):
            for q in range(-p,p+1,2):
                norm = self.get_norm_factor(p,q)
                # Get norm factor has a 1/pi term for the norm if self.sqrt_normed == true
                # We don't want that 1/pi term -- only the sqrt() term
                if (self.sqrt_normed):
                    norm = norm * math.pi
                if(q < 0):
                    val += zn_mat[p][abs(q)] * sin_phi_vec[p-1] * norm * self.coeffs[i]
                elif (q == 0):
                    val += zn_mat[p][q] * norm * self.coeffs[i]
                else:
                    val += zn_mat[p][q] * cos_phi_vec[p] * norm * self.coeffs[i]
                i = i + 1

        return val

    def get_poly_value_quick_vec(self, r, theta):
        ''' Compute the value of a polynomial at a point.

            Parameters
            ----------
            r : np.array
                 The radial points.  Not normalizated.
            theta : np.array
                 The theta values in radians.

        '''
        import math
        import sys
        from copy import deepcopy

        if(isinstance(r,float)):
            r_temp = np.empty([1])
            r_temp[0] = r
            r = deepcopy(r_temp)
            theta_temp = np.empty([1])
            theta_temp[0] = theta
            theta = deepcopy(theta_temp)

        r = r / self.radial_norm

        n_points = r.size * theta.size

        # Determine the vector of sin(theta) and cos(theta)
        sin_phi = np.sin(theta)
        cos_phi = np.cos(theta)

        # Determine the matrix of sin(n*theta) and cos(n*theta)
        sin_phi_vec = np.ones((self.order+1,theta.size))
        cos_phi_vec = np.ones((self.order+1,theta.size))

        sin_phi_vec[1] = 2.0 * cos_phi
        cos_phi_vec[1] = cos_phi

        for i in range(2,self.order+1):
            sin_phi_vec[i] = 2.0 * cos_phi * sin_phi_vec[i-1] - sin_phi_vec[i-2]
            cos_phi_vec[i] = 2.0 * cos_phi * cos_phi_vec[i-1] - cos_phi_vec[i-2]

        sin_phi_vec = sin_phi_vec * sin_phi

        # Calculate R_m_n(rho)
        zn_mat = np.zeros([self.order+1, self.order+1,r.size])

        # Fill out the main diagonal first
        for p in range(0,self.order+1):
            p_pow = r**p
            for ii in range(0,p_pow.size):
                zn_mat[p][p] = p_pow

        # Fill in the second diagonal
        for q in range(0,(self.order-2+1)):
            zn_mat[q+2][q] = (q+2) * zn_mat[q+2][q+2] - (q+1) * zn_mat[q][q]

        # Fill in the rest of the values using the original results
        for p in range(4,(self.order+1)):
            k2 = float(2 * p * (p - 1) * (p - 2))
            for q in range(p-4,-1,-2):
                k1 = float((p + q) * (p - q) * (p - 2.0) / 2.0)
                k3 = float(-q**2*(p - 1) - p * (p - 1) * (p - 2))
                k4 = float(-p * (p + q - 2) * (p - q - 2) / 2.0)
                zn_mat[p][q] = ((k2 * r**2 + k3) * zn_mat[p-2][q] + k4 * zn_mat[p-4][q]) / k1

        val = np.zeros(r.size)
        i = 0
        for p in range(0,(self.order+1)):
            for q in range(-p,p+1,2):
                norm = self.get_norm_factor(p,q)
                # Get norm factor has a 1/pi term for the norm if self.sqrt_normed == true
                # We don't want that 1/pi term -- only the sqrt() term
                if (self.sqrt_normed):
                    norm = norm * math.pi
                if(q < 0):
                    val += zn_mat[p][abs(q)] * sin_phi_vec[p-1] * norm * self.coeffs[i]
                elif (q == 0):
                    val += zn_mat[p][q] * norm * self.coeffs[i]
                else:
                    val += zn_mat[p][q] * cos_phi_vec[p] * norm * self.coeffs[i]
                i = i + 1

        return val
    
    def compute_integral(self, r_min, r_max, theta_min, theta_max):
        ''' Compute the integral of the zernike polynomial over some 
        subset of the unit disk

        Note that the normalization factor is included because
        it is assumed that the polynomials coefficients include
        the square root of the normalization constants but
        integrate_wedge() does not use this multiplicate constant.

        Parameters
        ----------
        r_min : float
            The inner radius
        r_max : float
            The outer radius
        theta_min : float
            The minimum theta value
        theta_max : float
            The maxmium theta value
        '''

        import math
        import openmc.zernike_int as zni
        import sys

        val = 0.0

        for n in range(0, self.order+1):
            for m in range(-n,(n+1),2):
                norm_factor = self.get_norm_factor(n,m)
                val += self.coeffs[self.order_to_index(n,m)] * norm_factor * \
                       zni.integrate_wedge(n,m,r_min,r_max,theta_min, theta_max)

        return val

    def plot_disk(self, n_rings, n_sectors, fname):
        ''' This function plots the volume averaged value of the
        Zernike polynonmials in radial rings and azimuthal
        sectors.

        Parameters
        ----------
        n_rings : int
             The number of rings to split the disk into.
        n_sectors : int
             The number of azimuthal sectors in each ring.
        fname : str
             The name of the file into which to save the plot.
        
        Returns
        -------
        patch_vals : List
            List of patch values used in plot.
        '''

        import math
        import numpy
        import matplotlib.pyplot as plt
        from matplotlib.patches import Circle, Wedge, Polygon
        from matplotlib.collections import PatchCollection
        
        vol_per_ring = self.radial_norm * self.radial_norm * math.pi / n_rings       
        ring_radii = [0.0]
        for i in range(0, n_rings):
            ring_radii.append( math.sqrt(vol_per_ring / math.pi \
                            + ring_radii[-1] * ring_radii[-1]) )

        theta_cuts = [0.0]
        for i in range(0, n_sectors):
            theta_cuts.append(theta_cuts[-1] + \
                              math.pi * 2.0 / n_sectors)

        patches = []
        patch_vals = []
        for i in range(0, n_rings):
            for j in range(0, n_sectors):
                
                thickness = ring_radii[i+1] - ring_radii[i]
                if(j == n_sectors-1):  # Loop around the sectors
                    patch_vals.append(self.compute_integral(ring_radii[i]/self.radial_norm, ring_radii[i+1]/self.radial_norm,\
                                                            theta_cuts[j], 2*math.pi))
                    wedge = Wedge( (0.0,0.0), ring_radii[i+1], theta_cuts[j]*180.0/math.pi, theta_cuts[0]*180.0/math.pi, width=thickness, linewidth=None)
                else:
                    patch_vals.append(self.compute_integral(ring_radii[i]/self.radial_norm, ring_radii[i+1]/self.radial_norm,\
                                                            theta_cuts[j], theta_cuts[j+1]))
                    wedge = Wedge( (0.0,0.0), ring_radii[i+1], theta_cuts[j]*180.0/math.pi, theta_cuts[j+1]*180.0/math.pi, width=thickness, linewidth=None)
                patches.append(wedge)
                
        fig, ax = plt.subplots()
        ax.set_xlim([-self.radial_norm,self.radial_norm])
        ax.set_ylim([-self.radial_norm,self.radial_norm])
        p = PatchCollection(patches, cmap=plt.cm.jet, linewidths=0.0)
        p.set_array(numpy.array(patch_vals))
        ax.add_collection(p)
        plt.colorbar(p)
        fig.savefig(fname)
        plt.close()

        return patch_vals

    def plot_over_line(self, theta, fname):
        ''' This funciton plots the polynomial over the entire disk
            along a specific theta value

        Parameters
        ----------
        theta : float
             The theta value over which to plot from [0, R_max]
        fname : str
             The filename to use when saving the figure
        
        '''

        import matplotlib.pyplot as plt

        # We will use a linearly spaced set of points for now to plot
        r_vals = np.linspace(0.0, self.radial_norm, num=1000)
        vals = []

        for r in r_vals:
            vals.append(self.get_poly_value(r,theta))

        fig = plt.figure()
        ax = fig.add_subplot(111)
        ax.plot(r_vals, vals)
        fig.savefig(fname)
        plt.close()

    def remove_fet_sqrt_normalization(self):
        ''' This function removes the sqrt(2(n+1)) or sqrt(n+1)
        normalization that is applied in OpenMC FETs.
        '''
        import math

        for n in range(0, self.order+1):
            for m in range(-n,(n+1),2):
                if (m == 0):
                    self.coeffs[self.order_to_index(n,m)] *= math.sqrt(n+1.0)
                else:
                    self.coeffs[self.order_to_index(n,m)] *= math.sqrt(2.0*n+2.0)

        self.sqrt_normed = False
        # Since we might have changed the normalization state, we need
        # to recompute the precomputed polynomial coefficients
        self._p_coeffs = self.precompute_zn_coeffs()

    def apply_fet_sqrt_normalization(self):
        ''' This function applies the sqrt(2(n+1)) or sqrt(n+1)
        normalization that is applied in OpenMC FETs.
        '''
        import math

        for n in range(0, self.order+1):
            for m in range(-n,(n+1),2):
                if (m == 0):
                    self.coeffs[self.order_to_index(n,m)] /= math.sqrt(n+1.0)
                else:
                    self.coeffs[self.order_to_index(n,m)] /= math.sqrt(2.0*n+2.0)

        self.sqrt_normed = True
        # Since we might have changed the normalization state, we need
        # to recompute the precomputed polynomial coefficients
        self._p_coeffs = self.precompute_zn_coeffs()

    def normalize_coefficients(self):
        ''' This function normalizes coefficients by (n+1) / pi or
        (2n+2) / pi.
        '''
        import math

        for n in range(0, self.order+1):
            for m in range(-n,(n+1),2):
                if (m == 0):
                    self.coeffs[self.order_to_index(n,m)] *= (n+1.0)
                else:
                    self.coeffs[self.order_to_index(n,m)] *= (2.0*n+2.0)

                self.coeffs[self.order_to_index(n,m)] /= math.pi

        self.sqrt_normed = False
        # Since we might have changed the normalization state, we need
        # to recompute the precomputed polynomial coefficients
        self._p_coeffs = self.precompute_zn_coeffs()

    def scale_coefficients(self, scale_value):
        ''' This function scales every coefficient in the expansion
        by a given value

        Parameters
        ----------
        scale_value : float
             The scaling value to apply
        '''

        for n in range(0, self.order+1):
            for m in range(-n,(n+1),2):
                self.coeffs[self.order_to_index(n,m)] *= scale_value


    def get_norm_factor(self, n, m):
        ''' This function determines the normalization factor
        to be applied

        Parameters
        ----------
        n : int
             The radial moment number
        m : int
             The azimuthal moment number
        '''
        import math
        
        if (m == 0 and self.sqrt_normed):
            return (1.0 / math.pi * math.sqrt(n + 1.0))
        elif (m != 0 and self.sqrt_normed):
            return (1.0 / math.pi * math.sqrt(2.0 * n + 2.0))
        elif (not self.sqrt_normed):
            return 1.0
        else:
            sys.exit("Invalid state when calculating normalization factor")

    def R_m_n_s(self, n, m, s, r_max=1.0):
        ''' This function calculates the R_{n,m}(r)
        coefficients.  Note that this funciton does not check if
        n, m, or s are valid.

        Parameters
        ----------
        n : int
             The azimuthal moment number
        m : int
             The azimuthal moment number
        s : int
             One of the summation terms that defines R_{n,m}
             s is defined as [0,(n-abs(m))/2]
        r_max : float
             The maximum radius of the disk

        '''

        import math

        norm_factor = self.get_norm_factor(n,m)

        return (1.0 / (r_max * r_max) * math.sqrt(norm_factor) * \
                math.pow(-1,s) * math.factorial(n-s) / \
                ( math.factorial(s) * math.factorial((n+m)//2 - s) * \
                  math.factorial( (n-m)//2 - s) ) )

    def openmc_form(self, radial_only=True):
        ''' This function returns the openmc poly_coeffs vector, which is
        radius followed by coefficients in barn/cm.

        Parameters
        ----------
        radial_only : bool
            Whether only Z_{n,0} terms are returned.
        
        '''
        import time
        
        start_time = time.time()
        self.force_positive()
        print('Force positive execution time  ' + str(time.time() - start_time))
        self.apply_fet_sqrt_normalization()

        if radial_only:
            n = int((self.order + 2)/2)
            form = np.zeros(n + 1)
            form[0] = self.radial_norm
            for i in range(n):
                form[i+1] = self.coeffs[zern_to_ind(2*i, 0)] / 1.0e24
        else:
            n = len(self.coeffs)
            form = np.zeros(n + 1)
            form[0] = self.radial_norm
            form[1::] = self.coeffs / 1.0e24
        
        self.remove_fet_sqrt_normalization()
        return form

    def product_integrate(self, other):
        ''' Integrates over the unit disk the convolution of two polynomials.

        If both polynomials are normalized (it is assumed that "other" is, as it
        is simply a tally from OpenMC), then the integral is the dot product of
        the vectors.

        Parameters
        ----------
        other : np.array(float)
            The second polynomial, normalized.
        
        '''

        self.apply_fet_sqrt_normalization()
        
        val = np.dot(self.coeffs, other)

        self.remove_fet_sqrt_normalization()
        return val

    def force_positive(self, N=200):
        ''' Computes the minimum of the function, and then shifts all but the
        first moment by a scaling parameter to guarantee positivity.
        '''

        fun = lambda x: self.get_poly_value_quick(x[0], x[1])

        r = np.linspace(0, self.radial_norm, N)
        theta = np.linspace(0, 2*np.pi, N)

        r, theta = np.meshgrid(r, theta)
        
        r = np.reshape(r, N**2)
        theta = np.reshape(theta, N**2)

        val = self.get_poly_value_quick_vec(r, theta)

        argmin = np.argmin(val)

        result = minimize(fun, [r.item(argmin), theta.item(argmin)], bounds=((0.0, self.radial_norm), (0.0, 2*np.pi)))

        y = result.fun

        #print(result) # modified by jiankai

        fudge_factor = 1.0e-6

        if y < fudge_factor * self.coeffs[0]:
            scaling = self.coeffs[0] / (self.coeffs[0] - y) / (1.0 + fudge_factor)
            print("scaling by = ", scaling)
            self.coeffs[1::] = self.coeffs[1::] * scaling

        
