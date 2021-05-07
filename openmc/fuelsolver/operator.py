"""
openmc.thsolver
==============

A thermal/hydraulic feedback front-end toola.
"""
from .materials import * 
from .variables import *  
#
class Operator(Variables):
    """
    """
    def __init__(self, power, conductance_gap, T_inlet, \
                       flow, zmeshes, radii, pitch, n_th_rings):
        """
        """
        super().__init__()
        self.th_setup(power, conductance_gap, T_inlet, \
                      flow, zmeshes, radii, pitch, n_th_rings)
    
    def thsolver(self, P_out, rel_power_dist):
        """
        """
        temp_dist = np.zeros((self.nrp5+1, self.nzth))
        den_dist = np.zeros(self.nzth)
        #
        for i in range(self.nzth):
            #self.relp[i] = np.sum(rel_power_dist[i])
            self.relp[i] = rel_power_dist[i]
        #
        self.plevel = 1.e+0
        tfmax = 0.e+0
        toutavg = 0.e+0
        dcoolavg = 0.e+0
        tcoolavg = 0.e+0
        tdoplmax = 0.e+0
        #
        fac = 1.0e+0/self.nzth
        #
        # coolant temperature calculation
        rhouin   = self.rhou[0]  # note of index from zero in fortran
        rhohuin  = self.rhohu[0] # note of index from zero in fortran
        #
        for kth in range(0, self.nzth):
           qprime    = self.plevel * self.powlin * self.relp[kth]
           qf        = self.fracdf * qprime / self.afp
           self.qvol[kth] = qf
           qc        = self.fracdc*qprime/self.acf
           qflux     = qf * self.afp / self.zeta
           qeffnew   = qflux*self.zetap+qc
           self.qeff[kth] = qeffnew
           #
           rhohuout  = rhohuin + qeffnew * self.hzth[kth]
           self.hcool[kth]= 0.5*(rhohuout + rhohuin)/rhouin
           self.tcool[kth]= ftemp(self.hcool[kth])
           self.dcool[kth]= fdensh(self.hcool[kth])
           rhohuin   = rhohuout
           self.rhohu[kth]= rhohuout
           dcoolavg  += self.dcool[kth]
           tcoolavg  += self.tcool[kth]
        #
        cool_out = [ftemp(rhohuout/rhouin) + self. CKELVIN, fdensh(rhohuout/rhouin)/1.0e+3] 
        #
        tcoolavg /= self.nzth
        dcoolavg /= self.nzth
        #
        hout = self.rhohu[self.nzth] / self.rhou[self.nzth]
        tout = ftemp(hout)
        #
        tfuelavg = 0
        is_converged = False
        #
        for k in range(0, self.nzth):
           qprime = self.plevel * self.powlin * self.relp[k]
           qf     = self.fracdf * qprime / self.afp
           #print(self.powlin)
           if(qf > 1e-11):
              # heat transfer coefficient
              # htcoef(k) = fhtcoef(tcool(k), self.deq, self.rhou(k)
              if(P_out > 15.513):
                  self.htcoef[k] = fhtcoef_if97(self.tcool[k], self.deq, self.rhou[k], P_out)
                  print('Warning: P_out > 15.513 MPa...')
              else:
                  #print(qf)
                  self.htcoef[k] = fhtcoef(self.tcool[k], self.deq, self.rhou[k])
              #
              is_converged = self.fuel_conduction(k, self.tcool[k], self.htcoef[k], qf)
           else:
              is_converged = False
           #
           if(not is_converged):
              return None
           #
        #
        tfuelavg *= fac
        #
        tmax = 0.0
        for k in range(self.nzth):
            temp_dist[0:self.nrp5,k] = self.tfuel[0:self.nrp5,k] + self.CKELVIN
            temp_dist[self.nrp5,k] =  self.tcool[k] + self.CKELVIN
            for i in range(0, self.nrp5):
               if(tmax <= temp_dist[i,k]):
                   tmax = temp_dist[i,k]
            den_dist[k] = self.dcool[k]/1.0e+3 # unit conversion 
        #
        #
        for k in range(0, self.nzth):
           if(den_dist[k] <= 0.0e+0 or tmax >= 3000.e+0):
              print("THSolver is converged, but the density is negative or fueal temperature is too high")
              is_converged = False
           #
        return is_converged, temp_dist, den_dist, cool_out

    def fuel_conduction(self, k, tcool1, htcoef1, qf):
        """
        """
        #
        is_converged = True
        #
        kf  =np.zeros(self.nrp5)
        kfm =np.zeros(self.nrp5)
        x   =np.zeros(self.nrp5)
        xd  =np.zeros(self.nrp5)
        ad  =np.zeros(self.nrp5)
        al  =np.zeros(self.nrp5)
        au  =np.zeros(self.nrp5)
        b   =np.zeros(self.nrp5)
        #
        errtf = 9.9e+30  # hunge() function in fortran 
        itr = 0
        #
        while (True):
            itr = itr + 1
            for i in range(0, self.nrp4):
                x[i] = self.tfuel[i,k] + self.CKELVIN
            #
            for i in range(0, self.nrp1):
                kf[i] = fkf(x[i])
            #
            for i in range(0, self.nr):
                kfm[i]=0.5*(kf[i]+kf[i + 1])
            #
            for i in range(self.nrp2-1, self.nrp4):
                kf[i] = fkc(x[i])
            #
            m = self.nrp3 - 1
            kmr = 0.5*(kf[m] + kf[m + 1])
            kml = 0.5*(kf[m] + kf[m - 1])
            # setup linear system
            for i in range(0, self.nrp4):
                x[i]  = self.tfuel[i,k]
                xd[i] = x[i]
            #
            qfd2 = qf * self.delr2 #*pp(1)
            m = 1 - 1 
            ad[m] = 4*kf[m]
            au[m] = -4*kf[m]
            b[m] = qf * self.radial_power_dist[m,k] # qfd2
            #
            i = m
            for m in range(1, self.nr):
                ri = 1/float(i+1) # fix a bug 
                ad[m] = kfm[i]+kfm[m] + 0.5*(kfm[m] - kfm[i])*ri
                al[m] = -kfm[i]*(1.0 - 0.5*ri)
                au[m] = -kfm[m]*(1.0 + 0.5*ri)
                b[m] = qf * self.radial_power_dist[m,k] #! qfd2
                i = m
            #
            #
            m = self.nrp1 - 1
            alpha = self.kgap*(1.0 - kf[m-1]/kf[m])
            ad[m] = 2.0*(kf[m]+self.kgap*(1.0 + 0.5/self.nr)) + alpha
            al[m] = -2.0*kf[m]
            au[m] = -2.0*self.kgap*(1.0 + 0.5/self.nr) - alpha
            b[m] = qf * self.radial_power_dist[m,k] # qfd2
            b[m] = qf * self.delr2
            #
            m = self.nrp2 - 1
            alpha = 2.0*self.kgap2*(kf[m + 1]/kf[m] - 1.0)
            ad[m] = 8.0*kf[m] + self.kgap4 - alpha
            al[m] = -self.kgap4 + alpha
            au[m] = -8.0*kf[m]
            b[m]  = 0.0
            m     = self.nrp3 - 1
            ad[m] = 4.0e+0*(kmr + kml) + self.tworm*(kmr - kml)
            al[m] = -kml*(4 - self.tworm)
            au[m] = -kmr*(4 + self.tworm)
            b[m]  = 0
            m = self.nrp4 - 1 
            kconv1= htcoef1 * self.tw
            alpha = 2.0e+0*kconv1*(1.0e+0 - kf[m-1]/kf[m])
            kconv = htcoef1 * self.tw * (4.0e+0 + self.tw/self.rw)
            ad[m] = 8.0e+0*kf[m] + kconv + alpha
            al[m] = -8.0e+0*kf[m]
            b[m]  = (kconv + alpha)*tcool1
            #
            # solve the tridiagonal system by gauss elimination
            im1 = 1 - 1
            for i in range(1, self.nrp4):
                aldi = al[i]/ad[im1]
                ad[i] = ad[i] - aldi * au[im1]
                b[i] = b[i] - aldi * b[im1]
                im1 = i
            #
            i    = self.nrp4 - 1
            ip1  = self.nrp4 - 1
            x[i] = b[i]/max(ad[i], 1.0e-30)
            for i in range(self.nrp3-1,-1,-1):
                x[i]=(b[i]-au[i]*x[ip1])/ad[i]
                ip1=i
            #
            errtf = 0.0e+0
            for i in range(0, self.nrp4):
                self.tfuel[i,k] = x[i]
                errtf = max(errtf, abs(x[i] - xd[i]))
            #
            if(errtf <= self.epstf):
                break
            #
            if(itr > 10000):
                is_converged = False
                break
            #
        #
        self.tfuel[self.nrp5-1,k] = ftfavg(x, self.nr, self.delr2)
        return is_converged
