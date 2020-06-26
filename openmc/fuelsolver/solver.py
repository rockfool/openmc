"""
openmc.thsolver
==============

A thermal/hydraulic feedback front-end tool.
"""
from .materials import * 


class THSolver(object):
    """
    """
    def th_solver_ss(P_out, rel_power_dist, temp_dist, den_dist, cool_out, is_converged):
        #use IF97PRO  ! for SCWR testing
        #real(8)     :: tfmax, toutavg, dcoolavg, tcoolavg, tdoplmax
        #real(8)     :: fac
        #integer     :: kth, k, i
        #real(8)     :: qprime
        #real(8)     :: qf
        #real(8)     :: rhohuin
        #real(8)     :: rhouin
        #real(8)     :: qc
        #real(8)     :: qeffnew
        #real(8)     :: qflux
        #real(8)     :: rhohuout
        #real(8)     :: hout
        #real(8)     :: tout
        #real(8)     :: tfuelavg
        #real(8)     :: rel_power_dist(:,:)
        #real(8)     :: temp_dist(:,:)
        #real(8)     :: den_dist(:)
        #real(8)     :: cool_out(:)
        #logical     :: is_converged
        #real(8)     :: tmax
        #real(8)     :: P_out
        #
        global nzth 
        global plevel  
        global tfmax   
        global toutavg 
        global dcoolavg
        global tcoolavg
        global tdoplmax
        #
        for i in range(0,nzth):
           relp[i] = sum(rel_power_dist[:][i])
        #
        plevel   = 1.e+0
        tfmax    = 0.e+0
        toutavg  = 0.e+0
        dcoolavg = 0.e+0
        tcoolavg = 0.e+0
        tdoplmax = 0.e+0
        #
        fac = 1.0e+0/nzth
        #
        # coolant temperature calculation
        rhouin   = rhou[0]
        rhohuin  = rhohu[0]
        #
        for kth in range(0, nzth)：
           qprime    = plevel*powlin*relp[kth]
           qf        = fracdf*qprime/afp
           qvol[kth] = qf
           qc        = fracdc*qprime/acf
           qflux     = qf*afp/zeta
           qeffnew   = qflux*zetap+qc
           qeff[kth] = qeffnew
           #
           rhohuout  = rhohuin+qeffnew*hzth[kth]
           hcool[kth]= 0.5*(rhohuout+rhohuin)/rhouin
           tcool[kth]= ftemp(hcool[kth])
           dcool[kth]= fdensh(hcool[kth])
           rhohuin   = rhohuout
           rhohu[kth]= rhohuout
           dcoolavg  = dcoolavg+dcool[kth]
           tcoolavg  = tcoolavg+tcool[kth]
        #
        cool_out[1-1] = ftemp(rhohuout/rhouin) + CKELVIN
        cool_out[2-1] = fdensh(rhohuout/rhouin)/1.0e+3
        #
        tcoolavg = tcoolavg/nzth
        dcoolavg = dcoolavg/nzth
        #
        hout = rhohu[nzth]/rhou[nzth]
        tout = ftemp[hout]
        #
        tfuelavg=0
        for k in range(0, nzth):
           qprime = plevel*powlin*relp[k]
           qf     = fracdf*qprime/afp
           if(qf > 1e-11)：
              # heat transfer coefficient
              # htcoef(k) = fhtcoef(tcool(k), deq, rhou(k)
              if(P_out > 15.513)： 
                 #htcoef[k] = fhtcoef_if97(tcool[k], deq, rhou[k],P_out)
                 print('Warning: P_out > 15.513 MPa...')
              else：
                 htcoef[k] = fhtcoef(tcool[k], deq, rhou[k])
              #
              fuel_t_dist(k, tcool[k], htcoef[k], qf, is_converged)
           else：
              is_converged = False
           #
           if(.not.is_converged)：
              return
           #
        #
        tfuelavg = tfuelavg*fac
        #
        tmax = 0.0
        for k in range(0, nzth):
           temp_dist(1:nrp5, k) =  tfuel(1:nrp5, k) + CKELVIN
           temp_dist(nrp5 + 1, k) =  tcool(k) + CKELVIN
           for i in range(0, nrp5):
              if(tmax <= temp_dist[i][k]):
                 tmax = temp_dist[i][k]
              #
           #
           den_dist[k] = dcool[k]/1.0e+3 # unit conversion 
        #
        #
        for k in range(0, nzth):
           if(den_dist[k] <= 0.0e+0 or tmax >= 3000.e+0):
              print("THSolver is converged, but the density is negative or fueal temperature is too high")
              is_converged = False
           #
        #

    def fuel_t_dist(k,tcool1,htcoef1,qf, is_converged):
        #integer     :: i, k, m, ip1
        #real(8)     :: tcool1, htcoef1, qf
        #real(8)     :: kmr         !  cladding conductivity eq. (5.5)
        #real(8)     :: kml         !  cladding conductivity eq. (5.5)
        #real(8)     :: qfd2        !
        #real(8)     :: ri          !
        #real(8)     :: alpha       !
        #real(8)     :: kconv1      !
        #real(8)     :: kconv       !
        #real(8)     :: aldi
        #integer     :: im1
        #logical     :: is_converged
        #real(8), allocatable          :: kf(:)   !thermal conductivity of uo2 at mesh pointss
        #real(8), allocatable          :: kfm(:)  ! kf at the middle points
        #real(8), allocatable          :: x(:)    ! fuel temp
        #real(8), allocatable          :: xd(:)
        #real(8), allocatable          :: ad(:)
        #real(8), allocatable          :: al(:)
        #real(8), allocatable          :: au(:)
        #real(8), allocatable          :: b(:)
        #real(8)                       :: errtf
        #integer                       :: itr
        #
        is_converged = True
        #
        kf  =[0.0]*nrp5
        kfm =[0.0]*nrp5
        x   =[0.0]*nrp5
        xd  =[0.0]*nrp5
        ad  =[0.0]*nrp5
        al  =[0.0]*nrp5
        au  =[0.0]*nrp5
        b   =[0.0]*nrp5
        #
        errtf = 9.9e+30  # hunge() function in fortran 
        itr = 0
        #
        while:
            itr = itr + 1
            for i in range(0, nrp4):
                x(i) = tfuel(i,k) + CKELVIN
            #
            for i in range(0, nrp1)
                kf(i) = fkf(x(i))
            #
            #
            for i in range(0, nr):
                kfm[i]=0.5*(kf[i]+kf[i + 1])
            #
            #
            for i in range(nrp2-1, nrp4):
                kf[i] = fkc(x[i])
            #
            #
            m = nrp3
            kmr = 0.5*(kf[m] + kf[m + 1])
            kml = 0.5*(kf[m] + kf[m - 1])
            # setup linear system
            for i in range(0, nrp4):
                x[i]  = tfuel[i][k]
                xd[i] = x[i]
            #
            qfd2 = qf*delr2 #*pp(1)
            m = 1
            ad[m] = 4*kf[m]
            au[m] = -4*kf[m]
            b[m] = qf*radial_power_dist[m][k] # qfd2
            #
            i = m
            for m in range(1, nr):
                ri = 1/float(i)
                ad[m] = kfm[i]+kfm[m] + 0.5*(kfm[m] - kfm[i])*ri
                al[m] = -kfm[i]*(1.0 - 0.5*ri)
                au[m] = -kfm[m]*(1.0 + 0.5*ri)
                b[m] = qf*radial_power_dist[m][k] #! qfd2
                i = m
            #
            #
            m = nrp1
            alpha = kgap*(1.0 - kf[m-1[/kf(m))
            ad(m) = 2.0*(kf(m)+kgap*(1.0 + 0.5/nr)) + alpha
            al(m) = -2.0*kf(m)
            au(m) = -2.0*kgap*(1.0 + 0.5/nr) - alpha
            b(m) = qf*radial_power_dist(m,k) # qfd2
            b(m) = qf*delr2
            #
            m = nrp2
            alpha = 2.0*kgap2*(kf(m + 1)/kf(m) - 1.0)
            ad[m] = 8.0*kf(m) + kgap4 - alpha
            al[m] = -kgap4+alpha
            au[m] = -8.0*kf(m)
            b[m]  = 0.0
            m     = nrp3
            ad[m] = 4.0e+0*(kmr + kml) + tworm*(kmr - kml)
            al[m] = -kml*(4 - tworm)
            au[m] = -kmr*(4 + tworm)
            b[m]  = 0
            m = nrp4
            kconv1= htcoef1*tw
            alpha = 2.0e+0*kconv1*(1.0e+0 - kf[m-1]/kf[m])
            kconv = htcoef1*tw*(4.0e+0 + tw/rw)
            ad[m] = 8.0e+0*kf[m]+kconv+alpha
            al[m] = -8.0e+0*kf[m]
            b[m]  = (kconv + alpha)*tcool1
            #
            # solve the tridiagonal system by gauss elimination
            im1 = 1
            for i in range(1, nrp4):
                aldi = al(i)/ad(im1)
                ad(i) = ad(i)-aldi*au(im1)
                b(i) = b(i)-aldi*b(im1)
                im1 = i
            #
            i    = nrp4
            ip1  = nrp4
            x[i] = b[i]/max(ad[i], 1.0e-30)
            for i in range(nrp3,0,-1):
                x[i]=(b[i]-au[i]*x[ip1])/ad[i]
                ip1=i
            #
            errtf = 0.0e+0
            for i in range(0, nrp4):
                tfuel[i][k] = x[i]
                errtf = max(errtf, abs(x[i] - xd[i]))
            #
            if(errtf <= epstf):
                break
            #
            if(itr > 10000):
                is_converged = False
                break
            #
        #
        #
        tfuel[nrp5][k] = ftfavg(x, nr, delr2)
        #
        #
        
        