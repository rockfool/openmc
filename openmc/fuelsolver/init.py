"""
openmc.thsolver
==============

A thermal/hydraulic feedback front-end tool.
"""
class Initial(object):
    """
    """
    def __init__(self):
        pass

    def th_setup(self, power, conductance_gap, T_inlet, flow, zmeshes, radii, pitch, n_th_rings):
        #input for TH 1D calculation
        #power             !  MW
        #conductance_gap   !  W/cm2-K
        #T_inlet           !  Celsius degree
        #zmeshes(:)        !  cm
        #radii(:)          !  cm
        #pitch             !  cm
        #flow              !  kg/sec
        #n_th_rings        !  #
        nzth = size(zmeshes)
        nrn = size(radii) - 2
        npint = 1
        ngt = 0
        powfa0 = power
        pfa0 = pitch
        rs  = radii(nrn)*10.                          # radius of fuel in mm
        rw  = radii(nrn + 2)*10.                      # radius of cladding in mm
        tw  = (radii(nrn + 2) - radii(nrn + 1))*10.   # thickness of gap in mm
        if(tw <= 1e-10):
            tw = 1e-10
        rgt = 0.
        tin = T_inlet           #  Celsius degree
        cfrperfa = flow         #  kg/sec
        hgap = conductance_gap  #  W/cm2-K
        nr = n_th_rings         #  #
        nrp1 = nr+1
        nrp2 = nr+2
        nrp3 = nr+3
        nrp4 = nr+4
        nrp5 = nr+5
        initbd()
        allocth0()
        hzth(1:nzth) = zmeshes(1:nzth)   #  cm
        initth()

    def allocth0(self):
        #
        global nrp5
        global nzth
        global hzth 
        global relp 
        global tfuel
        global tdopl
        global tcool
        global dcool
        global hcool
        global r
        global rhou 
        global rhohu
        global u 
        global ud
        global hflux
        global qvol 
        global htcoef
        global qeff
        global radial_power_dist
        global power_profile
        #
        if(len(hzth)>0):
            relp              =[]
            power_profile     =[]
            tdopl             =[]
            tcool             =[]
            dcool             =[]
            radial_power_dist =[]
            hcool             =[]
            rhou              =[]
            rhohu             =[]
            u                 =[]
            ud                =[]
            qeff              =[]
            tfuel             =[]
            hflux             =[]
            qvol              =[]
            htcoef            =[]
            hzth              =[]
            r                 =[]
        #
        hzth  = [0.0]*nzth            #t/h node height
        relp  = [0.0]*nzth            #relative power
        tfuel = [[0.0 for col in range(nzth)] for row in range(nrp5)]
        tdopl = [0.0]*nzth
        tcool = [0.0]*nzth
        dcool = [0.0]*nzth
        hcool = [0.0]*nzth
        r     = [0.0]*nrp5
        rhou  = [0.0]*(nzth+1)
        rhohu = [0.0]*(nzth+1)
        u     = [0.0]*nzth
        ud    = [0.0]*nzth
        hflux = [0.0]*nzth
        qvol  = [0.0]*nzth
        htcoef= [0.0]*nzth
        qeff  = [0.0]*nzth
        radial_power_dist = [[0.0 for col in range(nzth)] for row in range(nrp5)]
        power_profile = [[0.0 for col in range(nzth)] for row in range(nr)]
        #

    def initbd(self):
        global akfuel, akclad
        akfuel(1)   =  1.05
        akfuel(2)   =  0.
        akfuel(3)   =  0.
        akfuel(4)   =  0.
        akfuel(5)   =  2150.
        akfuel(6)   =  73.15
        akclad(1)   =  7.51
        akclad(2)   =  2.09e-2
        akclad(3)   = -1.45e-5
        akclad(4)   =  7.67e-9

    def initth(self):
        #
        global coreheight
        # 
        coreheight=0
        do kth = 1, nzth
            coreheight = coreheight + hzth(kth)
        enddo
        fracdc = 0. 
        #
        # initialize constants
        pfa = 0.01 *pfa0          #cm to m
        rs = rs*0.001             #mm to m
        rw = rw*0.001             #mm to m
        tw = tw*0.001             #mm to m
        rgt = rgt*0.001           #mm to m
        rs2 = rs*rs
        rw2 = rw*rw
        rg = rw - tw
        rg2 = rg*rg
        #
        powfa=powfa0*1e6          !Mw to w
        #
        acf = (pfa**2 - PI*(npint*rw2 + ngt*rgt**2))    !coolant flow area
        afp = npint*PI*rs2                              !total fuel pellet area
        xi = 2. *pi*(npint*rw + ngt*rgt)               !wetted perimeter
        zeta = npint*2*PI*rw                            !heated perimeter
        cfrperfa = cfrperfa                             !coolant flow rate
        zetap = zeta/acf                                !heated perimeter density
        deq = 4. *acf/xi                               !equivalent diameter
        hac = coreheight*0.01                          !active core height in meter
        powlin = powfa/hac                              !linear power density
        fracdf = 1.  - fracdc                          !fuel heat deposit fraction
        #
        # axial nodalization
        #nzthp1=nzth+1                           !number of junctions
        do kth = 1, nzth
            hzth(kth) = hzth(kth)*0.01            !cm to m
        enddo
        #
        # radial nodalization (in fuel pin)
        delr = rs/real(nr, 8)
        delr2 = delr*delr
        tw2 = tw*tw
        delrw = 0.5 *tw
        delrw2 = delrw*delrw
        tworm = tw/(rg+0.5 *tw)
        do i = 1, nrp1
            r(i) = delr*(i - 1)
        end do
        r(nrp2) = rg
        r(nrp3) = rg + delrw
        r(nrp4) = rw
        kgap = hgap*delr
        kgap2 = hgap*tw*rs/rg
        kgap4 = hgap*tw*(4.  - tw/rg)*rs/rg
        #
        # assign inlet condition to all nodes
        tdopin = dsqrt(tin+CKELVIN) !inlet doppler temperature
        din=fdens(tin) !inlet densit
        do k=1,nzth
            do ir=1,nrp1
                tfuel(ir,k)=tin
                tdopl(k)=tdopin
            enddo
            do ir=nrp2,nrp4
                tfuel(ir,k)=tin
            enddo
            tcool(k)=tin
            dcool(k)=din
        enddo
        # initialize vloume and junction variables
        hin = fenthal(tin)
        rhoin = fdens(tin)
        rhouin = cfrperfa/acf
        rhohuin = rhouin*hin
        uin = rhouin/rhoin
        rhou(:) = rhouin
        rhohu(:) = rhohuin
        u(:) = uin
        ud(:) = uin
        hcool(:) = hin
        radial_power_dist = delr2
        #












