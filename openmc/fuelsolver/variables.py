"""
openmc.thsolver
==============

A thermal/hydraulic feedback front-end tool.
"""
from .materials import * 
import numpy as np 
import math 

class Variables(object):
    """Define the common data block for thermal/hydraulic solver
    """
    def __init__(self):
        #
        self.PI = 3.1415926535898e+00
        self.CKELVIN = 273.15e+00
        #
        self.epstf      = 1.0e-03 
        self.npint      = 0.0
        self.akfuel     = np.array([0.,0.,0.,0.,0.,0.])
        self.alclad     = np.array([0.,0.,0.,0.])
        self.coreheight = 400.0        # cm 
        self.relp       = []
        self.power_profile = []      # two dimension 
        self.tdopl      = [] 
        self.tcool      = []
        self.dcool      = []
        self.radial_power_dist = []  # two dimension 
        self.plevel     = 0.
        self.powfa      = 0. #assembly power in watt
        self.powfa0     = 0. #assembly power in Mw
        self.tin        = 0. #inlet temperature
        self.cfrperfa   = 0. #coolant mass flow rate per assy in kg/sec
        self.hgap       = 0. #gap conductance in w/cm^2-C
        
        self.powlin     = 0. #linear power density
        self.fracdc     = 0. #fraction of heat deposited directly in coolant
        self.fracdf     = 0. #fraction of heat in fuel
        self.din        = 0. #inlet density
        self.tdopin     = 0. #inlet doppler temperature
        
        self.tcnom      = 0. #nominal coolant temperature
        self.tfnom      = 0. #nominal fuel temperature
        self.toutavg    = 0. #average outlet temperature
        self.tfmax      = 0. #maximum fuel centerline temperature
        self.fracdvb    = 0. #fraction of heat deposited directly into vessel bypass
        self.fracdwr    = 0. #fraction of heat deposited directly into water rods
        self.byp_dsat   = 0. #saturation density of coolant
        self.byp_a_frac = 0. #bypass area fraction
        self.wr_a_frac  = 0. #water rod area fraction 
        # ============================
        # from thgeom
        # ============================
        self.nperfa     = 0  #number of neutronic nodes per assy
        self.nchan      = 0  #number of total channels
        self.nzth       = 0  #axial t/h node number
        self.ngt        = 0  #number of guide tubes per assy
        self.nzthp1     = 0  #number of junctions
        self.nrn        = 0  #number of rings for neutronics calculation
        self.nr         = 0  #number of nodes in the pellet region
        self.nrp1       = 0  #nr+1, node number at the pellet surface
        self.nrp2       = 0  #nr+2, node number at the cladding inner surface
        self.nrp3       = 0  #nr+3, node number at the middle of the cladding
        self.nrp4       = 0  #nr+4, node number at the cladding outer surface
        self.nrp5       = 0  #nr+5, storage for vol. avg. fuel temperature
        # pellet radius in m (s for surface of pellet)
        self.rs         = 0. #pellet radius in m (s for surface of pellet)
        self.rg         = 0. #clad inner radius in m (s for gap)
        self.rw         = 0. #clad wall radius in m
        self.tw         = 0. #clad thickness in m
        self.rgt        = 0. #guide tube outer radius in m
        self.pfa0       = 0. #assembly pitch in cm
        self.pfa        = 0. #assembly pitch in m
        self.hac        = 0. #active core height
        self.acf        = 0. #coolant flow area
        self.afp        = 0. #fuel pellet area
        self.xi         = 0. #wetted preimeter
        self.zeta       = 0. #heated perimeter
        self.zetap      = 0. #heated perimeter density
        self.deq        = 0. #equivalent diameter
        self.delr       = 0. #radial mesh spacing in the pellet region
        self.delrw      = 0. #radial mesh spacing in the cladding region
        self.tworm      = 0. #tw over rm(=rg+0.5*tw)
        self.rs2        = 0. #rs^2
        self.rw2        = 0. #rw^2
        self.rg2        = 0. #rg^2
        self.tw2        = 0. #tw^2
        self.delr2      = 0. #delr^2
        self.delrw2     = 0. #delrw^2
        # 
        self.hzth       = []
        self.r          = []
        #              
        self.kgap       = 0.
        self.kgapb      = 0.
        self.kgap2      = 0.
        self.kgap4      = 0.
        self.kgap4b     = 0.
        # 
        self.hcool      = []
        self.rhou       = []
        self.rhohu      = []
        self.u          = []
        self.ud         = []
        self.qeff       = []
        #             
        self.tfuel      = [] #(nrp5,nzth)
        self.hflux      = [] #(nzth)
        self.qvol       = [] #(nzth)
        self.htcoef     = [] #(nzth)  


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
        self.nzth = len(zmeshes)
        self.nrn = len(radii) - 2
        self.npint = 1
        self.ngt = 0
        self.powfa0 = power
        self.pfa0 = pitch
        self.rs  = radii[self.nrn-1]*10.                                # radius of fuel in mm
        self.rw  = radii[self.nrn + 2 - 1]*10.                          # radius of cladding in mm
        self.tw  = (radii[self.nrn + 2 - 1] - radii[self.nrn + 1 - 1])*10.   # thickness of gap in mm
        if(self.tw <= 1e-10):
            self.tw = 1e-10
        self.rgt = 0.
        self.tin = T_inlet           #  Celsius degree
        self.cfrperfa = flow         #  kg/sec
        self.hgap = conductance_gap  #  W/cm2-K
        self.nr = int(n_th_rings)    #  #
        self.nrp1 = self.nr+1
        self.nrp2 = self.nr+2
        self.nrp3 = self.nr+3
        self.nrp4 = self.nr+4
        self.nrp5 = self.nr+5
        self.initbd()
        self.allocth0()
        self.hzth[:] = zmeshes[:]   #  cm
        self.initth()

    def allocth0(self):
        #
        if (self.nzth > 0):
            #
            self.hzth  = np.zeros(self.nzth)            #t/h node height
            self.relp  = np.zeros(self.nzth)            #relative power
            self.tfuel = np.zeros((self.nrp5, self.nzth))
            self.tdopl = np.zeros(self.nzth)
            self.tcool = np.zeros(self.nzth)
            self.dcool = np.zeros(self.nzth)
            self.hcool = np.zeros(self.nzth)
            self.r     = np.zeros(self.nrp5)
            self.rhou  = np.zeros(self.nzth+1)
            self.rhohu = np.zeros(self.nzth+1)
            self.u     = np.zeros(self.nzth)
            self.ud    = np.zeros(self.nzth)
            self.hflux = np.zeros(self.nzth)
            self.qvol  = np.zeros(self.nzth)
            self.htcoef= np.zeros(self.nzth)
            self.qeff  = np.zeros(self.nzth)
            self.radial_power_dist = np.zeros((self.nrp5, self.nzth))
            self.power_profile = np.zeros((self.nr, self.nzth))
        #
    
    def initbd(self):
        # define thermal conductivity for fuel and clad
        pass
        # 

    def initth(self):
        #
        #integer     :: kth, i, ir, k
        #real(8)     :: hin
        #real(8)     :: rhoin
        #real(8)     :: rhouin
        #real(8)     :: rhohuin
        #real(8)     :: uin
        #
        self.coreheight = np.sum(self.hzth[:])
        self.fracdc = 0. 
        #
        # initialize constants
        self.pfa = 0.01 *self.pfa0          #cm to m
        self.rs = self.rs*0.001             #mm to m
        self.rw = self.rw*0.001             #mm to m
        self.tw = self.tw*0.001             #mm to m
        self.rgt = self.rgt*0.001           #mm to m
        self.rs2 = self.rs*self.rs
        self.rw2 = self.rw*self.rw
        self.rg = self.rw - self.tw
        self.rg2 = self.rg*self.rg
        #
        self.powfa=self.powfa0*1e6               #MW to W
        #
        self.acf = (self.pfa**2 \
                   - self.PI*(self.npint*self.rw2 + \
                   self.ngt*self.rgt**2))        #coolant flow area
        self.afp = self.npint*self.PI*self.rs2   #total fuel pellet area
        self.xi = 2.*self.PI* \
                  (self.npint*self.rw + self.ngt*self.rgt) #wetted perimeter
        self.zeta = self.npint*2*self.PI*self.rw #heated perimeter
        self.cfrperfa = self.cfrperfa            #coolant flow rate
        self.zetap = self.zeta/self.acf          #heated perimeter density
        self.deq = 4. *self.acf/self.xi          #equivalent diameter
        self.hac = self.coreheight*0.01          #active core height in meter
        self.powlin = self.powfa/self.hac        #linear power density
        self.fracdf = 1. - self.fracdc           #fuel heat deposit fraction
        #
        # axial nodalization
        #nzthp1=nzth+1                           #number of junctions
        self.hzth[:] = self.hzth[:]*0.01         #cm to m
        #
        # radial nodalization (in fuel pin)
        self.delr = self.rs/float(self.nr)
        self.delr2 = self.delr*self.delr
        self.tw2 = self.tw*self.tw
        self.delrw = 0.5 * self.tw
        self.delrw2 = self.delrw*self.delrw
        self.tworm = self.tw/(self.rg+0.5 * self.tw)
        for i in range(1, self.nrp1+1):
            self.r[i-1] = self.delr*(i - 1)
        #
        self.r[self.nrp2] = self.rg
        self.r[self.nrp3] = self.rg + self.delrw
        self.r[self.nrp4] = self.rw
        self.kgap = self.hgap*self.delr
        self.kgap2 = self.hgap*self.tw*self.rs/self.rg
        self.kgap4 = self.hgap*self.tw*(4.- self.tw/self.rg)*self.rs/self.rg
        #
        # assign inlet condition to all nodes
        self.tdopin = math.sqrt(self.tin + self.CKELVIN) #inlet doppler temperature
        self.din=fdens(self.tin) #inlet densit
        for k in range(1, self.nzth+1):
            for ir in range(1, self.nrp1+1):
                self.tfuel[ir-1][k-1]=self.tin
                self.tdopl[k-1]=self.tdopin
            #
            for ir in range(self.nrp2, self.nrp4+1):
                self.tfuel[ir-1][k-1]=self.tin
            #
            self.tcool[k-1]=self.tin
            self.dcool[k-1]=self.din
        #
        # initialize vloume and junction variables
        hin = fenthal(self.tin)
        rhoin = fdens(self.tin)
        rhouin = self.cfrperfa/self.acf
        rhohuin = rhouin*hin
        uin = rhouin/rhoin
        self.rhou[:] = rhouin
        self.rhohu[:] = rhohuin
        self.u[:] = uin
        self.ud[:] = uin
        self.hcool[:] = hin
        self.radial_power_dist[:,:] = self.delr2
        #
