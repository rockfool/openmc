"""
openmc.thsolver
==============

A thermal/hydraulic feedback front-end tool.
"""
class Variables(object):
    """
    """

    def __init__(self):
        slef.epstf      = 1.0e-03 
        self.npint      = 0.0
        self.akfuel     = np.array(0.,0.,0.,0.,0.,0.)
        self.alclad     = mp.array(0.,0.,0.,0.)
        self.coreheight = 400.0
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
        self.hgap       = 0. #gap conductance in w/m^2-C
        self.           = 0. #
        self.powlin     = 0. #linear power density
        self.fracdc     = 0. #fraction of heat deposited directly in coolant
        self.fracdf     = 0. #fraction of heat in fuel
        self.din        = 0. #inlet density
        self.tdopin     = 0. #inlet doppler temperature
        self.           = 0. #
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

 
