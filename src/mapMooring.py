from openmdao.main.api import Component, Assembly,convert_units
from openmdao.lib.datatypes.api import Float, Array, Str, Int, Bool
from openmdao.lib.drivers.api import SLSQPdriver
from numpy import pi, array, arcsin, arcsinh, arc
from scipy.optimize import fmin, minimize, linspace
from sympy.solvers import solve
from sympy import Symbol
from map import InputMAP,

class MapMooring(Component):
    """Creates a mooring component that can be optimized using OpenMDAO.""" 
    water_density = Float(1025, iotype='in',units='kg/m**3',desc='density of water')
    water_depth = Float(iotype='in',units='m',desc='water depth')
    scope_ratio = Float(1.5, iotype='in',units='m',desc = 'scope to fairlead height ratio')
    pretension_percent = Float(5.0, iotype='in',desc='Pre-Tension Percentage of MBL (match PreTension)')
    mooring_diameter = Float(.09, iotype='in',units='m',desc='diameter of mooring chain')
    fairlead_depth = Float(13, iotype='in',units='m',desc = 'fairlead depth')
    number_of_mooring_lines = Int(3, iotype='in',desc='number of mooring lines')
    mooring_type = Str('CHAIN', iotype='in',desc='CHAIN, STRAND, IWRC, or FIBER')
    anchor_type = Str('PILE', iotype='in',desc='PILE or DRAG')
    fairlead_offset_from_shell = Float(.5, iotype='in',units='m',desc='fairlead offset from shell')
    user_MBL = Float(0.0,iotype='in',units='N',desc='user defined minimum breaking load ')
    user_WML = Float(0.0,iotype='in',units='kg/m',desc='user defined wet mass/length')
    user_AE_storm = Float(0.0,iotype='in',units='Pa',desc='user defined E modulus')
    user_MCPL = Float(0.0,iotype='in',units='USD/m',desc='user defined mooring cost per length')
    user_anchor_cost = Float(0.0,iotype='in',units='USD',desc='user defined cost per anchor')
    misc_cost_factor = Float(10.0,iotype='in',desc='miscellaneous cost factor in percent')
    number_of_discretizations = Int(20,iotype='in',desc='number of segments for mooring discretization')
    spar_elevations = Array(iotype='in', units='m',desc = 'end elevation of each section')
    spar_outer_diameter = Array(iotype='in',units='m',desc='top outer diameter')
    gravity = Float(9.806, iotype='in', units='m/s**2', desc='gravity')

    mooring_total_cost = Float(iotype='out',units='USD',desc='total cost for anchor + legs + miscellaneous costs')
    mooring_keel_to_CG = Float(iotype='out',units='m',desc='KGM used in spar.py')
    mooring_vertical_load = Float(iotype='out',units='N',desc='mooring vertical load in all mooring lines')
    mooring_horizontal_stiffness = Float(iotype='out',units='N/m',desc='horizontal stiffness of one single mooring line')
    mooring_vertical_stiffness = Float(iotype='out',units='N/m',desc='vertical stiffness of all mooring lines')
    sum_forces_x = Array(iotype='out',units='N',desc='sume of forces in x direction')
    offset_x = Array(iotype='out',units='m',desc='X offsets in discretization')
    damaged_mooring = Array(iotype='out',units='m',desc='range of damaged mooring')
    intact_mooring = Array(iotype='out',units='m',desc='range of intact mooring')
    mooring_mass = Float(iotype='out',units='kg',desc='total mass of mooring')

    def __init__(self):
        super(MapMooring,self).__init__()
    
    def execute(self):
        """Shows the relationship between each of the variables above."""
        g = self. gravity
        waterDepth = self.water_depth
        fairleadDepth = self.fairlead_depth
        mooringDiameter = self.mooring_diameter
        scopeRatio = self.scope_ratio
        pretensionPercent = self.pretension_percent
        mooringType = self.mooring_type
        numberMooringLines = self.number_of_mooring_lines
        fairleadOffset = self.fairlead_offset_from_shell
        sparOuterDiameter = self.spar_outer_diameter[-1]
        numberDicrestizations= self.number_of_discretizations
        waterDensity = self.water_density
        sparElevations = self.spar_elevations[1:]
        DRAFT = abs(min(sparElevations))
        FH = waterDepth-fairleadDepth 
        scope = FH*scopeRatio
        
        mooring_system = InputMAP(waterDepth, g, waterDensity, numberMooringLines)
        mooring_system.mooring_properties(mooringDiameter, mooringType)
        mooring_system.write_line_dictionary_header()
        #do I want to make a new variable that takes in air mass density and element axial stiffness?
        mooring_system.write_line_dictionary(77.7066, 384243000)

        mooring_system.write_node_properties_header()
        #firgure out how to find out the anchor radius and the fairlead radius
        mooring_system.write_node_properties(1, "FIX", 853.87, 0, -waterDepth, 0, 0)
        mooring_system.write_node_properties(2, "VESSEL", 5.2, 0, -fairleadDepth, 0, 0)
        
        mooring_system.write_line_properties_header()
        mooring_system.write_line_properties(1, mooringType, scope, 1, 2, " ")
        mooring_system.write_solver_options()
        mooring_system.main(2, 2)

        self.sum_forces_x = mooring_system.sum_of_fx()
        self.offset_x = mooring_system.offset_x()
        # COST
        each_leg = MCPL*S
        legs_total = each_leg*numberMooringLines
        if self.anchor_type =='DRAG':
            each_anchor = MBL/1000./9.806/20*2000
        elif self.anchor_type == 'PILE':
            each_anchor = 150000.*(MBL/1000./9.806/1250.)**0.5
        if self.user_anchor_cost != 0.0: 
            each_anchor = self.user_anchor_cost
        anchor_total = each_anchor*numberMooringLines
        misc_cost = (anchor_total+legs_total)*self.misc_cost_factor/100.
        self.mooring_total_cost = legs_total+anchor_total+misc_cost 
        # INITIAL CONDITIONS
        KGM = DRAFT - fairleadDepth 
        VTOP =  np.interp(PTEN,Ttop,Vtop)*numberMooringLines
        MHK = np.interp(PTEN,Ttop,mkh)
        MVK = np.interp(PTEN,Ttop,mkv)*numberMooringLines
        TMM = (WML+pi*mooringDiameter**2/4*waterDensity)*S*numberMooringLines
        self.mooring_keel_to_CG = KGM
        self.mooring_vertical_load = VTOP 
        self.mooring_horizontal_stiffness = MHK
        self.mooring_vertical_stiffness = MVK
        self.mooring_mass = (WML+pi*mooringDiameter**2/4*waterDensity)*S*numberMooringLines
