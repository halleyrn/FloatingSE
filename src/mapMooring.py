from openmdao.main.api import Component
from openmdao.lib.datatypes.api import Float, Array, Str, Int
from math import pi
from map import InputMAP


class MapMooring(Component):
    """Creates a mooring component that can be optimized using OpenMDAO.""" 
    water_density = Float(1025, iotype='in', units='kg/m**3', desc='density of water')
    water_depth = Float(iotype='in', units='m', desc='water depth')
    scope_ratio = Float(1.5, iotype='in', units='m', desc='scope to fairlead height ratio')
    mooring_diameter = Float(.09, iotype='in', units='m', desc='diameter of mooring chain')
    fairlead_depth = Float(13, iotype='in', units='m', desc='fairlead depth')
    number_of_mooring_lines = Int(3, iotype='in', desc='number of mooring lines')
    mooring_type = Str('CHAIN', iotype='in', desc='CHAIN, STRAND, IWRC, or FIBER')
    anchor_type = Str('PILE', iotype='in', desc='PILE or DRAG')
    fairlead_offset_from_shell = Float(.5, iotype='in', units='m', desc='fairlead offset from shell')
    user_MBL = Float(0.0, iotype='in', units='N', desc='user defined minimum breaking load ')
    user_WML = Float(0.0, iotype='in', units='kg/m', desc='user defined wet mass/length')
    user_AE_storm = Float(0.0, iotype='in', units='Pa', desc='user defined E modulus')
    user_MCPL = Float(0.0, iotype='in', units='USD/m', desc='user defined mooring cost per length')
    user_anchor_cost = Float(0.0, iotype='in', units='USD', desc='user defined cost per anchor')
    misc_cost_factor = Float(10.0, iotype='in', desc='miscellaneous cost factor in percent')
    spar_elevations = Array(iotype='in', units='m', desc='end elevation of each section')
    spar_outer_diameter = Array(iotype='in', units='m', desc='top outer diameter')
    gravity = Float(9.806, iotype='in', units='m/s**2', desc='gravity')
    user_mass_density_air = Float(0.0, iotype='in', units='kg/m', desc='user defined mass density in air')
    user_EA_stiffness = Float(0.0, iotype='in', units='N', desc='user defined elemental axial stiffness')
    anchor_radius = Float(iotype='in', units='m', desc='radius to anchors from platform center line')

    mooring_total_cost = Float(iotype='out', units='USD', desc='total cost for anchor + legs + miscellaneous costs')
    mooring_keel_to_CG = Float(iotype='out', units='m', desc='KGM used in spar.py')
    mooring_vertical_load = Float(iotype='out', units='N', desc='mooring vertical load in all mooring lines')
    mooring_horizontal_stiffness = Float(iotype='out', units='N/m', desc='horizontal stiffness of one single '
                                                                         'mooring line')
    mooring_vertical_stiffness = Float(iotype='out', units='N/m', desc='vertical stiffness of all mooring lines')
    sum_forces_x = Array(iotype='out', units='N', desc='sum of forces in x direction')
    offset_x = Array(iotype='out', units='m', desc='X offsets in discretization')
    damaged_mooring = Array(iotype='out', units='m', desc='range of damaged mooring')
    intact_mooring = Array(iotype='out', units='m', desc='range of intact mooring')
    mooring_mass = Float(iotype='out', units='kg', desc='total mass of mooring')

    def __init__(self):
        super(MapMooring, self).__init__()
    
    def execute(self):
        """Shows the relationship between each of the variables above."""
        print 'enter mooring'

        g = float(self. gravity)
        water_depth = float(self.water_depth)
        fairlead_depth = float(self.fairlead_depth)
        mooring_diameter = float(self.mooring_diameter)
        mooring_type = str(self.mooring_type)
        number_of_lines = int(self.number_of_mooring_lines)
        spar_outer_diameter = self.spar_outer_diameter[-1]
        water_density = float(self.water_density)
        spar_elevations = self.spar_elevations[1:]
        draft = abs(min(spar_elevations))
        fairlead_height = water_depth-fairlead_depth
        scope = fairlead_height*float(self.scope_ratio)
        fairlead_offset_from_shell = float(self.fairlead_offset_from_shell)
        misc_cost_factor = float(self.misc_cost_factor)
        fairlead_radius = (spar_outer_diameter/2) + fairlead_offset_from_shell
        
        mooring_system = InputMAP(water_depth, g, water_density, number_of_lines)
        mooring_system.mooring_properties(mooring_diameter, mooring_type, self.user_MBL, self.user_WML,
                                          self.user_AE_storm, self.user_MCPL)
        mooring_system.write_line_dictionary_header()
        mooring_system.write_line_dictionary(self.user_mass_density_air, self.user_EA_stiffness)
        mooring_system.write_node_properties_header()
        mooring_system.write_node_properties(1, "FIX", self.anchor_radius, 0, water_depth, 0, 0)
        mooring_system.write_node_properties(2, "VESSEL", fairlead_radius, 0, -fairlead_depth, 0, 0)
        mooring_system.write_line_properties_header()
        mooring_system.write_line_properties(1, mooring_type, scope, 1, 2, " ")
        mooring_system.write_solver_options()
        mooring_system.main(2, 2, "optimization")

        self.intact_mooring, self.damaged_mooring = mooring_system.intact_and_damaged_mooring()
        self.sum_forces_x, self.offset_x = mooring_system.sum_of_fx_and_offset()
        wml = mooring_system.wet_mass_per_length()
        mcpl = mooring_system.cost_per_length()
        mbl = mooring_system.minimum_breaking_load()
        # COST
        each_leg = mcpl*scope
        legs_total = each_leg*number_of_lines
        each_anchor = self.user_anchor_cost
        if self.anchor_type == 'DRAG':
            each_anchor = mbl/1000./9.806/20*2000
        if self.anchor_type == 'PILE':
            each_anchor = 150000.*(mbl/1000./9.806/1250.)**0.5
        anchor_total = each_anchor*number_of_lines
        misc_cost = (anchor_total+legs_total)*misc_cost_factor/100.
        self.mooring_total_cost = legs_total+anchor_total+misc_cost 
        # INITIAL CONDITIONS
        self.mooring_keel_to_CG = draft - fairlead_depth
        self.mooring_vertical_load, self.mooring_vertical_stiffness, self.mooring_horizontal_stiffness = \
            mooring_system.loads_and_stiffnesses()
        self.mooring_mass = (wml+pi*mooring_diameter**2/4*water_density)*scope*number_of_lines
        print 'end mooring'
