from openmdao.main.api import Component, Assembly, convert_units
from openmdao.main.datatypes.api import Float, Array, Enum, Str, Int, Bool
from openmdao.lib.drivers.api import COBYLAdriver,SLSQPdriver
from spar import Spar
from tower_RNA import Tower_RNA
from mapMooring import MapMooring
#from spar_discrete import spar_discrete
import numpy as np
import time
#from utils import filtered_stiffeners_table
from utils import sys_print

class sparAssembly(Assembly):
    """ Top level assembly """
   
    def configure(self):
        """Creates a new Assembly containing a chain of Tower_RNA, Spar and
        Mooring components, as well as a constained optimizer."""

        """Create optimizer instance."""
        self.add('driver',COBYLAdriver())
        self.driver.maxfun = 100000

        """Select component instances."""
        self.add('tower_RNA',Tower_RNA())
        self.add('spar',Spar())
        self.add('mapMooring',MapMooring())

        """Define iteration hierarchy."""
        self.driver.workflow.add(['tower_RNA', 'mapMooring', 'spar'])
        
        """Create a variable in the assembly and connects it to an internal
        component variable. If the variable is used again in a different 
        component instance, then it is manually connected."""
        self.create_passthrough('tower_RNA.base_outer_diameter','tower_base_outer_diameter')
        self.create_passthrough('tower_RNA.top_outer_diameter','tower_top_outer_diameter')
        self.create_passthrough('tower_RNA.length','tower_length')
        self.create_passthrough('tower_RNA.example_turbine_size','example_turbine_size')
        self.create_passthrough('tower_RNA.RNA_center_of_gravity_y', 'RNA_center_of_gravity_y')
        self.create_passthrough('spar.wall_thickness','wall_thickness')
        self.create_passthrough('tower_RNA.rotor_diameter','rotor_diameter')
        self.create_passthrough('tower_RNA.cut_out_speed','cut_out_speed')
        self.create_passthrough('tower_RNA.air_density','air_density')
        self.connect('air_density','spar.air_density')
        self.create_passthrough('spar.wind_reference_speed', 'wind_reference_speed')
        self.connect('wind_reference_speed','tower_RNA.wind_reference_speed')
        self.create_passthrough('spar.wind_reference_height','wind_reference_height')
        self.connect('wind_reference_height','tower_RNA.wind_reference_height')
        self.create_passthrough('spar.gust_factor','gust_factor')
        self.connect('gust_factor','tower_RNA.gust_factor')
        self.create_passthrough('spar.alpha', 'alpha')
        self.connect('alpha','tower_RNA.alpha')
        self.create_passthrough('spar.RNA_center_of_gravity_x','RNA_center_of_gravity_x')
        self.connect('RNA_center_of_gravity_x','tower_RNA.RNA_center_of_gravity_x')
        self.create_passthrough('spar.tower_mass','tower_mass')
        self.connect('tower_mass','tower_RNA.tower_mass')
        self.create_passthrough('spar.RNA_mass','RNA_mass')
        self.connect('RNA_mass','tower_RNA.RNA_mass')
        self.create_passthrough('spar.stiffener_index','stiffener_index')
        self.create_passthrough('spar.number_of_sections','number_of_sections')
        self.create_passthrough('spar.bulk_head','bulk_head')
        self.create_passthrough('spar.number_of_rings','number_of_rings')
        self.create_passthrough('spar.neutral_axis','neutral_axis')
        self.create_passthrough('spar.straight_col_cost','straight_col_cost')
        self.create_passthrough('spar.tapered_col_cost','tapered_col_cost')
        self.create_passthrough('spar.outfitting_cost','outfitting_cost')
        self.create_passthrough('spar.ballast_cost','ballast_cost')
        self.create_passthrough('spar.gravity','gravity')
        self.create_passthrough('spar.load_condition','load_condition')
        self.create_passthrough('spar.significant_wave_height','significant_wave_height')
        self.create_passthrough('spar.significant_wave_period','significant_wave_period')
        self.create_passthrough('spar.material_density','material_density')
        self.create_passthrough('spar.E','E')
        self.create_passthrough('spar.nu','nu')
        self.create_passthrough('spar.yield_stress','yield_stress')
        self.create_passthrough('spar.shell_mass_factor','shell_mass_factor')
        self.create_passthrough('spar.bulkhead_mass_factor','bulkhead_mass_factor')
        self.create_passthrough('spar.ring_mass_factor','ring_mass_factor')
        self.create_passthrough('spar.outfitting_factor','outfitting_factor')
        self.create_passthrough('spar.spar_mass_factor','spar_mass_factor')
        self.create_passthrough('spar.permanent_ballast_height','permanent_ballast_height')
        self.create_passthrough('spar.fixed_ballast_height','fixed_ballast_height')
        self.create_passthrough('spar.permanent_ballast_density','permanent_ballast_density')
        self.create_passthrough('spar.fixed_ballast_density','fixed_ballast_density')
        self.create_passthrough('spar.offset_amplification_factor','offset_amplification_factor')
        self.create_passthrough('spar.water_density','water_density')
        self.create_passthrough('spar.elevations', 'spar_elevations')
        self.create_passthrough('spar.outer_diameter','spar_outer_diameter')
        self.create_passthrough('spar.water_depth','water_depth')
        self.create_passthrough('spar.stiffener_curve_fit', 'stiffener_curve_fit')
        
        #mapMooring connections
        self.create_passthrough('mapMooring.fairlead_depth','fairlead_depth')
        self.connect('spar_elevations',['tower_RNA.spar_elevations','mapMooring.spar_elevations'])
        self.connect('spar_outer_diameter','mapMooring.spar_outer_diameter')
        self.create_passthrough('mapMooring.scope_ratio','scope_ratio')
        self.create_passthrough('mapMooring.mooring_diameter','mooring_diameter')
        self.create_passthrough('mapMooring.number_of_mooring_lines','number_of_mooring_lines')
        self.connect('water_depth','mapMooring.water_depth')
        self.create_passthrough('mapMooring.mooring_type','mooring_type')
        self.create_passthrough('mapMooring.anchor_type','anchor_type')
        self.create_passthrough('mapMooring.fairlead_offset_from_shell','fairlead_offset_from_shell')
        self.create_passthrough('mapMooring.user_MBL','user_MBL')
        self.create_passthrough('mapMooring.user_WML','user_WML')
        self.create_passthrough('mapMooring.user_AE_storm','user_AE_storm')
        self.create_passthrough('mapMooring.user_MCPL','user_MCPL')
        self.create_passthrough('mapMooring.user_anchor_cost','user_anchor_cost')
        self.create_passthrough('mapMooring.misc_cost_factor','misc_cost_factor')
        self.connect('water_density','mapMooring.water_density')
        self.connect('gravity','mapMooring.gravity')

        
        """Connect outputs to inputs."""
        self.connect('tower_RNA.RNA_keel_to_CG','spar.RNA_keel_to_CG')
        self.connect('tower_RNA.tower_center_of_gravity','spar.tower_center_of_gravity')
        self.connect('tower_RNA.tower_wind_force','spar.tower_wind_force')
        self.connect('tower_RNA.RNA_wind_force','spar.RNA_wind_force')

        #mapMooring connections
        self.connect('mapMooring.mooring_total_cost','spar.mooring_total_cost')
        self.connect('mapMooring.mooring_keel_to_CG','spar.mooring_keel_to_CG')
        self.connect('mapMooring.mooring_vertical_load','spar.mooring_vertical_load')
        self.connect('mapMooring.mooring_horizontal_stiffness','spar.mooring_horizontal_stiffness')
        self.connect('mapMooring.mooring_vertical_stiffness','spar.mooring_vertical_stiffness')
        self.connect('mapMooring.sum_forces_x','spar.sum_forces_x')
        self.connect('mapMooring.offset_x','spar.offset_x')
        self.connect('mapMooring.damaged_mooring','spar.damaged_mooring')
        self.connect('mapMooring.intact_mooring','spar.intact_mooring')
        self.connect('mapMooring.mooring_mass','spar.mooring_mass')

       
        """Design variables by adding a range of validity for certain variables."""
        self.driver.add_parameter('neutral_axis',low=10.,high=41.9,scaler=0.01)
        #self.driver.add_parameter('number_of_rings[0]',low=1,high=5)
        self.driver.add_parameter('number_of_rings[1]',low=1,high=10)
        self.driver.add_parameter('number_of_rings[2]',low=1,high=10)
        self.driver.add_parameter('number_of_rings[3]',low=1,high=50)
        self.driver.add_parameter('wall_thickness[0]',low=1.,high=10.,scaler=0.01)
        self.driver.add_parameter('wall_thickness[1]',low=1.,high=10.,scaler=0.01)
        self.driver.add_parameter('wall_thickness[2]',low=1.,high=10.,scaler=0.01)
        self.driver.add_parameter('wall_thickness[3]',low=10.,high=100.,scaler=0.001)
        self.driver.add_parameter('scope_ratio',low=15.,high=45.,scaler=0.1)
        self.driver.add_parameter('pretension_percent',low=2.5,high=10.)
        self.driver.add_parameter('mooring_diameter',low=30.,high=100.,scaler=0.001)
        self.driver.add_parameter('fixed_ballast_height',low=30.,high=100.,scaler=0.1)
        self.driver.add_parameter('permanent_ballast_height',low=30.,high=100.,scaler=0.1)

        """Specify objective function (what you want to minimize)."""
        self.driver.add_objective('spar.spar_mass', name='spar mass')

        """Add constraints to the driver."""
        self.driver.add_constraint('spar.water_ballast_height < 7.5')
        self.driver.add_constraint('spar.water_ballast_height > 5.5')
        self.driver.add_constraint('spar.flange_compactness < 1.')
        self.driver.add_constraint('spar.web_compactness < 1.')
        self.driver.add_constraint('spar.VAL[0] < 0.99')
        self.driver.add_constraint('spar.VAL[1] < 0.99')
        self.driver.add_constraint('spar.VAL[2] < 0.99')
        self.driver.add_constraint('spar.VAL[3] < 0.99')
        self.driver.add_constraint('spar.VAG[0] < 0.99')
        self.driver.add_constraint('spar.VAG[1] < 0.99')
        self.driver.add_constraint('spar.VAG[2] < 0.99')
        self.driver.add_constraint('spar.VAG[3] < 0.99')
        self.driver.add_constraint('spar.VEL[0] < 0.99')
        self.driver.add_constraint('spar.VEL[1] < 0.99')
        self.driver.add_constraint('spar.VEL[2] < 0.99')
        self.driver.add_constraint('spar.VEL[3] < 0.99')
        self.driver.add_constraint('spar.VEG[0] < 0.99')
        self.driver.add_constraint('spar.VEG[1] < 0.99')
        self.driver.add_constraint('spar.VEG[2] < 0.99')
        self.driver.add_constraint('spar.VEG[3] < 0.99')
        self.driver.add_constraint('spar.platform_stability_check < 1.')
        self.driver.add_constraint('spar.heel_angle <= 6.')
        self.driver.add_constraint('spar.min_offset_unity < 1.0')
        self.driver.add_constraint('spar.max_offset_unity < 1.0')

class sparAssemblyCalculation(sparAssembly):
    """This class inherits from the sparAssembly class. This means that it can
    inherit attributes and methods from sparAssembly. So the only difference is
    the configure, which we redefine below.""" 
    def configure(self):
        """Select component instances."""
        self.add('tower_RNA',Tower_RNA())
        self.add('spar',Spar())
        # self.add('mooring',Mooring())
        self.add('mapMooring',MapMooring())

        """Define iteration hierarchy."""
        # self.driver.workflow.add(['tower_RNA', 'mooring', 'spar'])
        self.driver.workflow.add(['tower_RNA', 'mapMooring', 'spar'])

        """Create a variable in the assembly and connects it to an internal
        component variable. If the variable is used again in a different 
        component instance, then it is manually connected."""
        self.create_passthrough('tower_RNA.base_outer_diameter','tower_base_outer_diameter')
        self.create_passthrough('tower_RNA.top_outer_diameter','tower_top_outer_diameter')
        self.create_passthrough('tower_RNA.length','tower_length')
        self.create_passthrough('tower_RNA.example_turbine_size','example_turbine_size')
        self.create_passthrough('tower_RNA.RNA_center_of_gravity_y', 'RNA_center_of_gravity_y')
        self.create_passthrough('spar.wall_thickness','wall_thickness')
        self.create_passthrough('tower_RNA.rotor_diameter','rotor_diameter')
        self.create_passthrough('tower_RNA.cut_out_speed','cut_out_speed')
        self.create_passthrough('tower_RNA.air_density','air_density')
        self.connect('air_density','spar.air_density')
        self.create_passthrough('spar.wind_reference_speed', 'wind_reference_speed')
        self.connect('wind_reference_speed','tower_RNA.wind_reference_speed')
        self.create_passthrough('spar.wind_reference_height','wind_reference_height')
        self.connect('wind_reference_height','tower_RNA.wind_reference_height')
        self.create_passthrough('spar.gust_factor','gust_factor')
        self.connect('gust_factor','tower_RNA.gust_factor')
        self.create_passthrough('spar.alpha', 'alpha')
        self.connect('alpha','tower_RNA.alpha')
        self.create_passthrough('spar.RNA_center_of_gravity_x','RNA_center_of_gravity_x')
        self.connect('RNA_center_of_gravity_x','tower_RNA.RNA_center_of_gravity_x')
        self.create_passthrough('spar.tower_mass','tower_mass')
        self.connect('tower_mass','tower_RNA.tower_mass')
        self.create_passthrough('spar.RNA_mass','RNA_mass')
        self.connect('RNA_mass','tower_RNA.RNA_mass')
        self.create_passthrough('spar.stiffener_index','stiffener_index')
        self.create_passthrough('spar.number_of_sections','number_of_sections')
        self.create_passthrough('spar.bulk_head','bulk_head')
        self.create_passthrough('spar.number_of_rings','number_of_rings')
        self.create_passthrough('spar.neutral_axis','neutral_axis')
        self.create_passthrough('spar.straight_col_cost','straight_col_cost')
        self.create_passthrough('spar.tapered_col_cost','tapered_col_cost')
        self.create_passthrough('spar.outfitting_cost','outfitting_cost')
        self.create_passthrough('spar.ballast_cost','ballast_cost')
        self.create_passthrough('spar.gravity','gravity')
        self.create_passthrough('spar.load_condition','load_condition')
        self.create_passthrough('spar.significant_wave_height','significant_wave_height')
        self.create_passthrough('spar.significant_wave_period','significant_wave_period')
        self.create_passthrough('spar.material_density','material_density')
        self.create_passthrough('spar.E','E')
        self.create_passthrough('spar.nu','nu')
        self.create_passthrough('spar.yield_stress','yield_stress')
        self.create_passthrough('spar.shell_mass_factor','shell_mass_factor')
        self.create_passthrough('spar.bulkhead_mass_factor','bulkhead_mass_factor')
        self.create_passthrough('spar.ring_mass_factor','ring_mass_factor')
        self.create_passthrough('spar.outfitting_factor','outfitting_factor')
        self.create_passthrough('spar.spar_mass_factor','spar_mass_factor')
        self.create_passthrough('spar.permanent_ballast_height','permanent_ballast_height')
        self.create_passthrough('spar.fixed_ballast_height','fixed_ballast_height')
        self.create_passthrough('spar.permanent_ballast_density','permanent_ballast_density')
        self.create_passthrough('spar.fixed_ballast_density','fixed_ballast_density')
        self.create_passthrough('spar.offset_amplification_factor','offset_amplification_factor')
        self.create_passthrough('spar.water_density','water_density')
        self.create_passthrough('spar.elevations', 'spar_elevations')
        self.create_passthrough('spar.outer_diameter','spar_outer_diameter')
        self.create_passthrough('spar.water_depth','water_depth')
        self.create_passthrough('spar.stiffener_curve_fit', 'stiffener_curve_fit')

        #mapMooring connections
        self.create_passthrough('mapMooring.fairlead_depth','fairlead_depth')
        self.connect('spar_elevations',['tower_RNA.spar_elevations','mapMooring.spar_elevations'])
        self.connect('spar_outer_diameter','mapMooring.spar_outer_diameter')
        self.create_passthrough('mapMooring.scope_ratio','scope_ratio')
        self.create_passthrough('mapMooring.mooring_diameter','mooring_diameter')
        self.create_passthrough('mapMooring.number_of_mooring_lines','number_of_mooring_lines')
        self.connect('water_depth','mapMooring.water_depth')
        self.create_passthrough('mapMooring.mooring_type','mooring_type')
        self.create_passthrough('mapMooring.anchor_type','anchor_type')
        self.create_passthrough('mapMooring.fairlead_offset_from_shell','fairlead_offset_from_shell')
        self.create_passthrough('mapMooring.user_MBL','user_MBL')
        self.create_passthrough('mapMooring.user_WML','user_WML')
        self.create_passthrough('mapMooring.user_AE_storm','user_AE_storm')
        self.create_passthrough('mapMooring.user_MCPL','user_MCPL')
        self.create_passthrough('mapMooring.user_anchor_cost','user_anchor_cost')
        self.create_passthrough('mapMooring.misc_cost_factor','misc_cost_factor')
        self.connect('water_density','mapMooring.water_density')
        self.connect('gravity','mapMooring.gravity')
        self.create_passthrough('mapMooring.user_mass_density_air','user_mass_density_air')
        self.create_passthrough('mapMooring.user_EA_stiffness','user_EA_stiffness')
        self.create_passthrough('mapMooring.anchor_radius','anchor_radius')

        
        """Connect outputs to inputs."""
        self.connect('tower_RNA.RNA_keel_to_CG','spar.RNA_keel_to_CG')
        self.connect('tower_RNA.tower_center_of_gravity','spar.tower_center_of_gravity')
        self.connect('tower_RNA.tower_wind_force','spar.tower_wind_force')
        self.connect('tower_RNA.RNA_wind_force','spar.RNA_wind_force')
        
        #mapMooring connections
        self.connect('mapMooring.mooring_total_cost','spar.mooring_total_cost')
        self.connect('mapMooring.mooring_keel_to_CG','spar.mooring_keel_to_CG')
        self.connect('mapMooring.mooring_vertical_load','spar.mooring_vertical_load')
        self.connect('mapMooring.mooring_horizontal_stiffness','spar.mooring_horizontal_stiffness')
        self.connect('mapMooring.mooring_vertical_stiffness','spar.mooring_vertical_stiffness')
        self.connect('mapMooring.sum_forces_x','spar.sum_forces_x')
        self.connect('mapMooring.offset_x','spar.offset_x')
        self.connect('mapMooring.damaged_mooring','spar.damaged_mooring')
        self.connect('mapMooring.intact_mooring','spar.intact_mooring')
        self.connect('mapMooring.mooring_mass','spar.mooring_mass')