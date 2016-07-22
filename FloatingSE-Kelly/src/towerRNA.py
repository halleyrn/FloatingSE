from openmdao.main.api import Component
from openmdao.lib.datatypes.api import Float, Array, Str
from math import pi
from utils import windPowerLaw, windDrag, thrust_table


class TowerRNA(Component):
    """Environmental factor inputs."""
    air_density = Float(1.198, iotype='in', units='kg/m**3', desc='density of air') 
    wind_reference_speed = Float(iotype='in', units='m/s', desc='reference wind speed')
    wind_reference_height = Float(iotype='in', units='m', desc='reference height')
    gust_factor = Float(1.0, iotype='in', desc='gust factor')
    alpha = Float(iotype='in', desc='power law exponent')
    """Additional inputs."""
    base_outer_diameter = Float(iotype='in', units='m', desc='outer diameter of tower base')
    top_outer_diameter = Float(iotype='in', units='m', desc='outer diameter of tower top')
    length = Float(iotype='in', units='m', desc='tower length')
    spar_elevations = Array(iotype='in', units='m', desc='elevations of each section')
    example_turbine_size = Str(iotype='in', desc='for example cases, 3MW, 6MW, or 10 MW')
    rotor_diameter = Float(iotype='in', units='m', desc='rotor diameter')
    RNA_center_of_gravity_x = Float(iotype='in', units='m', desc='rotor center of gravity')
    RNA_center_of_gravity_y = Float(iotype='in', units='m', desc='rotor center of gravity')
    cut_out_speed = Float(25.0, iotype='in', units='m/s', desc='cut-out speed of turbine')
    tower_mass = Float(iotype='in', units='kg', desc='tower mass')
    RNA_mass = Float(iotype='in', units='kg', desc='RNA mass')
    user_tower_cg = Float(0.0, iotype='in', units='m', desc='user defined tower CG ')
    """Outputs."""
    tower_center_of_gravity = Float(iotype='out', units='m', desc='tower center of gravity')
    tower_keel_to_CG = Float(iotype='out', units='m', desc='keel to tower center of gravity')
    tower_wind_force = Float(iotype='out', units='N', desc='wind force on tower')
    RNA_keel_to_CG = Float(iotype='out', units='m', desc='keel to RNA center of gravity')
    RNA_wind_force = Float(iotype='out', units='N', desc='wind force on RNA')
    
    def __init__(self):
        super(TowerRNA, self).__init__()
    
    def execute(self):   

        # tower
        freeboard = self.spar_elevations[0]
        tower_base_outer_diameter = float(self.base_outer_diameter)
        tower_top_outer_diameter = float(self.top_outer_diameter)
        tower_length = float(self.length)
        gust_factor = float(self.gust_factor)
        spar_end_elevations = self.spar_elevations[1:]
        draft = abs(min(spar_end_elevations))
        wind_reference_speed = float(self.wind_reference_speed)
        wind_reference_height = float(self.wind_reference_height)
        alpha = float(self.alpha)
        air_density = float(self.air_density)
        self.tower_center_of_gravity, self.tower_wind_force = \
            windDrag(tower_length, tower_base_outer_diameter, tower_top_outer_diameter, wind_reference_speed,
                     wind_reference_height, alpha, freeboard, air_density, gust_factor, self.user_tower_cg)
        self.tower_keel_to_CG = self.tower_center_of_gravity + freeboard + draft

        # RNA
        rotor_diameter = float(self.rotor_diameter)
        rna_center_of_gravity_y = float(self.RNA_center_of_gravity_y)
        cut_out_speed = float(self.cut_out_speed)
        self.RNA_keel_to_CG = draft + freeboard + tower_length + rna_center_of_gravity_y
        rotor_swept_area = (pi/4.)*rotor_diameter**2
        trust_coefficients, thrust = thrust_table(self.example_turbine_size, air_density, rotor_swept_area)
        max_thrust = max(thrust)
        max_index = thrust.index(max_thrust)
        max_thrust_coefficient = trust_coefficients[max_index]
        hh = rna_center_of_gravity_y + tower_length + freeboard  # what is this variable?
        wind_speed = windPowerLaw(wind_reference_speed, wind_reference_height, alpha, hh)
        if wind_speed < cut_out_speed:
            self.RNA_wind_force = 0.5*air_density*(wind_speed*gust_factor)**2*rotor_swept_area*max_thrust_coefficient
        else:
            self.RNA_wind_force = max_thrust*1000*gust_factor**2*0.75
