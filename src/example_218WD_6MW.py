import sys
import os
# just to temporarily change PYTHONPATH without installing
sys.path.append(os.path.expanduser('~') + '/Dropbox/NREL/NREL_WISDEM/src/twister/rotoraero')
from openmdao.main.api import Component, Assembly, convert_units
from openmdao.main.datatypes.api import Float, Array, Enum, Str, Int, Bool
from openmdao.lib.drivers.api import COBYLAdriver,SLSQPdriver
from sparAssemblyWithMAP import sparAssemblyCalculation
import numpy as np
import time
#from utils import filtered_stiffeners_table
from utils import sys_print

def example_218WD_6MW():
    """Executes an optimization of 218WD 6MW."""
    example = sparAssembly()
    example.tower_base_outer_diameter = 6.0
    example.tower_top_outer_diameter = 3.51
    example.tower_length = 80.5
    example.tower_mass =  366952.000
    example.wind_reference_speed = 11.
    example.wind_reference_height = 97.
    example.alpha = 0.110
    #example.spar_lengths = [6.,12.,15.,52.]
    example.spar_elevations = [13.,7.,-5.,-20.,-72.]
    example.example_turbine_size = '6MW'
    example.rotor_diameter = 154.
    example.RNA_mass = 365500.000
    example.RNA_center_of_gravity_x = 5.750
    example.RNA_center_of_gravity_y = 3.5
    example.fairlead_depth = 13.
    example.scope_ratio = 1.5
    example.pretension_percent = 5.0
    example.mooring_diameter = 0.090
    example.number_of_mooring_lines = 3
    example.water_depth = 218.
    example.mooring_type = 'CHAIN'
    example.anchor_type =  'PILE'
    example.fairlead_offset_from_shell = 0.5
    example.spar_outer_diameter= [7.000,8.000,8.000,13.000]
    example.wall_thickness=[0.05,0.05,0.05,0.05]
    example.spar.stiffener_curve_fit = False
    #example.neutral_axis = 0.22
    example.stiffener_index = 271
    example.fixed_ballast_height = 7.0
    example.permanent_ballast_height = 3.0
    #example.wall_thickness=[0.0263,0.0251,0.0262,0.038]
    example.number_of_rings = [1,4,4,19]
    example.scope_ratio = 1.5
    example.pretension_percent = 6.5
    example.mooring_diameter = 0.075
    example.number_of_sections = 4
    example.bulk_head = ['N', 'T', 'N', 'B']
    example.load_condition = 'N'
    example.significant_wave_height = 10.820
    example.significant_wave_period = 9.800
    example.run()
    print '----------218WD_6MW------------'
    sys_print(example)


if __name__ == "__main__":
    example_218WD_6MW()
    