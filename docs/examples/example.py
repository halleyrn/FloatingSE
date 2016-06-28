#!/usr/bin/env python
# encoding: utf-8

import sys
import os
# just to temporarily change PYTHONPATH without installing
sys.path.append(os.path.expanduser('~') + '/Dropbox/NREL/NREL_WISDEM/src/twister/rotoraero')
from openmdao.main.api import Component, Assembly, convert_units
from openmdao.main.datatypes.api import Float, Array, Enum, Str, Int, Bool
from openmdao.lib.drivers.api import COBYLAdriver,SLSQPdriver
from sparAssembly import sparAssemblyCalculation
import numpy as np
import time
#from utils import filtered_stiffeners_table
from utils import sys_print

def example_OC3():
    """Calculation with properties based mostly on the OC3."""
    example = sparAssemblyCalculation()
    example.tower_base_outer_diameter = 6.0
    example.tower_top_outer_diameter = 3.5
    example.tower_length = 77.6
    example.tower_mass =  249718.
    example.wind_reference_speed = 11.
    example.wind_reference_height = 75.
    example.alpha = 0.110
    example.spar_elevations = [10.,-4.,-12.,-120.]
    example.number_of_sections = 3
    example.example_turbine_size = '3MW'
    example.rotor_diameter = 101.0
    example.RNA_mass = 125000.
    example.RNA_center_of_gravity_x = 5.75
    example.RNA_center_of_gravity_y = 3.5
    example.fairlead_depth = 13. 
    example.scope_ratio = 1.5
    example.pretension_percent = 5.
    example.mooring_diameter = 0.09
    example.number_of_mooring_lines = 3
    example.user_WML = 698.094
    example.water_depth = 320.
    example.mooring_type = 'CHAIN'
    example.anchor_type =  'PILE'
    example.fairlead_offset_from_shell = 0.5
    example.spar_outer_diameter= [6.5,6.5,9.4]
    example.spar.stiffener_curve_fit = True
    example.neutral_axis = 0.21
    #example.stiffener_index = 232
    example.permanent_ballast_height = 3.
    example.fixed_ballast_height = 5.
    example.wall_thickness=[0.05,0.05,0.05]
    example.number_of_rings = [1,4,14]
    example.bulk_head = ['N', 'T', 'B']
    example.load_condition = 'N'
    example.significant_wave_height = 10.820
    example.significant_wave_period = 9.800
    example.run()
    print '-------------OC3---------------'
    sys_print(example)

if __name__ == "__main__":
    example_OC3()
    

