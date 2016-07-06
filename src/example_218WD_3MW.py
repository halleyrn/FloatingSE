#!/usr/bin/env python
# encoding: utf-8

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

def example_218WD_3MW():
    """Executes an optimization of 218WD 3MW."""
    example = sparAssembly()
    example.tower_base_outer_diameter = 4.890
    example.tower_top_outer_diameter = 2.5
    example.tower_length = 60.5
    example.tower_mass =  127877.
    example.wind_reference_speed = 11.
    example.wind_reference_height = 75.
    example.alpha = 0.110
    example.spar_elevations = [13.,7.,-5.,-20.,-67.]
    example.example_turbine_size = '3MW'
    example.rotor_diameter = 101.0
    example.RNA_mass = 125000.
    example.RNA_center_of_gravity_x = 4.1
    example.RNA_center_of_gravity_y = 1.5
    example.fairlead_depth = 13. 
    example.scope_ratio = 1.5
    example.pretension_percent = 5.
    example.mooring_diameter = 0.09
    example.number_of_mooring_lines = 3
    example.water_depth = 218.
    example.mooring_type = 'CHAIN'
    example.anchor_type =  'PILE'
    example.fairlead_offset_from_shell = 0.5
    example.spar_outer_diameter= [5.000,6.000,6.000,9.000]
    example.spar.stiffener_curve_fit = True
    example.neutral_axis = 0.21
    #example.stiffener_index = 232
    example.permanent_ballast_height = 3.
    example.fixed_ballast_height = 5.
    example.wall_thickness=[0.05,0.05,0.05,0.05]
    example.number_of_rings = [1,4,4,14]
    example.number_of_sections = 4
    example.bulk_head = ['N', 'T', 'N', 'B']
    example.load_condition = 'N'
    example.significant_wave_height = 10.820*1.5
    example.significant_wave_period = 9.800
    example.run()
    print '----------218WD_3MW------------'
    sys_print(example)



if __name__ == "__main__":
    example_218WD_3MW()
