#!/usr/bin/env python
# encoding: utf-8

import sys
import os
from src.sparAssemblyWithMAP import SparAssembly
from src.utils import sys_print
# just to temporarily change PYTHONPATH without installing
sys.path.append(os.path.expanduser('~') + '/Dropbox/NREL/NREL_WISDEM/src/twister/rotoraero')


def example_218wd_10mw():
    """Executes an optimization of 218WD 10MW."""
    example = SparAssembly()
    example.tower_base_outer_diameter = 7.72
    example.tower_top_outer_diameter = 4.050
    example.tower_length = 102.630
    example.tower_mass = 698235.000
    example.wind_reference_speed = 11.
    example.wind_reference_height = 119.
    example.alpha = 0.110
    # example.spar_lengths = [6., 12., 15., 72.]
    example.spar_elevations = [13., 7., -5., -20., -92.]
    example.example_turbine_size = '10MW'
    example.rotor_diameter = 194.
    example.RNA_mass = 677000.000
    example.RNA_center_of_gravity_x = 7.07
    example.RNA_center_of_gravity_y = 3.370
    example.fairlead_depth = 13.
    example.scope_ratio = 1.5
    example.pretension_percent = 5.0
    example.mooring_diameter = 0.090
    example.number_of_mooring_lines = 3
    example.water_depth = 218.
    example.mooring_type = 'CHAIN'
    example.anchor_type = 'PILE'
    example.fairlead_offset_from_shell = 0.5
    example.spar_outer_diameter = [8.0, 9.0, 9.0, 15.0]
    example.wall_thickness = [0.05, 0.05, 0.05, 0.05]
    example.spar.stiffener_curve_fit = False
    example.stiffener_index = 282
    example.fixed_ballast_height = 9.0
    example.permanent_ballast_height = 4.0
    # note: these are slightly off from excel: ie the stiffener AR, which resulted in diff in RGM 
    # example.wall_thickness = [0.0366, 0.035, 0.035, 0.059]
    # example.neutral_axis = 0.35
    example.number_of_rings = [1, 4, 4, 40]
    example.number_of_sections = 4
    example.bulk_head = ['N', 'T', 'N', 'B']
    example.load_condition = 'N'
    example.significant_wave_height = 10.820
    example.significant_wave_period = 9.800
    example.run()
    print '----------218WD_10MW------------'
    sys_print(example)

if __name__ == "__main__":
    example_218wd_10mw()
