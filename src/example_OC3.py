#!/usr/bin/env python
# encoding: utf-8

import sys
import os
from sparAssemblyWithMAP import SparAssemblyCalculation
from utils import sys_print
# just to temporarily change PYTHONPATH without installing
sys.path.append(os.path.expanduser('~') + '/Dropbox/NREL/NREL_WISDEM/src/twister/rotoraero')


def example_oc3():
    """Calculation with properties based mostly on the OC3."""
    example = SparAssemblyCalculation()
    example.example_turbine_size = '5MW'  # not sure if this is correct
    example.neutral_axis = .21  # not sure if this number is correct
    """Cost Variables"""
    example.straight_col_cost = 3492.
    example.tapered_col_cost = 4721.
    example.outfitting_cost = 18651.
    example.ballast_cost = 100.  # not sure if this number is correct
    """Environment Variables"""
    example.air_density = 1.198
    example.gravity = 9.806
    example.water_density = 1025.
    example.water_depth = 320.
    example.load_condition = 'N'
    example.significant_wave_height = 8.
    example.significant_wave_period = 10.
    example.wind_reference_height = 89.350
    example.wind_reference_speed = 11.
    example.alpha = .11
    example.gust_factor = 1.0
    """Material Variables"""
    example.material_density = 7850
    example.E = 200.e9
    example.nu = .3
    example.yield_stress = 345.e6
    example.permanent_ballast_density = 4000.
    example.fixed_ballast_density = 4492.48
    """RNA Variables"""
    example.RNA_mass = 347460.
    example.rotor_diameter = 126.
    example.RNA_center_of_gravity_y = 1.75
    example.RNA_center_of_gravity_x = 1.9
    # example.cut_out_speed 
    """Tower Variables"""
    example.tower_base_outer_diameter = 6.5
    example.tower_top_outer_diameter = 3.87
    example.tower_length = 77.6
    example.tower_mass = 249718.0
    """Stiffener Variables"""
    example.stiffener_index = 259
    example.stiffener_curve_fit = False  # not sure if this is correct
    """Section Variables"""
    example.spar_outer_diameter = [6.5, 6.5, 9.4, 9.4, 9.4]
    example.wall_thickness = [.055, .060, .040, .040, .040]
    example.spar_elevations = [10.0, -4.0, -12.0, -42., -71., -120.]
    example.bulk_head = ['N', 'T', 'N', 'B', 'B']
    example.number_of_rings = [3, 2, 10, 10, 32]
    example.number_of_sections = 5
    """Mooring Variables"""
    example.fairlead_depth = 70.
    example.scope_ratio = 3.609
    example.pretension_percent = 11.173  # map doesnt use
    example.number_of_mooring_lines = 3
    example.mooring_diameter = .09
    example.mooring_type = 'CHAIN'
    example.user_MBL = 8158000.
    example.user_WML = 71.186
    example.user_AE_storm = 384243000./.006
    example.user_MCPL = 0.
    example.anchor_type = 'PILE'
    example.user_anchor_cost = 0.
    example.misc_cost_factor = 10
    example.fairlead_offset_from_shell = .5
    example.number_of_discretizations = 20  # map doesnt use
    example.user_mass_density_air = 77.7066
    example.user_EA_stiffness = 384243000
    example.anchor_radius = 853.87
    """Platform Variables"""
    example.shell_mass_factor = 1
    example.bulkhead_mass_factor = 1.25
    # example.ring_mass_factor
    example.outfitting_factor = .06
    example.spar_mass_factor = 1.04
    example.permanent_ballast_height = 3.
    example.fixed_ballast_height = 7.
    example.offset_amplification_factor = 1
    example.run()
    print '-------------OC3---------------'
    sys_print(example)

if __name__ == "__main__":
    example_oc3()
