from openmdao.main.api import Component, convert_units
from openmdao.lib.datatypes.api import Float, Array, Str, Int, Bool
from math import cos, sinh, sin, cosh, log, exp, pi
from numpy import array, dot, minimum, asarray, round, isnan, cos, sinh, sin, cosh, interp
from scipy.optimize import fmin
from sympy.solvers import solve
from sympy import Symbol
from utils import full_stiffeners_table, plasticityRF, roots, calcPsi, calculateWindCurrentForces


class Spar(Component):
    """Environmental factor inputs."""
    gust_factor = Float(1.0, iotype='in', desc='gust factor')
    gravity = Float(9.806, iotype='in', units='m/s**2', desc='gravity')
    air_density = Float(1.198, iotype='in', units='kg/m**3', desc='density of air')
    water_density = Float(1025, iotype='in', units='kg/m**3', desc='density of water')
    water_depth = Float(iotype='in', units='m', desc='water depth')
    load_condition = Str(iotype='in', desc='Load condition - number_of_rings for normal or E for extreme')
    significant_wave_height = Float(iotype='in', units='m', desc='significant wave height')
    significant_wave_period = Float(iotype='in', units='m', desc='significant wave period')
    wind_reference_speed = Float(iotype='in', units='m/s', desc='reference wind speed')
    wind_reference_height = Float(iotype='in', units='m', desc='reference height')
    alpha = Float(iotype='in', desc='power law exponent')
    wall_thickness = Array(iotype='in', units='m', desc='wall thickness of each section')
    number_of_rings = Array(iotype='in', desc='number of stiffeners in each section')
    neutral_axis = Float(iotype='in', units='m', desc='neutral axis location')
    """Costs inputs."""
    straight_col_cost = Float(3490, iotype='in', units='USD', desc='cost of straight columns in $/ton')
    tapered_col_cost = Float(4720, iotype='in', units='USD', desc='cost of tapered columns in $/ton')
    outfitting_cost = Float(6980, iotype='in', units='USD', desc='cost of outfitting in $/ton')
    ballast_cost = Float(100, iotype='in', units='USD', desc='cost of ballast in $/ton')
    """Additional inputs."""
    stiffener_curve_fit = Bool(False, iotype='in', desc='flag for using optimized stiffener dimensions or '
                                                        'discrete stiffeners')
    stiffener_index = Int(iotype='in', desc='index of stiffener from filtered table')
    number_of_sections = Int(iotype='in', desc='number of sections in the spar')
    outer_diameter = Array(iotype='in', units='m', desc='outer diameter of each section')
    elevations = Array(iotype='in', units='m', desc='elevations of each section')
    bulk_head = Array(iotype='in', desc='number_of_rings for none, wall_thicknesses for top, B for bottom')
    material_density = Float(7850., iotype='in', units='kg/m**3', desc='density of spar material')
    E = Float(200.e9, iotype='in', units='Pa', desc='young"s modulus of spar material')
    nu = Float(0.3, iotype='in', desc='poisson"s ratio of spar material')
    yield_stress = Float(345000000., iotype='in', units='Pa', desc='yield stress of spar material')
    """Ballast stuff inputs."""
    shell_mass_factor = Float(1.0, iotype='in', desc='shell mass factor')
    bulkhead_mass_factor = Float(1.0, iotype='in', desc='bulkhead mass factor')
    ring_mass_factor = Float(1.0, iotype='in', desc='ring mass factor')
    outfitting_factor = Float(0.06, iotype='in', desc='outfitting factor')
    spar_mass_factor = Float(1.05, iotype='in', desc='spar mass factor')
    permanent_ballast_height = Float(3., iotype='in', units='m', desc='height of permanent ballast')
    fixed_ballast_height = Float(5., iotype='in', units='m', desc='height of fixed ballast')
    permanent_ballast_density = Float(4492., iotype='in', units='kg/m**3', desc='density of permanent ballast')
    fixed_ballast_density = Float(4000., iotype='in', units='kg/m**3', desc='density of fixed ballast')
    offset_amplification_factor = Float(1.0, iotype='in', desc='amplification factor for offsets')
    """Inputs from tower_RNA.py."""
    RNA_keel_to_CG = Float(iotype='in', units='m', desc='RNA keel to center of gravity')
    RNA_mass = Float(iotype='in', units='kg', desc='RNA mass')
    tower_mass = Float(iotype='in', units='kg', desc='tower mass')
    tower_center_of_gravity = Float(iotype='in', units='m', desc='tower center of gravity')
    tower_wind_force = Float(iotype='in', units='N', desc='wind force on tower')
    RNA_wind_force = Float(iotype='in', units='N', desc='wind force on RNA')
    RNA_center_of_gravity_x = Float(iotype='in', units='m', desc='RNA center of gravity in x-direction')
    """Inputs from MapMooring.py."""
    mooring_total_cost = Float(iotype='in', units='USD', desc='total cost for anchor + legs + miscellaneous costs')
    mooring_keel_to_CG = Float(iotype='in', units='m', desc='mooring_keel_to_cg used in spar.py')
    mooring_vertical_load = Float(iotype='in', units='N', desc='mooring vertical load in all mooring lines')
    mooring_horizontal_stiffness = Float(iotype='in', units='N/m', desc='horizontal stiffness of one single '
                                                                        'mooring line')
    mooring_vertical_stiffness = Float(iotype='in', units='N/m', desc='vertical stiffness of all mooring lines')
    sum_forces_x = Array(iotype='in', units='N', desc='sum of forces in x direction')
    offset_x = Array(iotype='in', units='m', desc='X offsets in discretization')
    damaged_mooring = Array(iotype='in', units='m', desc='range of damaged mooring')
    intact_mooring = Array(iotype='in', units='m', desc='range of intact mooring')
    mooring_mass = Float(iotype='in', units='kg', desc='total mass of mooring')
    """Outputs."""
    flange_compactness = Float(iotype='out', desc='check for flange compactness')
    web_compactness = Float(iotype='out', desc='check for web compactness')
    VAL = Array(iotype='out', desc='unity check for axial load - local buckling')
    VAG = Array(iotype='out', desc='unity check for axial load - general instability')
    VEL = Array(iotype='out', desc='unity check for external pressure - local buckling')
    VEG = Array(iotype='out', desc='unity check for external pressure - general instability')
    platform_stability_check = Float(iotype='out', desc='check for platform stability')
    heel_angle = Float(iotype='out', desc='heel angle unity check')
    min_offset_unity = Float(iotype='out', desc='minimum offset unity check')
    max_offset_unity = Float(iotype='out', desc='maximum offset unity check')
    total_cost = Float(iotype='out', units='USD', desc='cost of mooring and spar')
    water_ballast_height = Float(iotype='out', units='m', desc='height of water ballast')
    spar_cost = Float(iotype='out', units='USD', desc='cost of mooring and spar')
    outfit_cost = Float(iotype='out', units='USD', desc='cost of mooring and spar')
    ballasts_cost = Float(iotype='out', units='USD', desc='cost of mooring and spar')
    spar_mass = Float(iotype='out', units='kg', desc='mass of spar')
    ballast_mass = Float(iotype='out', units='kg', desc='ballasts mass')
    system_total_mass = Float(iotype='out', units='kg', desc='total mass of spar system')
    shell_mass = Float(iotype='out', units='kg', desc='total mass of spar system')
    bulkhead_mass = Float(iotype='out', units='kg', desc='total mass of spar system')
    stiffener_mass = Float(iotype='out', units='kg', desc='total mass of spar system')
    total_force = Float(iotype='out', units='N')

    def __init__(self):
        super(Spar, self).__init__()
    
    def execute(self):
        '''
        '''
        # assign all variables
        gravity = float(self.gravity)
        air_density = float(self.air_density)
        water_density = float(self.water_density)
        water_depth = float(self.water_depth)
        load_condition = str(self.load_condition)
        significant_wave_height = float(self.significant_wave_height)
        significant_wave_period = float(self.significant_wave_period)
        waveh = 0
        waven = 0
        if significant_wave_height != 0:
            waveh = 1.86*significant_wave_height
            wavep = 0.71*significant_wave_period
            wavel = gravity*wavep**2/(2*pi)
            waven = 2*pi/wavel
        wind_reference_speed = float(self.wind_reference_speed)
        wind_reference_height = float(self.wind_reference_height)
        alpha = float(self.alpha)
        material_density = float(self.material_density)
        e = float(self.E)
        poisson_ratio = float(self.nu)
        yield_stress = float(self.yield_stress)
        permanent_ballast_height = float(self.permanent_ballast_height)
        permanent_ballast_density = float(self.permanent_ballast_density)
        fixed_ballast_height = float(self.fixed_ballast_height)
        fixed_ballast_density = float(self.fixed_ballast_density)
        rna_mass = float(self.RNA_mass)
        rna_keel_to_cg = float(self.RNA_keel_to_CG)
        tower_mass = float(self.tower_mass)
        tower_center_of_gravity = float(self.tower_center_of_gravity)
        tower_wind_force = float(self.tower_wind_force)
        rna_wind_force = float(self.RNA_wind_force)
        rna_center_of_gravity_x = float(self.RNA_center_of_gravity_x)
        mooring_vertical_load = float(self.mooring_vertical_load)
        mooring_horizontal_stiffness = float(self.mooring_horizontal_stiffness)
        mooring_vertical_stiffness = float(self.mooring_vertical_stiffness)
        mooring_keel_to_cg = float(self.mooring_keel_to_CG)
        outer_diameters = array(self.outer_diameter)
        base_outer_diameter = outer_diameters[-1]
        wall_thicknesses = array(self.wall_thickness)
        end_elevations = array(self.elevations[1:])
        start_elevations = array(self.elevations[:-1])
        number_of_sections = int(self.number_of_sections)
        odtw = 0 
        offset_amplification_factor = float(self.offset_amplification_factor)
        for i in range(0, number_of_sections+1):
            if self.elevations[i] > 0:
                odtw = outer_diameters[i+1]  # What is this variable?
        section_lengths = start_elevations-end_elevations
        draft = abs(min(end_elevations))
        freeboard = start_elevations[0] 
        bulk_heads = self.bulk_head
        number_of_rings = array(self.number_of_rings)
        if self.stiffener_curve_fit:  # curve fits
            stiffener_yna = float(self.neutral_axis)
            stiffener_depth = 0.0029+1.3345977*stiffener_yna
            stiffener_moment_of_inertia = 0.051*stiffener_yna**3.7452
            stiffener_web_thickness = \
                exp(0.88132868+1.0261134*log(stiffener_moment_of_inertia) - 3.117086*log(stiffener_yna))
            stiffener_area = exp(4.6980391+0.36049717*stiffener_yna**0.5-2.2503113/(stiffener_web_thickness**0.5))
            stiffener_flange_thickness = 1.2122528*stiffener_yna**0.13430232*stiffener_yna**1.069737
            stiffener_flange_width = (0.96105249*stiffener_web_thickness**-0.59795001*stiffener_area**0.73163096)
            stiffener_moment_of_inertia = 0.47602202*stiffener_web_thickness**0.99500847*stiffener_yna**2.9938134    
        else:  # discrete, actual stiffener 
            all_stiffeners = full_stiffeners_table()
            stiffener = all_stiffeners[self.stiffener_index]
            stiffener_area = convert_units(stiffener[1], 'inch**2', 'm**2')
            stiffener_depth = convert_units(stiffener[2], 'inch', 'm')
            stiffener_web_thickness = convert_units(stiffener[3], 'inch', 'm')
            stiffener_flange_width = convert_units(stiffener[4], 'inch', 'm')
            stiffener_flange_thickness = convert_units(stiffener[5], 'inch', 'm')
            stiffener_yna = convert_units(stiffener[6], 'inch', 'm')
            self.neutral_axis = stiffener_yna
            stiffener_moment_of_inertia = convert_units(stiffener[7], 'inch**4', 'm**4')
        stiffener_web_height = stiffener_depth - stiffener_flange_thickness
        section_shell_masses, section_ring_masses, section_bulkhead_masses, section_buoyancies, section_wind_forces, \
            section_wind_moments, section_currentwave_forces, section_currentwave_moments, section_keel_to_cgs, \
            section_keel_to_center_buoyancies = \
            calculateWindCurrentForces(0., 0., number_of_rings, stiffener_area, bulk_heads, outer_diameters,
                                       number_of_sections, wall_thicknesses, section_lengths, material_density, draft,
                                       end_elevations, start_elevations, water_density, air_density, gravity,
                                       significant_wave_height, significant_wave_period, water_depth,
                                       wind_reference_speed, wind_reference_height, alpha)
        shell_buoyancy = sum(section_buoyancies)
        shell_mass = sum(section_shell_masses)*float(self.shell_mass_factor)
        bulkhead_mass = sum(section_bulkhead_masses)*float(self.bulkhead_mass_factor)
        ring_mass = sum(section_ring_masses)*float(self.ring_mass_factor)
        shell_bulkhead_ring_mass = shell_mass + bulkhead_mass + ring_mass
        shell_bulkhead_ring_masses = \
            array(section_shell_masses) + array(section_bulkhead_masses) + array(section_ring_masses)
        outfitting_mass = shell_bulkhead_ring_mass * float(self.outfitting_factor)
        mass = shell_bulkhead_ring_mass*float(self.spar_mass_factor) + outfitting_mass
        section_keel_to_cgs = dot(shell_bulkhead_ring_masses, array(section_keel_to_cgs)) / shell_bulkhead_ring_mass
        keel_to_center_of_buoyancy = \
            dot(array(section_buoyancies), array(section_keel_to_center_buoyancies)) / shell_buoyancy
        bm = ((pi / 64) * odtw**4) / (shell_buoyancy/water_density)
        shell_wind_force = sum(section_wind_forces)
        shell_current_force = sum(section_currentwave_forces)
        # NOTE: shell_current_force is inaccurate; setting an initial value and reruns later
        ballast_volume_per_length = pi/4.*(base_outer_diameter-2*wall_thicknesses[-1])**2.
        kgpb = (permanent_ballast_height/2.)+wall_thicknesses[-1] 
        permanent_ballast_mass = ballast_volume_per_length*permanent_ballast_height*permanent_ballast_density
        kgfb = (fixed_ballast_height/2.)+permanent_ballast_height+wall_thicknesses[-1] 
        fixed_ballast_mass = ballast_volume_per_length*fixed_ballast_height*fixed_ballast_density
        tower_keel_to_cg = tower_center_of_gravity+freeboard+draft
        water_ballast_mass = \
            shell_buoyancy - mass - rna_mass - tower_mass - mooring_vertical_load / gravity - fixed_ballast_mass -\
            permanent_ballast_mass
        water_ballast_height = water_ballast_mass / (water_density * ballast_volume_per_length)
        self.water_ballast_height = water_ballast_height
        kgwb = water_ballast_height / 2. + permanent_ballast_height + fixed_ballast_height + wall_thicknesses[-1]
        kgb = \
            (mass * section_keel_to_cgs + water_ballast_mass * kgwb + fixed_ballast_mass * kgfb +
             permanent_ballast_mass * kgpb + tower_mass * tower_keel_to_cg + rna_mass * rna_keel_to_cg) / \
            (mass + water_ballast_mass + fixed_ballast_mass + permanent_ballast_mass + tower_mass+rna_mass)
        kg_tension = \
            (mass * section_keel_to_cgs + water_ballast_mass * kgwb + fixed_ballast_mass * kgfb +
             permanent_ballast_mass * kgpb + tower_mass * tower_keel_to_cg + rna_mass * rna_keel_to_cg +
             mooring_vertical_load / gravity * mooring_keel_to_cg) / shell_buoyancy
        gm = keel_to_center_of_buoyancy + bm - kg_tension
        self.platform_stability_check = kg_tension / keel_to_center_of_buoyancy
        total_mass = mass + rna_mass + tower_mass + water_ballast_mass + fixed_ballast_mass + permanent_ballast_mass
        vd = \
            (rna_wind_force + tower_wind_force + shell_wind_force + shell_current_force) / \
            (mass + rna_mass + tower_mass + fixed_ballast_mass + permanent_ballast_mass + water_ballast_mass)
        section_shell_masses, section_ring_masses, section_bulkhead_masses, section_buoyancies, section_wind_forces, \
            section_wind_moments, section_currentwave_forces, section_currentwave_moments, section_keel_to_cgs, \
            section_keel_to_center_buoyancies = \
            calculateWindCurrentForces(kg_tension, vd, number_of_rings, stiffener_area, bulk_heads, outer_diameters,
                                       number_of_sections, wall_thicknesses, section_lengths, material_density, draft,
                                       end_elevations, start_elevations, water_density, air_density, gravity,
                                       significant_wave_height, significant_wave_period, water_depth,
                                       wind_reference_speed, wind_reference_height, alpha)
        # calculate moments
        rna_wind_moment = rna_wind_force * (rna_keel_to_cg - kg_tension)
        tower_wind_moment = tower_wind_force * (tower_keel_to_cg - kg_tension)
        # costs
        columns_mass = \
            sum(section_shell_masses[1::2]) + sum(section_ring_masses[1::2]) + sum(section_bulkhead_masses[1::2])
        tapered_mass = \
            sum(section_shell_masses[0::2]) + sum(section_ring_masses[0::2]) + sum(section_bulkhead_masses[0::2])
        straight_col_cost = float(self.straight_col_cost)
        tapered_col_cost = float(self.tapered_col_cost)
        outfitting_cost = float(self.outfitting_cost)
        ballast_cost = float(self.ballast_cost)
        mooring_cost = float(self.mooring_total_cost)
        self.spar_cost = straight_col_cost * columns_mass / 1000. + tapered_col_cost * tapered_mass / 1000.
        self.outfit_cost = outfitting_cost * outfitting_mass / 1000.
        self.ballasts_cost = ballast_cost * (fixed_ballast_mass + permanent_ballast_mass) / 1000.
        self.total_cost = self.spar_cost + self.outfit_cost + self.ballasts_cost + mooring_cost

        """SIZING TAB"""
        # [TOP MASS(RNA+TOWER)]
        top_mass = rna_mass + tower_mass 
        kg_top = (rna_mass * rna_keel_to_cg + tower_mass * tower_keel_to_cg)/top_mass
        # [INERTIA PROPERTIES - LOCAL]
        inertia_top_local = (1./12.) * top_mass * kg_top**2
        inertia_hull_local = (1./12.) * mass * (draft + freeboard)**2
        inertia_variable_ballast_local = (1./12.) * water_ballast_mass * water_ballast_height**2
        inertia_fixed_ballast_local = (1./12.) * fixed_ballast_mass * fixed_ballast_height**2
        inertia_permanent_ballast_local = (1./12.) * permanent_ballast_mass * permanent_ballast_height**2
        # [INERTIA PROPERTIES - SYSTEM]
        inertia_top_system = inertia_top_local + top_mass * (kg_top - kg_tension)**2
        inertia_hull_system = inertia_hull_local + mass * (section_keel_to_cgs - kg_tension)**2
        inertia_variable_ballast_system = inertia_variable_ballast_local + water_ballast_mass * (kgwb - kgb)**2 
        inertia_fixed_ballast_system = inertia_fixed_ballast_local + fixed_ballast_mass * (kgfb - kgb)**2
        inertia_permanent_ballast_system = inertia_permanent_ballast_local + permanent_ballast_mass*(kgpb-kgb)**2
        inertia_total_system = \
            inertia_top_system + inertia_hull_system + inertia_variable_ballast_system + \
            inertia_fixed_ballast_system + inertia_permanent_ballast_system
        # [ADDED MASS]
        surge = (pi/4.) * (base_outer_diameter**2) * draft * water_density
        heave = (1/6.) * water_density * base_outer_diameter**3
        pitch = \
            (surge * ((kg_tension - draft) - (keel_to_center_of_buoyancy - draft))**2 +
             surge * draft**2/12.) + inertia_total_system
        # [PLATFORM STIFFNESS]
        k33_heave = water_density * gravity * ((pi/4.) * odtw**2) + mooring_vertical_stiffness  # heave
        # [PERIOD]
        surge_t = 2 * pi * ((total_mass + surge) / mooring_horizontal_stiffness)**0.5
        heave_t = 2 * pi * ((total_mass + heave) / k33_heave)**0.5
        pitch_k = gm * shell_buoyancy * gravity
        pitch_t = 2 * pi * (pitch / pitch_k)**0.5
        total_force = rna_wind_force + tower_wind_force + sum(section_wind_forces) + sum(section_currentwave_forces)
        self.total_force = total_force
        sum_forces_x = array(self.sum_forces_x)
        offset_x = self.offset_x
        if isnan(sum_forces_x).any() or sum_forces_x[-1] > (-total_force/1000.):
            self.max_offset_unity = 10.
            self.min_offset_unity = 10.
        else:
            max_offset = interp(array(-total_force/1000), sum_forces_x[::-1], offset_x[::-1])
            min_offset = interp(array(total_force/1000), sum_forces_x[::-1], offset_x[::-1])
        # unity checks! 
            if self.load_condition == 'E':
                self.max_offset_unity = max_offset/self.damaged_mooring[1]
                self.min_offset_unity = min_offset/self.damaged_mooring[0]
            elif self.load_condition == 'N':
                self.max_offset_unity = max_offset/self.intact_mooring[1]
                self.min_offset_unity = min_offset/self.intact_mooring[0]
        m_total = \
            rna_wind_moment + tower_wind_moment + sum(section_wind_moments) + sum(section_currentwave_moments) + \
            (-total_force * (mooring_keel_to_cg - kg_tension)) + (rna_mass * gravity * -rna_center_of_gravity_x)
        self.heel_angle = (m_total / pitch_k) * 180./pi
        """ API BULLETIN """
        # shell data
        outer_radii = array(outer_diameters / 2.)  # outer radius
        R = array(outer_radii-wall_thicknesses / 2.)  # radius to centerline of wall/mid fiber radius
        # ring data
        LR = array(section_lengths / (number_of_rings + 1.))  # number of ring spacing
        # shell and ring data
        radii_to_flange = outer_radii - stiffener_web_height  # radius to flange
        MX = array(LR / (R * wall_thicknesses)**0.5)  # geometry parameter
        # effective width of shell plate in longitudinal direction
        shell_longitudinal_effective_width = array([0.] * number_of_sections)
        for i in range(0, number_of_sections):
            if MX[i] <= 1.56:
                shell_longitudinal_effective_width[i] = LR[i]
            else:
                shell_longitudinal_effective_width = 1.1 * (2 * R * wall_thicknesses)**0.5 + stiffener_web_thickness
        # ring properties with effective shell plate
        shell_effective_area = stiffener_area + shell_longitudinal_effective_width * wall_thicknesses
        YENA = \
            (shell_longitudinal_effective_width * wall_thicknesses * wall_thicknesses / 2 + stiffener_web_height *
             stiffener_web_thickness * (stiffener_web_height / 2 + wall_thicknesses) + stiffener_flange_thickness *
             stiffener_flange_width * (stiffener_flange_thickness / 2 + stiffener_web_height + wall_thicknesses)) /\
            shell_effective_area
        IER = \
            array(stiffener_moment_of_inertia + stiffener_area * (stiffener_yna + wall_thicknesses / 2.)**2 *
                  shell_longitudinal_effective_width * wall_thicknesses / shell_effective_area +
                  shell_longitudinal_effective_width * (wall_thicknesses**3) / 12.)  # moment of inertia
        RC = array(outer_radii - YENA - wall_thicknesses / 2.)  # radius to centroid of ring stiffener
        # set loads (0 mass loads for external pressure)
        ballast_mass = permanent_ballast_mass + fixed_ballast_mass + water_ballast_mass  # sum of all ballast masses
        weight = (rna_mass + tower_mass + ballast_mass + mass) * gravity
        pressure = water_density * gravity * abs(end_elevations)  # hydrostatic pressure at depth of section bottom
        if significant_wave_height != 0:  # dynamic head
            DH = waveh/2*(cosh(waven*(water_depth-abs(end_elevations)))/cosh(waven*water_depth))
        else:
            DH = 0
        pressure += water_density * gravity * DH  # hydrostatic pressure + dynamic head
        """RING SECTION COMPACTNESS (SECTION 7)"""
        self.flange_compactness = (0.5*stiffener_flange_width/stiffener_flange_thickness)/(0.375*(e/yield_stress)**0.5)
        self.web_compactness = (stiffener_web_height/stiffener_web_thickness)/((e/yield_stress)**0.5)
        """PLATE AND RING STRESS (SECTION 11)"""
        # shell hoop stress at ring midway
        Dc = e * (wall_thicknesses**3) / (12 * (1 - poisson_ratio**2))  # parameter stiffener_depth
        BETAc = (e * wall_thicknesses / (4 * (outer_radii**2) * Dc))**0.25  # parameter beta
        TWS = stiffener_area / stiffener_web_height
        dum1 = BETAc * LR
        KT = 8 * BETAc**3 * Dc * (cosh(dum1) - cos(dum1)) / (sinh(dum1) + sin(dum1))
        KD = \
            e * TWS * (outer_radii**2 - radii_to_flange**2) / \
            (outer_radii * ((1 + poisson_ratio) * outer_radii**2 + (1 - poisson_ratio) * radii_to_flange**2))
        dum = dum1 / 2.
        PSIK = array(2 * (sin(dum) * cosh(dum) + cos(dum) * sinh(dum)) / (sinh(dum1) + sin(dum1)))
        PSIK = PSIK.clip(min=0)  # psik >= 0; set all negative values of psik to zero
        SIGMAXA = -weight / (2 * pi * R * wall_thicknesses)
        PSIGMA = pressure + (poisson_ratio * SIGMAXA * wall_thicknesses) / outer_radii
        PSIGMA = minimum(PSIGMA, pressure)  # PSIGMA has to be <= pressure
        dum = KD / (KD + KT)
        KTHETAL = 1 - PSIK * PSIGMA / pressure * dum
        FTHETAS = KTHETAL * pressure * outer_radii / wall_thicknesses
        # shell hoop stress at ring
        KTHETAG = 1 - (PSIGMA / pressure * dum)
        FTHETAR = KTHETAG * pressure * outer_radii / wall_thicknesses
        """LOCAL BUCKLING (SECTION 4)"""
        # axial compression and bending
        ALPHAXL = 9 / (300 + (2 * R) / wall_thicknesses)**0.4
        CXL = (1 + (150 / ((2 * R) / wall_thicknesses)) * (ALPHAXL**2) * (MX**4))**0.5
        FXEL = array(CXL * (pi**2 * e / (12 * (1 - poisson_ratio**2))) * (wall_thicknesses / LR)**2)  # elastic
        FXCL = array(number_of_sections * [0.])
        for i in range(0, len(FXEL)):
            FXCL[i] = plasticityRF(FXEL[i], yield_stress)  # inelastic
        # external pressure
        BETA = array([0.] * number_of_sections)
        ALPHATHETAL = array([0.] * number_of_sections)
        ZM = array(12 * (MX**4 * (1 - poisson_ratio**2)) / pi**4)
        for i in range(0, number_of_sections):
            f = lambda x: x**2 * (1 + x**2)**4 / (2 + 3 * x**2) - ZM[i]
            ans = roots(f, 0., 15.)
            ans_array = array(asarray(ans))
            is_scalar = True
            if ans_array.ndim > 0:
                is_scalar = False
            if is_scalar:
                BETA[i] = ans
            else:
                BETA[i] = float(min(ans_array))
            if MX[i] < 5:
                ALPHATHETAL[i] = 1
            elif MX[i] >= 5:
                ALPHATHETAL[i] = 0.8
        n = round(BETA * pi * R / LR)  # solve for smallest whole number n
        BETA = LR / (pi * R / n)
        left = (1 + BETA**2)**2 / (0.5 + BETA**2)
        right = 0.112 * MX**4 / ((1 + BETA**2)**2 * (0.5 + BETA**2))
        CTHETAL = (left + right) * ALPHATHETAL
        FREL = array(CTHETAL * pi**2 * e * (wall_thicknesses / LR)**2 / (12 * (1 - poisson_ratio**2)))  # elastic
        FRCL = array(number_of_sections * [0.])
        for i in range(0, len(FREL)):
            FRCL[i] = plasticityRF(FREL[i], yield_stress)  # inelastic
        """GENERAL INSTABILITY (SECTION 4)"""
        # axial compression and bending
        AC = array(stiffener_area / (LR * wall_thicknesses))  # Ar bar
        ALPHAX = array(0.85 / (1 + 0.0025 * (outer_diameters / wall_thicknesses)))
        ALPHAXG = array([0.] * number_of_sections)
        for i in range(0, number_of_sections):
            if AC[i] >= 0.2:
                ALPHAXG[i] = 0.72
            elif 0.2 > AC[i] > 0.06:
                ALPHAXG[i] = (3.6-0.5*ALPHAX[i])*AC[i]+ALPHAX[i]
            else:
                ALPHAXG[i] = ALPHAX[i]
        FXEG = array(ALPHAXG * 0.605 * e * wall_thicknesses / R * (1 + AC)**0.5)  # elastic
        FXCG = array(number_of_sections * [0.])
        for i in range(0, len(FXEG)):
            FXCG[i] = plasticityRF(FXEG[i], yield_stress)  # inelastic
        # external pressure
        LAMBDAG = array(pi * R / section_lengths)
        k = 0.5
        PEG = array([0.] * number_of_sections)
        for i in range(0, number_of_sections):
            t = wall_thicknesses[i]
            r = R[i]
            lambdag = LAMBDAG[i]
            ier = IER[i]
            rc = RC[i]
            ro = outer_radii[i]
            lr = LR[i]

            def f(x, e, t, r, lambdag, k, ier, rc, ro, lr):
                return e*(t/r)*lambdag**4/((x**2+k*lambdag**2-1)*(x**2+lambdag**2)**2)+e*ier*(x**2-1)/(lr*rc**2*ro)
            x0 = [2]
            # solve for n that gives minimum P_eG
            m = float(fmin(f, x0, xtol=1e-3, args=(e, t, r, lambdag, k, ier, rc, ro, lr)))
            PEG[i] = f(m, e, t, r, lambdag, k, ier, rc, ro, lr)
        ALPHATHETAG = 0.8  # adequate for ring stiffeners
        FREG = array(ALPHATHETAG * PEG * outer_radii * KTHETAG / wall_thicknesses)  # elastic
        FRCG = array(number_of_sections * [0.])
        for i in range(0, len(FREG)):
            FRCG[i] = plasticityRF(FREG[i], yield_stress)  # inelastic
        # General Load Case
        NPHI = weight / (2 * pi * R)
        NTHETA = pressure * outer_radii 
        K = NPHI/NTHETA
        """Local Buckling (SECTION 6) - Axial Compression and bending"""
        C = array((FXCL + FRCL) / yield_stress - 1)
        KPHIL = 1
        CST = array(K * KPHIL / KTHETAL)
        FTHETACL = array([0.] * number_of_sections)  # external pressure
        for i in range(0, number_of_sections):
            cst = CST[i]
            fxcl = FXCL[i]
            frcl = FRCL[i]
            c = C[i]
            x = Symbol('x')
            ans = solve((cst * x / fxcl)**2 - c * (cst * x / fxcl) * (x / frcl) + (x / frcl)**2 - 1, x)
            FTHETACL[i] = float(min([a for a in ans if a > 0]))
        FPHICL = CST * FTHETACL
        """General Instability (SECTION 6) - Axial Compression and bending"""
        C = array((FXCG + FRCG) / yield_stress - 1)
        KPHIG = 1
        CST = array(K * KPHIG / KTHETAG)
        FTHETACG = array([0.]*number_of_sections)
        for i in range(0, number_of_sections):
            cst = CST[i]
            fxcg = FXCG[i]
            frcg = FRCG[i]
            c = C[i]
            x = Symbol('x', real=True)
            ans = solve((cst * x / fxcg)**2 - c * (cst * x / fxcg) * (x / frcg) + (x / frcg)**2 - 1, x)
            FTHETACG[i] = float(min([a for a in ans if a > 0]))
        FPHICG = CST * FTHETACG
        """Allowable Stresses"""
        # factor of safety
        FOS = 1.25
        if load_condition.upper() == 'N':
            FOS = 1.65
        FAL = array([0.] * number_of_sections)
        FAG = array([0.] * number_of_sections)
        FEL = array([0.] * number_of_sections)
        FEG = array([0.] * number_of_sections)
        for i in range(0, number_of_sections):
            # axial load    
            FAL[i] = FPHICL[i] / (FOS * calcPsi(FPHICL[i], yield_stress))
            FAG[i] = FPHICG[i] / (FOS * calcPsi(FPHICG[i], yield_stress))
            # external pressure
            FEL[i] = FTHETACL[i] / (FOS * calcPsi(FTHETACL[i], yield_stress))
            FEG[i] = FTHETACG[i] / (FOS * calcPsi(FTHETACG[i], yield_stress))
        # unity check 
        self.VAL = abs(SIGMAXA / FAL)
        self.VAG = abs(SIGMAXA / FAG)
        self.VEL = abs(FTHETAS / FEL)
        self.VEG = abs(FTHETAS / FEG)

        print 'surge period: ', surge_t
        print 'heave period: ', heave_t
        print 'pitch stiffness: ', pitch_k
        print 'pitch period: ', pitch_t
        print 'stiffener_yna: ', stiffener_yna
        print 'number of stiffeners: ', self.number_of_rings
        print 'wall thickness: ', self.wall_thickness
        print 'VAL: ', self.VAL
        print 'VAG: ', self.VAG
        print 'VEL: ', self.VEL
        print 'VEG: ', self.VEG
        print 'web compactness: ', self.web_compactness
        print 'flange compactness: ', self.flange_compactness
        print 'heel angle: ', self.heel_angle
        print 'outer diameters: ', self.outer_diameter
        self.spar_mass = shell_bulkhead_ring_mass
        self.ballast_mass = ballast_mass
        self.system_total_mass = \
            shell_bulkhead_ring_mass + permanent_ballast_mass + fixed_ballast_mass \
            + water_ballast_mass + float(self.mooring_mass)
        self.shell_mass = shell_mass 
        self.bulkhead_mass = bulkhead_mass
        self.stiffener_mass = ring_mass
        print 'spar mass: ', self.spar_mass
        print 'shell mass: ', self.shell_mass
        print 'bulkhead mass: ', self.bulkhead_mass
        print 'stiffener mass: ', self.stiffener_mass
        print 'end spar'
