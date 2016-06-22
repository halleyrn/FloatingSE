from openmdao.main.api import Component, Assembly,convert_units
from openmdao.lib.datatypes.api import Float, Array, Str, Int, Bool
from openmdao.lib.drivers.api import SLSQPdriver
from numpy import pi, array, arcsin, arcsinh, arccosh, linspace, interp, arccosh, sin, empty, object, power, cos
from scipy.optimize import fmin, minimize
from sympy.solvers import solve
from sympy import Symbol
from input_map import InputMAP
import os

class Map(Component):
    """Creates a mooring component that can be optimized using MAP and 
    OpenMDAO.""" 
    water_density = Float(1025, iotype='in',units='kg/m**3',desc='density of water')
    water_depth = Float(iotype='in',units='m',desc='water depth')
    scope_ratio = Float(1.5, iotype='in',units='m',desc = 'scope to fairlead height ratio')
    pretension_percent = Float(5.0, iotype='in',desc='Pre-Tension Percentage of minBreakingLoad (match PreTension)')
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
        super(Map,self).__init__()
    
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
        fairlead2seafloor = waterDepth-fairleadDepth 
        scope = fairlead2seafloor*scopeRatio
        
        if mooringType == 'CHAIN':    
            minBreakingLoad = 27600.*(mooringDiameter**2)*(44.-80.*mooringDiameter)*(10**3)
            wetMassPerLength = 18070.*(mooringDiameter**2)
            AE_storm = (1.3788*(mooringDiameter**2)-4.93*(mooringDiameter**3))*(10**11)
            area = 2.64*(mooringDiameter**2)
            costPerLength = 0.58*(minBreakingLoad/1000./g)-87.6
            massInWater = 18070*(mooringDiameter**2)
        elif mooringType == 'STRAND':
            minBreakingLoad = (937600*(mooringDiameter**2)-1408.3*mooringDiameter)*(10**3)
            wetMassPerLength = 4110*(mooringDiameter**2)
            AE_storm = 9.28*(mooringDiameter**2)*(10**10)
            area = 0.58*(mooringDiameter**2)
            costPerLength = 0.42059603*(minBreakingLoad/1000./g)+109.5
            massInWater = 4110*(mooringDiameter**2)
        elif mooringType == 'IWRC':
            minBreakingLoad = 648000*(mooringDiameter**2)*(10**3)
            wetMassPerLength = 3670*(mooringDiameter**2)
            AE_storm = 6.01*(mooringDiameter**2)*(10**10)
            area = 0.54*(mooringDiameter**2)
            costPerLength = 0.33*(minBreakingLoad/1000./g)+139.5
            massInWater = 3670*(mooringDiameter**2)
        elif mooringType == 'FIBER': 
            minBreakingLoad = (274700*(mooringDiameter**2)+7953.9*mooringDiameter-879.24)*(10**3)
            wetMassPerLength = 160.9*(mooringDiameter**2)+5.522*mooringDiameter-0.04798
            AE_storm = (10120*(mooringDiameter**2)+320.7*mooringDiameter-35.47)*(10**6)
            AE_drift = (5156*(mooringDiameter**2)+142.7*mooringDiameter-16)*(10**6)
            area = (pi/4)*(mooringDiameter**2)
            costPerLength = 0.53676471*(minBreakingLoad/1000./g)
            massInWater = 160.9*(mooringDiameter**2)+5.522*(mooringDiameter)-0.04798
        else: 
            print "PLEASE PICK AVAILABLE MOORIN' TYPE M8"
        # if user defined values available
        if self.user_MBL != 0.0:
            minBreakingLoad = self.user_MBL
        if self.user_WML != 0.0: 
            wetMassPerLength = self.user_WML
        if self.user_AE_storm != 0.0: 
            AE_storm = self.user_AE_storm 
        if self.user_MCPL != 0.0: 
            costPerLength = self.user_MCPL
        pretension = minBreakingLoad*pretensionPercent/100.
        wetWeightPerLength = wetMassPerLength*g
        
        #I dont think this is necessary with MAP++
        max_a =(minBreakingLoad-wetWeightPerLength*fairlead2seafloor)/wetWeightPerLength
        empty_array = array(linspace(0.0, max_a, num=755))
        H = wetWeightPerLength*empty_array #wet weight per length array
        Ttop = H + wetWeightPerLength*fairlead2seafloor
        Vtop = (Ttop**2-H**2)**0.5
        ang = 90. - arcsin(H/Ttop)*180/pi
        sp_temp = -(scope/2.)+(fairlead2seafloor/2.)*(1.+(4*empty_array**2/(scope**2.-fairlead2seafloor**2)))**0.5 
        sp_temp = [0 if i < 0 else i for i in sp_temp]
        # INITIALIZE ARRanchor_yS
        Vbot = array([0.]*755)
        Tbot = array([0.]*755)
        Tave = array([0.]*755)
        stretch = array([0.]*755)
        x = array([0.]*755)
        s = array([0.]*755)
        yp = array([0.]*755)
        sp = array([0.]*755)
        for i in range (0,755):
            if sp_temp[i] == 0.: 
                Vbot[i] = 0.
                Tbot[i] = H[i]
                Tave[i] = 0.5*(Ttop[i]+Tbot[i])
                stretch[i] = 1.+Tave[i]/AE_storm
                if i == 0:
                    x[i] = scope*stretch[i] -fairlead2seafloor
                else: 
                    x[i] = scope*stretch[i-1] -fairlead2seafloor*(1.+2.*empty_array[i]/fairlead2seafloor)**0.5+empty_array[i]*arccosh(1.+fairlead2seafloor/empty_array[i])
                s[i] = (fairlead2seafloor**2+2.*fairlead2seafloor*empty_array[i])**0.5
                sp [i] = 0.
            else: 
                s[i] = scope*stretch[i-1]
                Vbot[i] = Vtop[i]-wetWeightPerLength*s[i]
                Tbot[i] = (H[i]**2+Vbot[i]**2)**0.5
                Tave[i] =  0.5*(Ttop[i]+Tbot[i])
                stretch[i] = 1.+Tave[i]/AE_storm
                sp[i] = sp_temp[i]*stretch[i-1]
                x[i] = empty_array[i]*(arcsinh((scope+sp[i])/empty_array[i])-arcsinh(sp[i]/empty_array[i]))*stretch[i-1]
                
            if i == 0: 
                yp[i] = -empty_array[i]+(empty_array[i]**2+sp[i]**2)**0.5
            else: 
                yp[i] = (-empty_array[i]+(empty_array[i]**2+sp[i]**2)**0.5)*stretch[i-1]
        x0 = interp(pretension,Ttop,x)
        offset = x - x0
        mkh = array([0.]*755)
        mkv = array([0.]*755)
        for i in range(1,755):
            mkh[i] = (H[i]-H[i-1])/(offset[i]-offset[i-1])
            mkv[i] = (Vtop[i]-Vtop[i-1])/(offset[i]-offset[i-1])

        if interp(pretension,Ttop,sp)>0.:
            cat_type = 'semi-taut'
        else: 
            cat_type = 'catenary'
        XANG = interp(pretension,Ttop,ang)
        XMAX = max(x)
        TALL = 0.6*minBreakingLoad 
        XALL = interp(TALL,Ttop,x)
        TEXT = 0.8*minBreakingLoad 
        XEXT = interp(TEXT,Ttop,x)
        # OFFSETS 
        direction =  array(linspace(0.0,360-360/numberMooringLines,num=numberMooringLines))
        self.intact_mooring = [-(XALL-x0), XALL/sin(direction[1]/180.*pi)*sin(pi-direction[1]/180*pi-arcsin(x0/XALL*sin(direction[1]/180.*pi)))]
        self.damaged_mooring = [-(XEXT-x0), XEXT/sin(direction[1]/180.*pi)*sin(pi-direction[1]/180*pi-arcsin(x0/XEXT*sin(direction[1]/180.*pi)))]
        survival_mooring = [-(XMAX-x0), XMAX/sin(direction[1]/180.*pi)*sin(pi-direction[1]/180*pi-arcsin(x0/XMAX*sin(direction[1]/180.*pi)))]

        # fairlead
        fairlead_x = array(((sparOuterDiameter/2)+fairleadOffset)*cos(direction*pi/180.))
        fairlead_y = array(((sparOuterDiameter/2)+fairleadOffset)*sin(direction*pi/180.))
        fairlead_z = array([-fairleadDepth]*len(fairlead_x))
        # anchor 
        anchor_x = array(fairlead_x+x0*cos(direction*pi/180.))
        anchor_y = array(fairlead_y+x0*sin(direction*pi/180.))
        anchor_z = array([-waterDepth]*len(anchor_x))
        # delta 
        delta = (survival_mooring[1]-survival_mooring[0])/numberDicrestizations
        X_Offset = array(linspace(survival_mooring[0],survival_mooring[1],num=numberDicrestizations+1))
        # initialize arrays
        X_Fairlead = empty((1,numberDicrestizations+1),dtype=object)
        anchor_distance = empty((1,numberDicrestizations+1),dtype=object)
        Ttop_tension = empty((1,numberDicrestizations+1),dtype=object)
        H_Force = empty((1,numberDicrestizations+1),dtype=object)
        FX = empty((1,numberDicrestizations+1),dtype=object)
        FY = empty((1,numberDicrestizations+1),dtype=object)
        sum_FX = array([0.]*(numberDicrestizations+1))
        stiffness = array([0.]*(numberDicrestizations+1))
        # fill arrays
        for i in range(0,numberDicrestizations+1):
            X_Fairlead[0,i] = array(X_Offset[i]+fairlead_x)
            anchor_distance[0,i] = array((((X_Fairlead[0,i])-anchor_x)**2+(fairlead_y-anchor_y)**2)**0.5-0.00001)
            Ttop_vect = array([0.]*numberMooringLines)
            H_vect = array([0.]*numberMooringLines)
            for j in range(0,numberMooringLines):
                Ttop_vect[j]=interp(anchor_distance[0,i][j],x,Ttop)/1000.
                H_vect[j]=interp(anchor_distance[0,i][j],x,H)/1000.
            Ttop_tension[0,i] = Ttop_vect
            H_Force[0,i] = H_vect
            FX[0,i] = array((anchor_x-X_Fairlead[0,i])/anchor_distance[0,i]*H_Force[0,i])
            FY[0,i] = array((anchor_y-fairlead_y)/anchor_distance[0,i]*H_Force[0,i])
            sum_FX[i] = sum(FX[0,i])
        for i in range(0,numberDicrestizations+1):    
            if i == numberDicrestizations:
                stiffness[i] = abs(sum_FX[i]/X_Offset[i])
            else: 
                stiffness[i] = abs((sum_FX[i]-sum_FX[i+1])/(X_Offset[i]-X_Offset[i+1]))    
        FR = power((power(FY,2)+power(FX,2)),0.5)
        # pack some things 
        fairlead_loc = [fairlead_x,fairlead_y,fairlead_z]
        anchor_loc = [anchor_x,anchor_y,anchor_z]


        massInAir = massInWater + waterDensity*pi*(mooringDiameter**2)/4
        for_MAP = InputMAP(waterDepth, g, waterDensity)
        for_MAP.write_line_dictionary_header()
        for_MAP.write_line_dictionary(mooringType, mooringDiameter, massInAir, AE_storm)
        for_MAP.write_node_properties_header()
        for_MAP.write_node_properties(1, "FIX", x0, anchor_y[0], anchor_z[0], 0, 0)
        for_MAP.write_node_properties(2, "VESSEL", fairlead_x[0], fairlead_y[0], fairlead_z[0], 0, 0) 
        for_MAP.write_line_properties_header()
        for_MAP.write_line_properties(1, mooringType, scope, 1, 2)
        for_MAP.write_solver_options(numberMooringLines)
        #execfile(os.path.abspath("../../src/main.py"))
        # read the stiffness_matrix.txt


        self.sum_forces_x = sum_FX
        self.offset_x = X_Offset
        # COST
        each_leg = costPerLength*scope
        legs_total = each_leg*numberMooringLines
        if self.anchor_type =='DRAG':
            each_anchor = minBreakingLoad/1000./9.806/20*2000
        elif self.anchor_type == 'PILE':
            each_anchor = 150000.*(minBreakingLoad/1000./9.806/1250.)**0.5
        if self.user_anchor_cost != 0.0: 
            each_anchor = self.user_anchor_cost
        anchor_total = each_anchor*numberMooringLines
        misc_cost = (anchor_total+legs_total)*self.misc_cost_factor/100.
        self.mooring_total_cost = legs_total+anchor_total+misc_cost 
        # INITIAL CONDITIONS
        keel2CG = DRAFT - fairleadDepth 
        mooringVerticalLoad =  interp(pretension,Ttop,Vtop)*numberMooringLines
        mooringHorizontalStiffness = interp(pretension,Ttop,mkh)
        mooringVerticalStiffness = interp(pretension,Ttop,mkv)*numberMooringLines
        TMM = (wetMassPerLength+pi*mooringDiameter**2/4*waterDensity)*scope*numberMooringLines
        self.mooring_keel_to_CG = keel2CG
        self.mooring_vertical_load = mooringVerticalLoad 
        self.mooring_horizontal_stiffness = mooringHorizontalStiffness
        self.mooring_vertical_stiffness = mooringVerticalStiffness
        self.mooring_mass = (wetMassPerLength+pi*mooringDiameter**2/4*waterDensity)*scope*numberMooringLines