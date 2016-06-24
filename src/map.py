# This Python file uses the following encoding: utf-8

'''
  Copyright (C) 2014 mdm                                     
  marco[dot]masciola[at]gmail                                
                                                             
Licensed to the Apache Software Foundation (ASF) under one   
or more contributor license agreements.  See the NOTICE file 
distributed with this work for additional information        
regarding copyright ownership.  The ASF licenses this file   
to you under the Apache License, Version 2.0 (the            
"License"); you may not use this file except in compliance   
with the License.  You may obtain a copy of the License at   
                                                             
  http://www.apache.org/licenses/LICENSE-2.0                 
                                                             
Unless required by applicable law or agreed to in writing,   
software distributed under the License is distributed on an  
"AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY       
KIND, either express or implied.  See the License for the    
specific language governing permissions and limitations            
under the License.                                             
'''  

from mapsys import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from numpy import array, set_printoptions, diag, cos, pi, sin
import os
import sys

class InputMAP(object):
    """The InputMAP class takes everything from FloatingSE and puts it in the
    correct format for MAP++. MAP++ then outputs a linearlized stiffness matrix
    that FloatingSE uses in its optimization analysis."""
    def __init__(self, water_depth, gravity, water_density, 
        number_of_mooring_lines):
        """Initalizes x, y, and z stiffness variables to 0. Initalizes line
        type, fixed nodes, connect node, and vessel nodes as empyt lists. Sets
        water depth, gravity, and water density as WATER_DEPTH, GRAVITY, and
        WATER_DENSITY, respectively."""
        super(InputMAP, self).__init__()
        self.x_stiffness = 0
        self.y_stiffness = 0
        self.z_stiffness = 0
        self.line_type = 0
        self.water_depth = water_depth
        self.gravity = gravity
        self.water_density = water_density
        self.number_of_mooring_lines = number_of_mooring_lines
        self.MBL = 0
        self.WML = 0
        self.AE_storm = 0
        self.AREA = 0
        self.MCPL = 0
        self.diameter = 0

    def mooring_properties(self, mooringDiameter, line_type):
        """Minimun breaking load (MBL), wet mass per length (WML), element axial
        stiffness (AE_storm and AE_drift), area, and MCPL are calcuated here
        with MOORING DIAMETER and LINE TYPE"""
        self.diameter = mooringDiameter
        self.line_type = line_type
        if line_type == 'CHAIN':
            # print "I am a chain"
            self.MBL = (27600.*(mooringDiameter**2)*(44.-80.*mooringDiameter))*10**3
            self.WML = 18070.*(mooringDiameter**2)
            self.AE_storm = (1.3788*mooringDiameter**2-4.93*mooringDiameter**3)*10**11
            self.AREA = 2.64*(mooringDiameter**2)
            self.MCPL = 0.58*self.MBL/1000./self.gravity-87.6
            # print "MBL: %d WML: %d AE_storm: %d area: %.4f MCPL: %d" % (self.MBL, self.WML, self.AE_storm, self.AREA, self.MCPL)
        elif line_type == 'STRAND':
            self.MBL = (937600*mooringDiameter**2-1408.3*mooringDiameter)*10**3
            self.WML = 4110*(mooringDiameter**2)
            self.AE_storm = 9.28*(mooringDiameter**2)*(10**10)
            self.AREA = 0.58*(mooringDiameter**2)
            self.MCPL = 0.42059603*(self.MBL/1000./self.gravity)+109.5
        elif line_type == 'IWRC':
            self.MBL = 648000*(mooringDiameter**2)*(10**3)
            self.WML = 3670*(mooringDiameter**2)
            self.AE_storm = 6.01*(mooringDiameter**2)*(10**10)
            self.AREA = 0.54*(mooringDiameter**2)
            self.MCPL = 0.33*(self.MBL/1000./self.gravity)+139.5
        elif line_type == 'FIBER': 
            self.MBL = (274700*(mooringDiameter**2)+7953.9*mooringDiameter-879.24)*(10**3)
            self.WML = 160.9*(mooringDiameter**2)+5.522*mooringDiameter-0.04798
            self.AE_storm = (10120*(mooringDiameter**2)+320.7*mooringDiameter-35.47)*(10**6)
            self.AE_drift = (5156*(mooringDiameter**2)+142.7*mooringDiameter-16)*(10**6)
            self.AREA = (pi/4)*(mooringDiameter**2)
            self.MCPL = 0.53676471*(self.MBL/1000./self.gravity)
        else: 
            print "PLEASE PICK AVAILABLE MOORIN' TYPE M8"

    def write_line_dictionary_header(self):
        """Writes the first three lines of the input.map file:
---------------------- LINE DICTIONARY ---------------------------------------
LineType  Diam      MassDenInAir   EA            CB   CIntDamp  Ca   Cdn    Cdt
(-)       (m)       (kg/m)        (N)           (-)   (Pa-s)    (-)  (-)    (-)
        """
        file = open(os.path.abspath("../src/input.map"), "wb")
        file.write("----------------------")
        file.write(" LINE DICTIONARY ---------------------------------------\n")
        file.write("LineType  ")
        file.write("Diam      ")
        file.write("MassDenInAir   ")
        file.write("EA            ")
        file.write("CB   ")
        file.write("CIntDamp  ")
        file.write("Ca   ")
        file.write("Cdn    ")
        file.write("Cdt\n")
        file.write("(-)       ")
        file.write("(m)       ")
        file.write("(kg/m)        ")
        file.write("(N)           ")
        file.write("(-)   ")
        file.write("(Pa-s)    ")
        file.write("(-)  ")
        file.write("(-)    ")
        file.write("(-)\n")
        file.close()

    def write_line_dictionary(self, air_mass_density = "#", 
        element_axial_stiffness = "#", cable_sea_friction_coefficient=1):
        """Writes the forth line of the input.map file. This is where 
        LINE_TYPE, DIAMETER, AIR_MASS_DENSITY, ELEMENT_AXIAL_STIFFNESS, and
        CABLE_SEA_FRICTION_COEFFICIENT is inputted. 
        CABLE_SEA_FRICTION_COEFFICIENT defaults to 1 when none is given. If "#"
        is inputted into AIR_MASS_DENSITY and/or ELEMENT_AXIAL_STIFFNESS, then
        their respected calculated values are used."""
        file = open(os.path.abspath("../src/input.map"), "ab")
        file.write("%s   " % self.line_type)
        file.write("%.5f   " % self.diameter)
        if air_mass_density == "#":
            air_mass_density = self.WML+(self.water_density*self.AREA)
        file.write("%.5f   " % air_mass_density)
        if element_axial_stiffness == "#":
            element_axial_stiffness = self.AE_storm
        file.write("%.5f   " % element_axial_stiffness)
        file.write("%.5f   " % cable_sea_friction_coefficient)
        file.write("1.0E8   ")
        file.write("0.6   ")
        file.write("-1.0   ")
        file.write("0.05\n")
        file.close()

    def write_node_properties_header(self):
        """Writes the node properties header:
---------------------- NODE PROPERTIES ---------------------------------------
Node  Type       X       Y       Z      M     B     FX      FY      FZ
(-)   (-)       (m)     (m)     (m)    (kg)  (mˆ3)  (N)     (N)     (N)
        """
        file = open(os.path.abspath("../src/input.map"), "ab")
        file.write("----------------------")
        file.write(" NODE PROPERTIES ---------------------------------------\n")
        file.write("Node  ")
        file.write("Type       ")
        file.write("X       ")
        file.write("Y       ")
        file.write("Z      ")
        file.write("M     ")
        file.write("B     ")
        file.write("FX      ")
        file.write("FY      ")
        file.write("FZ\n")
        file.write("(-)   ")
        file.write("(-)       ")
        file.write("(m)     ")
        file.write("(m)     ")
        file.write("(m)    ")
        file.write("(kg)  ")
        file.write("(mˆ3)  ")
        file.write("(N)     ")
        file.write("(N)     ")
        file.write("(N)\n")
        file.close()


    def write_node_properties(self, number, node_type, x_coordinate, y_coordinate,
        z_coordinate, point_mass_appl, displaced_volume_appl, x_force_appl ="#",
        y_force_appl = "#", z_force_appl = "#"):
        """Writes the input information for a node based on NODE_TYPE. X_FORCE_APPL, 
        Y_FORCE_APP, Z_FORCE_APP defaults to '#' if none is given."""
        file = open(os.path.abspath("../src/input.map"), "ab")
        file.write("%d   " % number)
        if node_type.lower() == "fix":
            file.write("%s   " % node_type)
        elif node_type.lower() == "connect":
            file.write("%s   " % node_type)
        elif node_type.lower() == "vessel":
            file.write("%s   " % node_type)
        else:
            raise ValueError("%s is not a valid node type for node %d" 
                % (node_type, number))
        
        if node_type.lower() == "connect":
            file.write("#%.5f   " % x_coordinate)
            file.write("#%.5f   " % y_coordinate)
            file.write("#%.5f   " % z_coordinate)
            file.write("%.5f   " % point_mass_appl)
            file.write("%.5f   " % displaced_volume_appl)
            if not x_force_appl.isdigit() or not y_force_appl.isdigit() or not z_force_appl.isdigit():
                raise ValueError("%s must have numerical force applied values."
                    % node_type)
            file.write("%.5f   " % x_force_appl)
            file.write("%.5f   " % y_force_appl)
            file.write("%.5f\n" % z_force_appl)
        else:
            file.write("%.5f   " % x_coordinate)
            file.write("%.5f   " % y_coordinate)
            if z_coordinate == self.water_depth:
                file.write("depth   ")
            else:
                file.write("%.5f   " % z_coordinate)
            file.write("%.5f   " % point_mass_appl)
            file.write("%.5f   " % displaced_volume_appl)
            if str(x_force_appl).isdigit() or str(y_force_appl).isdigit() or str(z_force_appl).isdigit():
                raise ValueError("%s can only have '#' force applied values."
                    % node_type)
            file.write("%s   " % x_force_appl)
            file.write("%s   " % y_force_appl)
            file.write("%s\n" % z_force_appl)
        file.close()


    def write_line_properties_header(self):
        """Writes the line properties header:
---------------------- LINE PROPERTIES ---------------------------------------
Line    LineType  UnstrLen  NodeAnch  NodeFair  Flags
(-)      (-)       (m)       (-)       (-)       (-)
        """
        file = open(os.path.abspath("../src/input.map"), "ab")
        file.write("----------------------")
        file.write(" LINE PROPERTIES ---------------------------------------\n")
        file.write("Line    ")
        file.write("LineType  ")
        file.write("UnstrLen  ")
        file.write("NodeAnch  ")
        file.write("NodeFair  ")
        file.write("Flags\n")
        file.write("(-)      ")
        file.write("(-)       ")
        file.write("(m)       ")
        file.write("(-)       ")
        file.write("(-)       ")
        file.write("(-)\n")
        file.close()

    def write_line_properties(self, line_number, line_type, unstretched_length,
        anchor_node_number, fairlead_node_number, control_output_text_stream = " "):
        """Writes the input information for the line properties. This explains
        what node number is the ANCHOR and what node number is the FAIRLEAD, 
        as well as the UNSTRETCHED_LENGTH between the two nodes."""
        file = open(os.path.abspath("../src/input.map"), "ab")
        file.write("%d   " % line_number)
        file.write("%s   " % line_type)
        # if unstretched_length == "#":
        #   unstretched_length = 
        file.write("%.5f   " % unstretched_length)
        file.write("%d   " % anchor_node_number)
        file.write("%d   " % fairlead_node_number)
        file.write("%s\n" % control_output_text_stream)
        file.close()

    def write_solver_options(self):
        """Writes the solver options at the end of the input file, as well as 
        takes the self.NUMBER_OF_MOORING_LINES and places them evenly within 360
        degrees. For NUMBER_OF_MOORING_LINES = 3:
---------------------- SOLVER OPTIONS-----------------------------------------
Option
(-)
help
 integration_dt 0
 kb_default 3.0e6
 cb_default 3.0e5
 wave_kinematics 
inner_ftol 1e-6
inner_gtol 1e-6
inner_xtol 1e-6
outer_tol 1e-4
 pg_cooked 10000 1
 outer_fd 
 outer_bd 
 outer_cd
 inner_max_its 100
 outer_max_its 500
repeat 120 240 
 krylov_accelerator 3
 ref_position 0.0 0.0 0.0
        """
        file = open(os.path.abspath("../src/input.map"), "ab")
        file.write("----------------------")
        file.write(" SOLVER OPTIONS-----------------------------------------\n")
        file.write("Option\n")
        file.write("(-)\n")
        file.write("help\n")
        file.write(" integration_dt 0\n")
        file.write(" kb_default 3.0e6\n")
        file.write(" cb_default 3.0e5\n")
        file.write(" wave_kinematics \n")
        file.write("inner_ftol 1e-6\n")
        file.write("inner_gtol 1e-6\n")
        file.write("inner_xtol 1e-6\n")
        file.write("outer_tol 1e-4\n")
        file.write(" pg_cooked 10000 1\n")
        file.write(" outer_fd \n")
        file.write(" outer_bd \n")
        file.write(" outer_cd\n")
        file.write(" inner_max_its 100\n")
        file.write(" outer_max_its 500\n")
        file.write("repeat ")
        n = 360/self.number_of_mooring_lines
        degree = n
        while degree + n <= 360:
            file.write("%d " % degree)
            degree += n
        file.write("\n")
        file.write(" krylov_accelerator 3\n")
        file.write(" ref_position 0.0 0.0 0.0\n")


def main(depth, gravity, density, total_lines, min_breaking_load, doffset, dangle):
    """This runs MAP given the water DEPTH [m], GRAVITY[m/s^2], water DENSITY 
    [kg/m^3], the number of TOTAL_LINES, the MIN_BREAKING_LOAD [N], the DOFFSET
    [m], and DANGLE [degrees]. The vessel is displaced by DOFFSET until the
    maxium tension on at least one mooring line is over MIN_BREAKING_LOAD. Then
    the angle is changed by DANGLE and displaced by the offset once again. All
    tensions and the diagonals of the stiffness maxtrix are saved in a .txt.""" 

    float_formatter = lambda x: "%.3f" % x
    set_printoptions(precision=0)
    set_printoptions(suppress=True)
    set_printoptions(formatter={'float_kind':float_formatter})

    mooring_1 = Map( )

    mooring_1.map_set_sea_depth(depth)
    mooring_1.map_set_gravity(gravity)
    mooring_1.map_set_sea_density(density)


    offset = doffset
    dangle = pi*dangle/180
    angle = 0

    list_of_system_T = []
    T = []

    mooring_1.read_file(os.path.abspath("../src/input.map")) # 100 m depth
    mooring_1.init( )

    epsilon = 1e-3
    K = mooring_1.linear(epsilon)    
    print "\nHere is the linearized stiffness matrix with zero vessel displacement:"
    print array(K)
    file = open(os.path.abspath("../src/stiffness_diagonals.txt"),"w")
    file.write(str(diag(array(K))) + "\n")
    file.close

    for line_number in range(0, total_lines):
        H,V = mooring_1.get_fairlead_force_2d(line_number)
        T.append((H**2 + V**2)**.5)
        # print "Line %d: H = %2.2f [N]  V = %2.2f [N] T = %2.2f [N]"%(line_number, H, V, T[line_number])     
    
    # fig = plt.figure()
    # ax = Axes3D(fig)
    # for i in range(0,mooring_1.size_lines()):
    #     x = mooring_1.plot_x( i, 10 )
    #     y = mooring_1.plot_y( i, 10 )
    #     z = mooring_1.plot_z( i, 10 )        
    #     ax.plot(x,y,z,'b-')
     
    # ax.set_xlabel('X [m]')
    # ax.set_ylabel('Y [m]')
    # ax.set_zlabel('Z [m]')
     
    # plt.show()
    file = open(os.path.abspath("../src/stiffness_diagonals.txt"),"a")
    # print "lines: %s MBL: %s doffset: %s dangle: %s angle: %s" % (total_lines, min_breaking_load, doffset, dangle, angle) #delete
    angle_changed = 0 #delete
    offset_changed = 0 #delete

    intact_mooring = min_breaking_load*.60
    damaged_mooring = min_breaking_load*.80
    red_x = []
    red_y = []
    yellow_x = []
    yellow_y = []
    blue_x = []
    blue_y = []
    green_x = [0]
    green_y = [0]

    while angle < 2*pi:
        max_tension = 0
        while max_tension < min_breaking_load:
            list_of_system_T.append(T[:])
            surge = offset*cos(angle)
            sway = offset*sin(angle)
            mooring_1.displace_vessel(surge,sway,0,0,0,0)
            mooring_1.update_states(0.0,0)
         
            K = mooring_1.linear(epsilon)    
            # print "\nLinearized stiffness matrix with %2.2f surge and %2.2f sway vessel displacement:\n"%(surge, sway)
            # print array(K)
            file.write(str(diag(array(K))) + "\n" )
            for line_number in range(0, total_lines):
                H,V = mooring_1.get_fairlead_force_2d(line_number)
                T[line_number] = (H**2 + V**2)**.5    
                # print "Line %d: H = %2.2f [N]  V = %2.2f [N] T = %2.2f [N]"%(line_number, H, V, T[line_number]) 
            offset += doffset
            offset_changed +=1 #delete
            max_tension = max(T)   
            if max_tension > min_breaking_load:
                red_x.append(surge)
                red_y.append(sway)
            elif max_tension > damaged_mooring:
                yellow_x.append(surge)
                yellow_y.append(sway)
            elif max_tension > intact_mooring:
                blue_x.append(surge)
                blue_y.append(sway)
            else:
                green_x.append(surge)
                green_y.append(sway)
        angle += dangle
        offset = doffset
        angle_changed += 1 #delete
    print "angle changed: %d offset changed: %d" %(angle_changed, offset_changed) #delete
    file.write(str(list_of_system_T) + "\n" )
    file.close

    plt.plot(red_x, red_y, 'ro', yellow_x, yellow_y, 'yo', blue_x, blue_y, 'bo', green_x, green_y, 'go')
    plt.axis([-60, 80, -70, 70])
    plt.show()

    mooring_1.end( )

if __name__ == '__main__':
    """Testing the interface using homogeneous line OC3 mooring information."""
    OC3 = InputMAP(320.0, 9.806, 1025.0, 3)
    OC3.mooring_properties(0.09, "CHAIN")
    OC3.write_line_dictionary_header()
    OC3.write_line_dictionary(77.7066, 384243000)
    OC3.write_node_properties_header()
    OC3.write_node_properties(1, "FIX", 853.87, 0, 320.0, 0, 0)
    OC3.write_node_properties(2, "VESSEL", 5.2, 0, -70.0, 0, 0)
    OC3.write_line_properties_header()
    OC3.write_line_properties(1, "CHAIN", 902.2, 1, 2, " ")
    OC3.write_solver_options()
    depth = OC3.water_depth
    gravity = OC3.gravity
    density = OC3.water_density
    line = OC3.number_of_mooring_lines
    MBL = OC3.MBL
    main(depth, gravity, density, line, MBL, 2, 2)
