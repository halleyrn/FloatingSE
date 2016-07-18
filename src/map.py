# This Python file uses the following encoding: utf-8

"""
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
"""

from mapsys import *
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
from math import cos, pi, sin
from numpy import array, set_printoptions, diag, interp
import os


class InputMAP(object):
    """The InputMAP class takes everything from FloatingSE and puts it in the
    correct format for MAP++. MAP++ then outputs a linearized stiffness matrix
    that FloatingSE uses in its optimization analysis."""
    def __init__(self, water_depth, gravity, water_density,
                 number_of_mooring_lines):
        """Initializes x, y, and z stiffness variables to 0. Initializes line
        type, fixed nodes, connect node, and vessel nodes as empty lists. Sets
        water depth, gravity, and water density as WATER_DEPTH, GRAVITY, and
        WATER_DENSITY, respectively."""
        super(InputMAP, self).__init__()
        # self.x_stiffness = 0
        # self.y_stiffness = 0
        # self.z_stiffness = 0
        self.line_type = 0
        self.water_depth = water_depth
        self.gravity = gravity
        self.water_density = water_density
        self.number_of_mooring_lines = number_of_mooring_lines
        self.MBL = 0
        self.WML = 0
        self.AE_storm = 0
        self.AE_drift = 0
        self.AREA = 0
        self.MCPL = 0
        self.diameter = 0
        self.sum_fx = []
        self.offset_x = []
        self.V_initial = 0
        self.horizontal_stiffness = 0
        self.vertical_stiffness = 0
        # self.damaged_mooring_bounds = [0, 0]
        # self.intact_mooring_bounds = [0, 0]
        self.max_tension = []

    def mooring_properties(self, mooring_diameter, line_type, mbl=0, wml=0,
                           ae_storm=0, mcpl=0):
        """Minimum breaking load (MBL), wet mass per length (WML), element axial
        stiffness (AE_storm and AE_drift), area, and MCPL are calculated here
        with MOORING DIAMETER and LINE TYPE"""
        self.diameter = mooring_diameter
        self.line_type = line_type
        if line_type == 'CHAIN':
            self.MBL = (27600.*(mooring_diameter**2)*(44.-80.*mooring_diameter))*10**3
            self.WML = 18070.*(mooring_diameter**2)
            self.AE_storm = (1.3788*mooring_diameter**2-4.93*mooring_diameter**3)*10**11
            self.AREA = 2.64*(mooring_diameter**2)
            self.MCPL = 0.58*self.MBL/1000./self.gravity-87.6
        elif line_type == 'STRAND':
            self.MBL = (937600*mooring_diameter**2-1408.3*mooring_diameter)*10**3
            self.WML = 4110*(mooring_diameter**2)
            self.AE_storm = 9.28*(mooring_diameter**2)*(10**10)
            self.AREA = 0.58*(mooring_diameter**2)
            self.MCPL = 0.42059603*(self.MBL/1000./self.gravity)+109.5
        elif line_type == 'IWRC':
            self.MBL = 648000*(mooring_diameter**2)*(10**3)
            self.WML = 3670*(mooring_diameter**2)
            self.AE_storm = 6.01*(mooring_diameter**2)*(10**10)
            self.AREA = 0.54*(mooring_diameter**2)
            self.MCPL = 0.33*(self.MBL/1000./self.gravity)+139.5
        elif line_type == 'FIBER': 
            self.MBL = (274700*(mooring_diameter**2)+7953.9*mooring_diameter-879.24)*(10**3)
            self.WML = 160.9*(mooring_diameter**2)+5.522*mooring_diameter-0.04798
            self.AE_storm = (10120*(mooring_diameter**2)+320.7*mooring_diameter-35.47)*(10**6)
            self.AE_drift = (5156*(mooring_diameter**2)+142.7*mooring_diameter-16)*(10**6)
            self.AREA = (pi/4)*(mooring_diameter**2)
            self.MCPL = 0.53676471*(self.MBL/1000./self.gravity)
        else: 
            print "PLEASE PICK AVAILABLE MOORING TYPE M8"
        if mbl != 0.0:
            self.MBL = mbl
        if wml != 0.0:
            self.WML = wml
        if ae_storm != 0.0:
            self.AE_storm = ae_storm
        if mcpl != 0.0:
            self.MCPL = mcpl

    def write_line_dictionary_header(self):
        """Writes the first three lines of the input.map file:
---------------------- LINE DICTIONARY ---------------------------------------
LineType  Diam      MassDenInAir   EA            CB   CIntDamp  Ca   Cdn    Cdt
(-)       (m)       (kg/m)        (N)           (-)   (Pa-s)    (-)  (-)    (-)
        """
        opened = open(os.path.abspath("../src/input.map"), "wb")
        opened.write("----------------------")
        opened.write(" LINE DICTIONARY ---------------------------------------\n")
        opened.write("LineType  ")
        opened.write("Diam      ")
        opened.write("MassDenInAir   ")
        opened.write("EA            ")
        opened.write("CB   ")
        opened.write("CIntDamp  ")
        opened.write("Ca   ")
        opened.write("Cdn    ")
        opened.write("Cdt\n")
        opened.write("(-)       ")
        opened.write("(m)       ")
        opened.write("(kg/m)        ")
        opened.write("(N)           ")
        opened.write("(-)   ")
        opened.write("(Pa-s)    ")
        opened.write("(-)  ")
        opened.write("(-)    ")
        opened.write("(-)\n")
        opened.close()

    def write_line_dictionary(self, air_mass_density, element_axial_stiffness, cable_sea_friction_coefficient=1):
        """Writes the forth line of the input.map file. This is where 
        LINE_TYPE, DIAMETER, AIR_MASS_DENSITY, ELEMENT_AXIAL_STIFFNESS, and
        CABLE_SEA_FRICTION_COEFFICIENT is inputted. 
        CABLE_SEA_FRICTION_COEFFICIENT defaults to 1 when none is given. If "#"
        is inputted into AIR_MASS_DENSITY and/or ELEMENT_AXIAL_STIFFNESS, then
        their respected calculated values are used."""
        opened = open(os.path.abspath("../src/input.map"), "ab")
        opened.write("%s   " % self.line_type)
        opened.write("%.5f   " % self.diameter)
        if air_mass_density == 0:
            air_mass_density = self.WML+(self.water_density*self.AREA)
        opened.write("%.5f   " % air_mass_density)
        if element_axial_stiffness == 0:
            element_axial_stiffness = self.AE_storm
        opened.write("%.5f   " % element_axial_stiffness)
        opened.write("%.5f   " % cable_sea_friction_coefficient)
        opened.write("1.0E8   ")
        opened.write("0.6   ")
        opened.write("-1.0   ")
        opened.write("0.05\n")
        opened.close()

    def write_node_properties_header(self):
        """Writes the node properties header:
---------------------- NODE PROPERTIES ---------------------------------------
Node  Type       X       Y       Z      M     B     FX      FY      FZ
(-)   (-)       (m)     (m)     (m)    (kg)  (mˆ3)  (N)     (N)     (N)
        """
        opened = open(os.path.abspath("../src/input.map"), "ab")
        opened.write("----------------------")
        opened.write(" NODE PROPERTIES ---------------------------------------\n")
        opened.write("Node  ")
        opened.write("Type       ")
        opened.write("X       ")
        opened.write("Y       ")
        opened.write("Z      ")
        opened.write("M     ")
        opened.write("B     ")
        opened.write("FX      ")
        opened.write("FY      ")
        opened.write("FZ\n")
        opened.write("(-)   ")
        opened.write("(-)       ")
        opened.write("(m)     ")
        opened.write("(m)     ")
        opened.write("(m)    ")
        opened.write("(kg)  ")
        opened.write("(mˆ3)  ")
        opened.write("(N)     ")
        opened.write("(N)     ")
        opened.write("(N)\n")
        opened.close()

    def write_node_properties(self, number, node_type, x_coordinate, y_coordinate, z_coordinate, point_mass_appl,
                              displaced_volume_appl, x_force_appl="#", y_force_appl="#", z_force_appl="#"):
        """Writes the input information for a node based on NODE_TYPE. X_FORCE_APPL, 
        Y_FORCE_APP, Z_FORCE_APP defaults to '#' if none is given."""
        opened = open(os.path.abspath("../src/input.map"), "ab")
        opened.write("%d   " % number)
        if node_type.lower() == "fix":
            opened.write("%s   " % node_type)
        elif node_type.lower() == "connect":
            opened.write("%s   " % node_type)
        elif node_type.lower() == "vessel":
            opened.write("%s   " % node_type)
        else:
            raise ValueError("%s is not a valid node type for node %d" % (node_type, number))
        if node_type.lower() == "connect":
            opened.write("#%.5f   " % x_coordinate)
            opened.write("#%.5f   " % y_coordinate)
            opened.write("#%.5f   " % z_coordinate)
            opened.write("%.5f   " % point_mass_appl)
            opened.write("%.5f   " % displaced_volume_appl)
            if not x_force_appl.isdigit() or not y_force_appl.isdigit() or not z_force_appl.isdigit():
                raise ValueError("%s must have numerical force applied values." % node_type)
            opened.write("%.5f   " % x_force_appl)
            opened.write("%.5f   " % y_force_appl)
            opened.write("%.5f\n" % z_force_appl)
        else:
            opened.write("%.5f   " % x_coordinate)
            opened.write("%.5f   " % y_coordinate)
            if z_coordinate == self.water_depth:
                opened.write("depth   ")
            else:
                opened.write("%.5f   " % z_coordinate)
            opened.write("%.5f   " % point_mass_appl)
            opened.write("%.5f   " % displaced_volume_appl)
            if str(x_force_appl).isdigit() or str(y_force_appl).isdigit() or str(z_force_appl).isdigit():
                raise ValueError("%s can only have '#' force applied values." % node_type)
            opened.write("%s   " % x_force_appl)
            opened.write("%s   " % y_force_appl)
            opened.write("%s\n" % z_force_appl)
        opened.close()

    def write_line_properties_header(self):
        """Writes the line properties header:
---------------------- LINE PROPERTIES ---------------------------------------
Line    LineType  UnstrLen  NodeAnch  NodeFair  Flags
(-)      (-)       (m)       (-)       (-)       (-)
        """
        opened = open(os.path.abspath("../src/input.map"), "ab")
        opened.write("----------------------")
        opened.write(" LINE PROPERTIES ---------------------------------------\n")
        opened.write("Line    ")
        opened.write("LineType  ")
        opened.write("UnstrLen  ")
        opened.write("NodeAnch  ")
        opened.write("NodeFair  ")
        opened.write("Flags\n")
        opened.write("(-)      ")
        opened.write("(-)       ")
        opened.write("(m)       ")
        opened.write("(-)       ")
        opened.write("(-)       ")
        opened.write("(-)\n")
        opened.close()

    def write_line_properties(self, line_number, line_type, unstretched_length, anchor_node_number,
                              fairlead_node_number, control_output_text_stream=" "):
        """Writes the input information for the line properties. This explains
        what node number is the ANCHOR and what node number is the FAIRLEAD, 
        as well as the UNSTRETCHED_LENGTH between the two nodes."""
        opened = open(os.path.abspath("../src/input.map"), "ab")
        opened.write("%d   " % line_number)
        opened.write("%s   " % line_type)
        # if unstretched_length == "#":
        #   unstretched_length = 
        opened.write("%.5f   " % unstretched_length)
        opened.write("%d   " % anchor_node_number)
        opened.write("%d   " % fairlead_node_number)
        opened.write("%s\n" % control_output_text_stream)
        opened.close()

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
        opened = open(os.path.abspath("../src/input.map"), "ab")
        opened.write("----------------------")
        opened.write(" SOLVER OPTIONS-----------------------------------------\n")
        opened.write("Option\n")
        opened.write("(-)\n")
        opened.write("help\n")
        opened.write(" integration_dt 0\n")
        opened.write(" kb_default 3.0e6\n")
        opened.write(" cb_default 3.0e5\n")
        opened.write(" wave_kinematics \n")
        opened.write("inner_ftol 1e-6\n")
        opened.write("inner_gtol 1e-6\n")
        opened.write("inner_xtol 1e-6\n")
        opened.write("outer_tol 1e-4\n")
        opened.write(" pg_cooked 10000 1\n")
        opened.write(" outer_fd \n")
        opened.write(" outer_bd \n")
        opened.write(" outer_cd\n")
        opened.write(" inner_max_its 100\n")
        opened.write(" outer_max_its 500\n")
        opened.write("repeat ")
        n = 360/self.number_of_mooring_lines
        degree = n
        while degree + n <= 360:
            opened.write("%d " % degree)
            degree += n
        opened.write("\n")
        opened.write(" krylov_accelerator 3\n")
        opened.write(" ref_position 0.0 0.0 0.0\n")

    def main(self, doffset, dangle, objective):
        """This runs MAP given the water DEPTH [m], GRAVITY[m/s^2], water DENSITY 
        [kg/m^3], the number of TOTAL_LINES, the MIN_BREAKING_LOAD [N], the DOFFSET
        [m], and DANGLE [degrees]. The vessel is displaced by DOFFSET until the
        maximum tension on at least one mooring line is over MIN_BREAKING_LOAD. Then
        the angle is changed by DANGLE and displaced by the offset once again. All
        tensions and the diagonals of the stiffness matrix are saved in a .txt."""

        float_formatter = lambda x: "%.3f" % x
        set_printoptions(precision=0)
        set_printoptions(suppress=True)
        set_printoptions(formatter={'float_kind': float_formatter})

        mooring_1 = Map()

        mooring_1.map_set_sea_depth(self.water_depth)
        mooring_1.map_set_gravity(self.gravity)
        mooring_1.map_set_sea_density(self.water_density)

        offset = doffset
        dangle = pi*dangle/180
        angle = 0
        intact_mooring = self.MBL*.60
        damaged_mooring = self.MBL*.80

        list_of_system_t = []
        tension = []

        mooring_1.read_file(os.path.abspath("../src/input.map"))  # 100 m depth
        mooring_1.init()

        epsilon = 1e-3
        K = mooring_1.linear(epsilon)    
        # print "\nHere is the linearized stiffness matrix with zero vessel displacement:"
        # print array(K)
        diagonal = array(diag(array(K)))
        self.horizontal_stiffness = diagonal.item(0)
        self.vertical_stiffness = diagonal.item(3)

        for line_number in range(0, self.number_of_mooring_lines):
            horizontal, vertical = mooring_1.get_fairlead_force_2d(line_number)
            tension.append((horizontal**2 + vertical**2)**.5)
            self.V_initial += vertical
            # print "Line %d: horizontal = %2.2f [N]  vertical = %2.2f [N] " \
            #       "tension = %2.2f [N]" % (line_number, horizontal, vertical, tension[line_number])
        
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

        if objective.lower() == "find full area" or objective == True:

            red_x = []
            red_y = []
            yellow_x = []
            yellow_y = []
            blue_x = []
            blue_y = []
            green_x = [0.]
            green_y = [0.]

            # finds the linearized stiffness matrix diagonal and tension in each line
            # as the vessel is displaced around 360 degrees
            while angle < 2*pi:
                max_tension = 0
                while max_tension <= self.MBL:
                    list_of_system_t.append(tension[:])
                    surge = offset*cos(angle)
                    sway = offset*sin(angle)
                    mooring_1.displace_vessel(surge, sway, 0, 0, 0, 0)
                    mooring_1.update_states(0.0, 0)
                 
                    # K = mooring_1.linear(epsilon)
                    # print "\nLinearized stiffness matrix with %2.2f surge and %2.2f sway vessel displacement:\n" % \
                    #       (surge, sway)
                    # print array(K)
                    for line_number in range(0, self.number_of_mooring_lines):
                        # fx ,fy ,fz = mooring_1.get_fairlead_force_3d(line_number)
                        horizontal, vertical = mooring_1.get_fairlead_force_2d(line_number)
                        tension[line_number] = (horizontal**2 + vertical**2)**.5
                        # print "Line %d: horizontal = %2.2f [N]  vertical = %2.2f [N] " \
                        #       "tension = %2.2f [N]" % (line_number, horizontal, vertical, tension[line_number])
                    offset += doffset
                    max_tension = max(tension)
                    if max_tension >= self.MBL:
                        red_x.append(surge)
                        red_y.append(sway)
                    elif max_tension >= damaged_mooring:
                        yellow_x.append(surge)
                        yellow_y.append(sway)
                    elif max_tension >= intact_mooring:
                        blue_x.append(surge)
                        blue_y.append(sway)
                    else:
                        green_x.append(surge)
                        green_y.append(sway)
                angle += dangle
                offset = doffset

            # uncomment if you want to see offsets plotted
            plt.plot(red_x, red_y, 'ro', yellow_x, yellow_y, 'yo', blue_x, blue_y, 'bo', green_x, green_y, 'go')
            plt.axis([-60, 80, -70, 70])
            plt.show()

        if objective.lower() == "optimization" or objective == True:
            # find sum of fx at each displacement along the x-axis (list from neg to
            # pos offset) and the offset in x as well
            max_tension = 0 
            surge = 0
            while max_tension <= self.MBL:
                tot_fx = 0
                mooring_1.displace_vessel(surge, 0, 0, 0, 0, 0)
                mooring_1.update_states(0.0, 0)
                # K = mooring_1.linear(epsilon)
                for line_number in range(0, self.number_of_mooring_lines):
                    fx, fy, fz = mooring_1.get_fairlead_force_3d(line_number)
                    tot_fx += -fx
                    horizontal, vertical = mooring_1.get_fairlead_force_2d(line_number)
                    tension[line_number] = (horizontal**2 + vertical**2)**.5
                max_tension = max(tension[:])
                self.max_tension.append(max_tension)
                self.sum_fx.append(tot_fx)
                self.offset_x.append(surge)
                surge += doffset
            max_tension = 0
            surge = -doffset
            while max_tension <= self.MBL:
                tot_fx = 0
                mooring_1.displace_vessel(surge, 0, 0, 0, 0, 0)
                mooring_1.update_states(0.0, 0)
                # K = mooring_1.linear(epsilon)
                for line_number in range(0, self.number_of_mooring_lines):
                    fx, fy, fz = mooring_1.get_fairlead_force_3d(line_number)
                    tot_fx += -fx
                    horizontal, vertical = mooring_1.get_fairlead_force_2d(line_number)
                    tension[line_number] = (horizontal**2 + vertical**2)**.5
                max_tension = max(tension[:])
                self.max_tension.insert(0, max_tension)
                self.sum_fx.insert(0, tot_fx)
                self.offset_x.insert(0, surge)
                surge -= doffset
            x = self.offset_x.index(0)
            positive_x = self.offset_x[x:]
            negative_x = list(reversed(self.offset_x[:x]))
            positive_x_max_t = self.max_tension[x:]
            negative_x_max_t = list(reversed(self.max_tension[:x]))
            max_x = interp(self.MBL, positive_x_max_t, positive_x)
            mooring_1.displace_vessel(max_x, 0, 0, 0, 0, 0)
            mooring_1.update_states(0.0, 0)
            tot_fx = 0
            for line_number in range(0, self.number_of_mooring_lines):
                fx, fy, fz = mooring_1.get_fairlead_force_3d(line_number)
                tot_fx += -fx
                horizontal, vertical = mooring_1.get_fairlead_force_2d(line_number)
                tension[line_number] = (horizontal ** 2 + vertical ** 2) ** .5
            max_tension = max(tension[:])
            self.max_tension[-1] = max_tension
            self.sum_fx[-1] = tot_fx
            self.offset_x[-1] = max_x
            min_x = interp(self.MBL, negative_x_max_t, negative_x)
            mooring_1.displace_vessel(min_x, 0, 0, 0, 0, 0)
            mooring_1.update_states(0.0, 0)
            tot_fx = 0
            for line_number in range(0, self.number_of_mooring_lines):
                fx, fy, fz = mooring_1.get_fairlead_force_3d(line_number)
                tot_fx += -fx
                horizontal, vertical = mooring_1.get_fairlead_force_2d(line_number)
                tension[line_number] = (horizontal ** 2 + vertical ** 2) ** .5
            max_tension = max(tension[:])
            self.max_tension[0] = max_tension
            self.sum_fx[0] = tot_fx
            self.offset_x[0] = min_x
        # print "angle changed: %d offset changed: %d" %(angle_changed, offset_changed) #delete
        # opened.write(str(list_of_system_t) + "\n" )
        # opened.close

        mooring_1.end()

    def sum_of_fx_and_offset(self):
        return array(self.sum_fx), array(self.offset_x)

    def wet_mass_per_length(self):
        return self.WML

    def cost_per_length(self):
        return self.MCPL

    def minimum_breaking_load(self):
        return self.MBL

    def loads_and_stiffnesses(self):
        return self.V_initial, self.vertical_stiffness, self.horizontal_stiffness

    def intact_and_damaged_mooring(self):
        x = self.offset_x.index(0)
        positive_x = self.offset_x[x:]
        negative_x = list(reversed(self.offset_x[:x]))
        positive_x_max_t = self.max_tension[x:]
        negative_x_max_t = list(reversed(self.max_tension[:x]))
        intact_mooring = float(self.MBL*.60)
        damaged_mooring = float(self.MBL*.80)
        intact_mooring_bounds = list([interp(intact_mooring, negative_x_max_t, negative_x)])
        intact_mooring_bounds.append(interp(intact_mooring, positive_x_max_t, positive_x))
        damaged_mooring_bounds = list([interp(damaged_mooring, negative_x_max_t, negative_x)])
        damaged_mooring_bounds.append(interp(damaged_mooring, positive_x_max_t, positive_x))
        return intact_mooring_bounds, damaged_mooring_bounds

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
    OC3.main(2, 2, "optimization")
    intact_mooring1, damaged_mooring1 = OC3.intact_and_damaged_mooring()
    sum_fx1, offset_x1 = OC3.sum_of_fx_and_offset()
    print sum_fx1
    print offset_x1
    print OC3.max_tension
    print intact_mooring1, damaged_mooring1
