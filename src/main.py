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
import numpy as np

# def run(sea_depth, gravity, sea_density):
#     """This function is called within mooring.py in order to run MAP++. It takes
#     in SEA_DEPTH, GRAVITY, and SEA_DENSITY that is used in the MAP++ program."""
#     np.set_printoptions(precision=0)
#     np.set_printoptions(suppress=True)

#     mooring_1 = Map( )

#     mooring_1.map_set_sea_depth(sea_depth)
#     mooring_1.map_set_gravity(gravity)
#     mooring_1.map_set_sea_density(sea_density)

#     mooring_1.read_file("./input.map") # 100 m depth
#     mooring_1.init( )

#     epsilon = 1e-3
#     K = mooring_1.linear(epsilon)    
#     print "\nHere is the linearized stiffness matrix with zero vessel displacement:"
#     print np.array(K)
#     file = open("./stiffness_matrix.txt","w")
#     file.write(str(np.array(K)))
#     file.close
    
#     fig = plt.figure()
#     ax = Axes3D(fig)
#     for i in range(0,mooring_1.size_lines()):
#         x = mooring_1.plot_x( i, 10 )
#         y = mooring_1.plot_y( i, 10 )
#         z = mooring_1.plot_z( i, 10 )        
#         ax.plot(x,y,z,'b-')
     
#     ax.set_xlabel('X [m]')
#     ax.set_ylabel('Y [m]')
#     ax.set_zlabel('Z [m]')           
     
#     plt.show()

#     mooring_1.end( )


if __name__ == '__main__':
    """This runs if the file is run directly through terminal or command prompt."""      
    from mapsys import *
    import matplotlib.pyplot as plt
    from mpl_toolkits.mplot3d import Axes3D
    from numpy import array, set_printoptions
    import os

    file = open(os.path.abspath("../src/main_input.txt"),"r")
    data = file.readline().split()
    file.close()


    set_printoptions(precision=0)
    set_printoptions(suppress=True)

    mooring_1 = Map( )
    
    print data[1]
    print data[3]
    print data[5]
    mooring_1.map_set_sea_depth(float(data[1]))
    mooring_1.map_set_gravity(float(data[3]))
    mooring_1.map_set_sea_density(float(data[5]))

    mooring_1.read_file(os.path.abspath("../src/input.map")) # 100 m depth
    mooring_1.init( )

    epsilon = 1e-3
    K = mooring_1.linear(epsilon)    
    print "\nHere is the linearized stiffness matrix with zero vessel displacement:"
    print array(K)
    file = open(os.path.abspath("../src/stiffness_matrix.txt"),"w")
    file.write(str(array(K)))
    file.close
    
    fig = plt.figure()
    ax = Axes3D(fig)
    for i in range(0,mooring_1.size_lines()):
        x = mooring_1.plot_x( i, 10 )
        y = mooring_1.plot_y( i, 10 )
        z = mooring_1.plot_z( i, 10 )        
        ax.plot(x,y,z,'b-')
     
    ax.set_xlabel('X [m]')
    ax.set_ylabel('Y [m]')
    ax.set_zlabel('Z [m]')           
     
    plt.show()
    
    surge = 5.0 # 5 meter surge displacements
    mooring_1.displace_vessel(surge,0,0,0,0,0)
    mooring_1.update_states(0.0,0)
     
    K = mooring_1.linear(epsilon)    
    print "\nLinearized stiffness matrix with %2.2f surge vessel displacement:\n"%(surge)
    print np.array(K)

    fig = plt.figure()
    ax = Axes3D(fig)
    for i in range(0,mooring_1.size_lines()):
        x = mooring_1.plot_x( i, 10 )
        y = mooring_1.plot_y( i, 10 )
        z = mooring_1.plot_z( i, 10 )        
        ax.plot(x,y,z,'b-')
     
    ax.set_xlabel('X [m]')
    ax.set_ylabel('Y [m]')
    ax.set_zlabel('Z [m]')           
     
    plt.show()
    
    mooring_1.end( )
