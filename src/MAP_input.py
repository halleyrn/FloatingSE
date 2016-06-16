# This Python file uses the following encoding: utf-8
import os, sys



class MainFile(object):
	"""docstring for MainFile"""
	def __init__(self):
		super(MainFile, self).__init__()
		self.x_stiffness = 0
		self.y_stiffness = 0
		self.z_stiffness = 0
		self.water_depth = 0
		self.line_types = []
		self.fix_nodes = []
		self.connect_nodes = []
		self.vessel_nodes = []

	def write_to_main_input(self, water_depth, gravity, water_density):
		self.water_depth = water_depth
		#print "water depth: "+ str(self.water_depth)
		file = open("./MAP/main_input.txt", "wb")
		file.write("water_depth:\t%s\t" % str(water_depth))
		file.write("gravity:\t%s\t" % str(gravity))
		file.write("water_density:\t%s\n" % str(water_density))
		file.close()

	def write_line_dictionary_header(self):
		"""Writes the first three lines of the input.map file."""
		file = open("./MAP/input.map", "wb")
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

	def write_line_dictionary(self, line_type, diameter, air_mass_density,
		element_axial_stiffness, cable_sea_friction_coefficient):
		"""Writes the forth line of the input.map file. This is where line type
		properties are inputted."""
		file = open("./MAP/input.map", "ab")
		if line_type in self.line_types:
			raise ValueError("Already have %s line type." % line_type)
		else:
			file.write("%s   " % line_type)
			self.line_types.append(line_type)
		file.write("%f   " % diameter)
		file.write("%f   " % air_mass_density)
		file.write("%f   " % element_axial_stiffness)
		file.write("%f   " % cable_sea_friction_coefficient)
		file.write("1.0E8   ")
		file.write("0.6   ")
		file.write("-1.0   ")
		file.write("0.05\n")
		file.close()

	def write_node_properties_header(self):
		"""Writes the node properties header"""
		file = open("./MAP/input.map", "ab")
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
		z_coordinate, point_mass_appl, displaced_volume_appl, x_force_appl,
		y_force_appl, z_force_appl):
		"""Writes the input information for all the nodes."""
		file = open("./MAP/input.map", "ab")
		file.write("%d   " % number)
		if node_type.lower() == "fix":
			file.write("%s   " % node_type)
			self.fix_nodes.append(number)
		elif node_type.lower() == "connect":
			file.write("%s   " % node_type)
			self.connect_nodes.append(number)
		elif node_type.lower() == "vessel":
			file.write("%s   " % node_type)
			self.vessel_nodes.append(number)
		else:
			raise ValueError("%s is not a valid node type for node %d" 
				% (node_type, number))
		
		if node_type.lower() == "connect":
			file.write("#%f   " % x_coordinate)
			file.write("#%f   " % y_coordinate)
			file.write("#%f   " % z_coordinate)
			file.write("%f   " % point_mass_appl)
			file.write("%f   " % displaced_volume_appl)
			if not x_force_appl.isdigit() or not y_force_appl.isdigit() or not z_force_appl.isdigit():
				raise ValueError("%s must have numerical force applied values."
					% node_type)
			file.write("%f   " % x_force_appl)
			file.write("%f   " % y_force_appl)
			file.write("%f\n" % z_force_appl)
		else:
			file.write("%f   " % x_coordinate)
			file.write("%f   " % y_coordinate)
			if z_coordinate == self.water_depth:
				file.write("depth   ")
			else:
				file.write("%f   " % z_coordinate)
			file.write("%f   " % point_mass_appl)
			file.write("%f   " % displaced_volume_appl)
			if x_force_appl.isdigit() or y_force_appl.isdigit() or z_force_appl.isdigit():
				raise ValueError("%s can only have '#' force applied values."
					% node_type)
			file.write("%s   " % x_force_appl)
			file.write("%s   " % y_force_appl)
			file.write("%s\n" % z_force_appl)
		file.close()


	def write_line_properties_header(self):
		file = open("./MAP/input.map", "ab")
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
		anchor_node_number, fairlead_node_number, control_output_text_stream):
		file = open("./MAP/input.map", "ab")
		file.write("%d   " % line_number)
		if line_type in self.line_types:
			file.write("%s   " % line_type)
			
		else:
			raise ValueError("%s is not a previously defined linetype within the" 
				+ " Line Dictionary." % line_type)
		file.write("%f   " % unstretched_length)
		# if anchor_node_number in self.connect_nodes or anchor_node_number in self.fix_nodes:
		# 	if fairlead_node_number in self.connect_nodes or fairlead_node_number in self.vessel_nodes:
		# 		file.write("%d   " % anchor_node_number)
		# 		file.write("%d   " % fairlead_node_number)
		# 	else:
		# 		raise ValueError("%d cannot be an fairlead node" % fairlead_node_number)
		# else:
		# 	raise ValueError("%d cannot be an anchor node" % anchor_node_number)
		file.write("%d   " % anchor_node_number)
		file.write("%d   " % fairlead_node_number)
		file.write("%s\n" % control_output_text_stream)
		file.close()

	def write_solver_options(self):
		file = open("./MAP/input.map", "ab")
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
		file.write("repeat 120 240\n")
		file.write(" krylov_accelerator 3\n")
		file.write(" ref_position 0.0 0.0 -6.0\n")
		file.write(" -172.772029\n")
		file.write("    0.000000\n")
		file.write("  -82.007818\n")
		file.write("  -63.224371\n")
		file.write("    0.000000\n")
		file.write("  -28.179102\n")
		file.close()

if __name__ == '__main__':
	"""Testing the interface using homogeneous line OC3 mooring information."""
	OC3 = MainFile()
	OC3.write_to_main_input(320.0, 9.806, 1025.0)
	OC3.write_line_dictionary_header()
	OC3.write_line_dictionary("CHAIN", 0.09, 77.7066, 384243000, 1.0)
	OC3.write_node_properties_header()
	OC3.write_node_properties(1, "FIX", 853.87, 0, 320.0, 0, 0, "#", "#", "#")
	OC3.write_node_properties(2, "VESSEL", 5.2, 0, -70.0, 0, 0, "#", "#", "#")
	OC3.write_line_properties_header()
	OC3.write_line_properties(1, "CHAIN", 902.2, 1, 2, "gy_pos  gz_a_pos")
	OC3.write_solver_options()

