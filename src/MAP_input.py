# This Python file uses the following encoding: utf-8
import os, sys
from main import run


class MainFile(object):
	"""The MainFile class takes everything from FloatingSE and puts it in the
	correct format for MAP++. MAP++ then outputs a linearlized stiffness matrix
	that FloatingSE uses in its optimization analysis."""
	def __init__(self, water_depth, gravity, water_density):
		"""Initalizes x, y, and z stiffness variables to 0. Initalizes line
		type, fixed nodes, connect node, and vessel nodes as empyt lists. Sets
		water depth, gravity, and water density as WATER_DEPTH, GRAVITY, and
		WATER_DENSITY, respectively."""
		super(MainFile, self).__init__()
		self.x_stiffness = 0
		self.y_stiffness = 0
		self.z_stiffness = 0
		self.line_types = []
		self.fix_nodes = []
		self.connect_nodes = []
		self.vessel_nodes = []
		self.water_depth = water_depth
		self.gravity = gravity
		self. water_density = water_density

	def write_line_dictionary_header(self):
		"""Writes the first three lines of the input.map file:
---------------------- LINE DICTIONARY ---------------------------------------
LineType  Diam      MassDenInAir   EA            CB   CIntDamp  Ca   Cdn    Cdt
(-)       (m)       (kg/m)        (N)           (-)   (Pa-s)    (-)  (-)    (-)
		"""
		file = open("./input.map", "wb")
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
		element_axial_stiffness, cable_sea_friction_coefficient=1):
		"""Writes the forth line of the input.map file. This is where LINE_TYPE,
		DIAMETER, AIR_MASS_DENSITY, ELEMENT_AXIAL_STIFFNESS, and CABLE_SEA_
		FRICTION_COEFFICIENT is inputted. CABLE_SEA_FRICTION_COEFFICIENT defaults
		to 1 when none is given."""
		file = open("./input.map", "ab")
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
		"""Writes the node properties header:
---------------------- NODE PROPERTIES ---------------------------------------
Node  Type       X       Y       Z      M     B     FX      FY      FZ
(-)   (-)       (m)     (m)     (m)    (kg)  (mˆ3)  (N)     (N)     (N)
		"""
		file = open("./input.map", "ab")
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
		file = open("./input.map", "ab")
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
		file = open("./input.map", "ab")
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
		file = open("./input.map", "ab")
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

		#check the line number
		file.write("%d   " % anchor_node_number)
		file.write("%d   " % fairlead_node_number)
		file.write("%s\n" % control_output_text_stream)
		file.close()

	def write_solver_options(self, number_of_mooring_lines):
		"""Writes the solver options at the end of the input file, as well as 
		takes the NUMBER_OF_MOORING_LINES and places them evenly within 360
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
 ref_position 0.0 0.0 -6.0
 -172.772029
    0.000000
  -82.007818
  -63.224371
    0.000000
  -28.179102

		"""
		file = open("./input.map", "ab")
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
		n = 360/number_of_mooring_lines
		degree = n
		while degree + n <= 360:
			file.write("%d " % degree)
			degree += n
		file.write("\n")
		file.write(" krylov_accelerator 3\n")
		file.write(" ref_position 0.0 0.0 -6.0\n")
		file.write(" -172.772029\n")
		file.write("    0.000000\n")
		file.write("  -82.007818\n")
		file.write("  -63.224371\n")
		file.write("    0.000000\n")
		file.write("  -28.179102\n")
		file.close()

	def run_MAP(self):
		"""This function runs MAP++."""
		run(self.water_depth, self.gravity, self.water_density)


if __name__ == '__main__':
	"""Testing the interface using homogeneous line OC3 mooring information."""
	OC3 = MainFile(320.0, 9.806, 1025.0)
	# OC3.write_to_main_input(320.0, 9.806, 1025.0)
	OC3.write_line_dictionary_header()
	OC3.write_line_dictionary("CHAIN", 0.09, 77.7066, 384243000, 1.0)
	OC3.write_node_properties_header()
	OC3.write_node_properties(1, "FIX", 853.87, 0, 320.0, 0, 0)
	OC3.write_node_properties(2, "VESSEL", 5.2, 0, -70.0, 0, 0)
	OC3.write_line_properties_header()
	OC3.write_line_properties(1, "CHAIN", 902.2, 1, 2, "gy_pos  gz_a_pos")
	OC3.write_solver_options(3)
	OC3.run_MAP()

