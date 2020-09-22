
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'
__author__ = 'Michael_Lidia_Jiaxing_Elif'

# Import defect formation energy diagram object
from vtandem.visualization.defects_diagram import DefectFormationEnergy_Diagram

class Quaternary_Defects_Diagram(DefectFormationEnergy_Diagram):
	
	def __init__(self, parent = None, main_compound = None, first_element = None, second_element = None, third_element = None, fourth_element = None):
		
		# Inherit all variables (plot object, etc.) from parent object (DefectFormationEnergy_Diagram)
		super().__init__()
		
		# Establish the first, second, third, and fourth species of the quaternary compound.
		# Note that this list is subject to change, depending on what the user chooses.
		self.main_compound  = main_compound
		self.first_element	= first_element
		self.second_element	= second_element
		self.third_element	= third_element
		self.fourth_element	= fourth_element
		self.elements_list   = [self.first_element, self.second_element, self.third_element, self.fourth_element]
		
		# Keep track of chemical potential values
		self.mu_elements = {self.first_element: {"mu0": 0.0, "deltamu": 0.0},
							self.second_element: {"mu0": 0.0, "deltamu": 0.0},
							self.third_element: {"mu0": 0.0, "deltamu": 0.0},
							self.fourth_element: {"mu0": 0.0, "deltamu": 0.0} }
		
		# Store all extracted DFT data
		self.main_compound_total_energy = 0.0
		self.first_element_mu0 = 0.0
		self.second_element_mu0 = 0.0
		self.third_element_mu0 = 0.0
		self.fourth_element_mu0 = 0.0








