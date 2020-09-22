
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'
__author__ = 'Michael_Lidia_Jiaxing_Elif'

# Import carrier concentration plot object
from vtandem.visualization.carrier_concentration_plot import CarrierConcentration_Plot


class Ternary_Carrier_Concentration(CarrierConcentration_Plot):
	
	def __init__(self, parent = None, main_compound = None, first_element = None, second_element = None, third_element = None):
		
		# Inherit all variables (plot object, etc.) from parent object (CarrierConcentration_Plot)
		super().__init__()
		
		# Establish the first, second, third, and fourth species of the quaternary compound.
		# Note that this list is subject to change, depending on what the user chooses.
		self.main_compound  = main_compound
		self.first_element	= first_element
		self.second_element	= second_element
		self.third_element	= third_element
		self.elements_list  = [self.first_element, self.second_element, self.third_element]
		
		
		# Keep track of chemical potential values
		self.mu_elements = {self.first_element: {"mu0": 0.0, "deltamu": 0.0},
							self.second_element: {"mu0": 0.0, "deltamu": 0.0},
							self.third_element: {"mu0": 0.0, "deltamu": 0.0} }
		
		# Number of each specie in the main ternary compound
		self.main_compound_number_first_specie  = 0
		self.main_compound_number_second_specie = 0
		self.main_compound_number_third_specie  = 0
		self.number_species = {	self.first_element: self.main_compound_number_first_specie, 
								self.second_element: self.main_compound_number_second_specie, 
								self.third_element: self.main_compound_number_third_specie }
		
		# Store all extracted DFT data
		self.main_compound_total_energy = 0.0
		self.first_element_mu0 = 0.0
		self.second_element_mu0 = 0.0
		self.third_element_mu0 = 0.0







