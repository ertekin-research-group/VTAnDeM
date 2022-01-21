
__author__ = 'Michael_Lidia_Jiaxing_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'


import numpy as np

from vtandem.visualization.quaternary.quaternary_plots.plot_composition_quaternary_phase_diagram import Plot_Composition_Quaternary_PhaseDiagram
from vtandem.visualization.quaternary.quaternary_plots.plot_quaternary_defects_diagram import Plot_Quaternary_DefectsDiagram
from vtandem.visualization.quaternary.quaternary_plots.plot_quaternary_carrier_concentration import Plot_Quaternary_Carrier_Concentration

from vtandem.visualization.tabs.tab_phasediagram_composition import Tab_Compositional_PhaseDiagram


class Tab_Quaternary_Compositional_PhaseDiagram3D(Tab_Compositional_PhaseDiagram):
	
	def __init__(self, main_compound = None, first_element = None, second_element = None, third_element = None, fourth_element = None, compounds_info = None, defects_data = None, main_compound_info = None, dos_data = None, show_defects_diagram = True, show_carrier_concentration = True):	# User specifies the main compound and its constituents
		
		###############################################################################################
		########################### Initialize materials-related variables ############################
		###############################################################################################
		
		# Initialize the main quaternary compound
		self.main_compound = main_compound
		
		# Label the first, second, third, and fourth species of the atoms in the quaternary compound
		self.first_element = first_element
		self.second_element = second_element
		self.third_element = third_element
		self.fourth_element = fourth_element
		self.elements_list = [self.first_element, self.second_element, self.third_element, self.fourth_element]					# Species list (order MAY change)
		
		print(main_compound_info)

		# Compositional phase diagram
		self.Compositional_PhaseDiagram = Plot_Composition_Quaternary_PhaseDiagram(main_compound = self.main_compound, first_element = self.first_element, second_element = self.second_element, third_element = self.third_element, fourth_element = self.fourth_element, compounds_info = compounds_info, main_compound_info = main_compound_info)
		
		# Defects diagram
		if show_defects_diagram:
			self.DefectsDiagram = Plot_Quaternary_DefectsDiagram(main_compound = self.main_compound, first_element = self.first_element, second_element = self.second_element, third_element = self.third_element, fourth_element = self.fourth_element)
		
		# Carrier concentration
		if show_carrier_concentration:
			self.CarrierConcentration = Plot_Quaternary_Carrier_Concentration(main_compound = self.main_compound, first_element = self.first_element, second_element = self.second_element, third_element = self.third_element, fourth_element = self.fourth_element)
		
		
		
		Tab_Compositional_PhaseDiagram.__init__(	self, \
													type = "quaternary", \
													compounds_info = compounds_info, \
													defects_data = defects_data, \
													main_compound_info = main_compound_info, \
													dos_data = dos_data, \
													show_defects_diagram = show_defects_diagram, \
													show_carrier_concentration = show_carrier_concentration	)




