
__author__ = 'Michael_Lidia_Jiaxing_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'


import numpy as np
import matplotlib.pyplot as plt

from vtandem.visualization.ternary.ternary_plots.plot_composition_ternary_phase_diagram import Plot_Composition_Ternary_PhaseDiagram
from vtandem.visualization.ternary.ternary_plots.plot_ternary_defects_diagram import Plot_Ternary_DefectsDiagram

from vtandem.visualization.tabs.tab_phasediagram_composition import Tab_Compositional_PhaseDiagram



class Tab_Ternary_Compositional_PhaseDiagram(Tab_Compositional_PhaseDiagram):
	
	def __init__(self, main_compound=None, first_element = None, second_element = None, third_element = None, compounds_info = None, defects_data = None, main_compound_info = None, dos_data = None, show_defects_diagram = True, show_carrier_concentration = True):	# User specifies the main compound and its constituents
		
		###############################################################################################
		########################### Initialize materials-related variables ############################
		###############################################################################################
		
		# Initialize the main ternary compound
		self.main_compound = main_compound
		
		# Label the first, second, and third species of the atoms in the ternary compound
		self.first_element = first_element
		self.second_element = second_element
		self.third_element = third_element
		self.elements_list = [self.first_element, self.second_element, self.third_element]					# Species list (order MAY change)
		
		
		###############################################################################################
		################################# Compositional phase diagram #################################
		###############################################################################################
		
		self.Compositional_PhaseDiagram = Plot_Composition_Ternary_PhaseDiagram(main_compound = self.main_compound, first_element = self.first_element, second_element = self.second_element, third_element = self.third_element, compounds_info = compounds_info, main_compound_info = main_compound_info)
		
		# Defects diagram
		if show_defects_diagram:
			self.DefectsDiagram = Plot_Ternary_DefectsDiagram(main_compound = self.main_compound, first_element = self.first_element, second_element = self.second_element, third_element = self.third_element)
		
		
		Tab_Compositional_PhaseDiagram.__init__(	self, \
													type = "ternary", \
													compounds_info = compounds_info, \
													defects_data = defects_data, \
													dos_data = dos_data, \
													main_compound_info = main_compound_info, \
													show_defects_diagram = show_defects_diagram,
													show_carrier_concentration = show_carrier_concentration	)




