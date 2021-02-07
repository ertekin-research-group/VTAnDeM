
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'
__author__ = 'Michael_Lidia_Jiaxing_Elif'

import numpy as np

from vtandem.visualization.windows.window_defectsdiagram import Window_DefectsDiagram




class Window_Binary_DefectsDiagram(Window_DefectsDiagram):
	
	def __init__(self,	first_element, \
						second_element, \
						main_compound, \
						DefectsDiagram_Plot, \
						show_carrier_concentration	):
		
		super().__init__(main_compound, DefectsDiagram_Plot, show_carrier_concentration)
		
		
		###############################################################################################
		########################### Initialize materials-related variables ############################
		###############################################################################################
		
		# Label the first and second species of the atoms in the binary compound
		self.first_element = first_element
		self.second_element = second_element
		self.elements_list = [self.first_element, self.second_element]
		
		# Keep track of mu values of the species in the binary compound
		self.deltamu_values = {}
		self.deltamu_values[first_element] = 0.0
		self.deltamu_values[second_element] = 0.0
	
	
	
	
	###############################################################################################
	################################# Generate Defects Diagram ####################################
	###############################################################################################
	
	def Generate_DefectsDiagram_Plot_Function(self):
		
		# This function specifies what happens when the user clicks the "Generate Defects Diagram" button.
		
		
		# Reset defects diagram
		self.DefectsDiagram_Plot.defects_diagram_plot_drawing.remove()
		self.DefectsDiagram_Plot.defects_diagram_plot_drawing = self.DefectsDiagram_Plot.defects_diagram_plot_figure.add_subplot(111)
		self.DefectsDiagram_Plot.Activate_DefectsDiagram_Plot_Axes()
		
		# Calculate defect formation energies
		self.DefectsDiagram_Plot.Calculate_DefectFormations()
		
		# Plot defect formation energies
		self.DefectsDiagram_Plot.intrinsic_defect_plots = {}
		self.DefectsDiagram_Plot.extrinsic_defect_plots = {}
		self.DefectsDiagram_Plot.Initialize_Intrinsic_DefectsDiagram_Plot()
		if self.DefectsDiagram_Plot.dopant != "None":
			self.DefectsDiagram_Plot.Initialize_Extrinsic_DefectsDiagram_Plot()








