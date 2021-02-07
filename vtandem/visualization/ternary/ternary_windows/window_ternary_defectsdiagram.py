
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'
__author__ = 'Michael_Lidia_Jiaxing_Elif'

import numpy as np

from vtandem.visualization.windows.window_defectsdiagram import Window_DefectsDiagram




class Window_Ternary_DefectsDiagram(Window_DefectsDiagram):
	
	def __init__(self,	first_element, \
						second_element, \
						third_element, \
						main_compound, \
						DefectsDiagram_Plot, \
						show_carrier_concentration	):
		
		super().__init__(main_compound, DefectsDiagram_Plot, show_carrier_concentration)
		
		
		###############################################################################################
		########################### Initialize materials-related variables ############################
		###############################################################################################
		
		# Label the first, second, and third species of the atoms in the ternary compound
		self.first_element = first_element
		self.second_element = second_element
		self.third_element = third_element
		self.elements_list = [self.first_element, self.second_element, self.third_element]
		
		"""
		# Information about main ternary compound
		self.main_compound_number_first_specie = self.compounds_info[self.main_compound][self.first_element]	# Number of first species in ternary compound
		self.main_compound_number_second_specie = self.compounds_info[self.main_compound][self.second_element]	# Number of second species in ternary compound
		self.main_compound_number_third_specie = self.compounds_info[self.main_compound][self.third_element]	# Number of third species in ternary compound
		self.main_compound_enthalpy = self.compounds_info[self.main_compound]["enthalpy"]						# Enthalpy of ternary compound
		"""
		
		# Keep track of mu values of the species in the ternary compound
		self.deltamu_values = {}
		self.deltamu_values[first_element] = 0.0
		self.deltamu_values[second_element] = 0.0
		self.deltamu_values[third_element] = 0.0
		
		"""
		# Extrinsic dopant
		self.dopant = "None"
		self.dopant_mu0 = 0.0
		self.dopant_deltamu = 0.0
		self.extrinsic_defects = []
		"""
	
	
	
	
	
	###############################################################################################
	################################# Generate Defects Diagram ####################################
	###############################################################################################
	
	def Generate_DefectsDiagram_Plot_Function(self):
		
		# This function specifies what happens when the user clicks the "Generate Defects Diagram" button.
		
		"""
		# Check whether the chosen elements are unique initially
		if len(self.elements_list) > len(set([self.mu1_species_selection_box.currentText(), self.mu2_species_selection_box.currentText(), self.mu3_species_selection_box.currentText()])):
			QMessageBox.about(self, "WARNING", "Pick UNIQUE elements!")
			return
		"""
		
		"""
		# Update elements and chemical potentials
		self.DefectsDiagram_Plot.first_element	= self.first_element
		self.DefectsDiagram_Plot.second_element	= self.second_element
		self.DefectsDiagram_Plot.third_element	= self.third_element
		self.DefectsDiagram_Plot.mu_elements[self.first_element]["mu0"] = self.compounds_info[self.first_element]["mu0"]
		self.DefectsDiagram_Plot.mu_elements[self.second_element]["mu0"] = self.compounds_info[self.second_element]["mu0"]
		self.DefectsDiagram_Plot.mu_elements[self.third_element]["mu0"] = self.compounds_info[self.third_element]["mu0"]
		self.DefectsDiagram_Plot.mu_elements[self.first_element]["deltamu"] = self.deltamu_values[self.first_element]
		self.DefectsDiagram_Plot.mu_elements[self.second_element]["deltamu"] = self.deltamu_values[self.second_element]
		self.DefectsDiagram_Plot.mu_elements[self.third_element]["deltamu"] = self.deltamu_values[self.third_element]
		"""
		
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








