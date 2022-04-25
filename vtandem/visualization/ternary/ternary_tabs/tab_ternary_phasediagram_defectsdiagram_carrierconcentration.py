
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'
__author__ = 'Michael_Lidia_Jiaxing_Elif'

import numpy as np

import PyQt5
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

from vtandem.visualization.ternary.ternary_plots.plot_ternary_phase_diagram import ChemicalPotential_Ternary_PhaseDiagramProjected2D
from vtandem.visualization.ternary.ternary_plots.plot_ternary_defects_diagram import Plot_Ternary_DefectsDiagram
from vtandem.visualization.ternary.ternary_plots.plot_ternary_carrier_concentration import Plot_Ternary_Carrier_Concentration

from vtandem.visualization.tabs.tab_phasediagram_defectsdiagram_carrierconcentration import Tab_PhaseDiagram_DefectsDiagram_CarrierConcentration



class Tab_Ternary_PhaseDiagram_DefectsDiagram_CarrierConcentration(Tab_PhaseDiagram_DefectsDiagram_CarrierConcentration):
	
	def __init__(self, main_compound=None, first_element=None, second_element=None, third_element=None, compounds_info=None, defects_data=None, main_compound_info=None, dos_data=None, show_defects_diagram = True, show_carrier_concentration = True):
		
		###############################################################################################
		########################### Initialize materials-related variables ############################
		###############################################################################################
		
		# Initialize the main ternary compound
		self.main_compound = main_compound
		
		# Label the first, second, and third species of the atoms in the ternary compound
		self.first_element = first_element
		self.second_element = second_element
		self.third_element = third_element
		self.elements_list = [self.first_element, self.second_element, self.third_element]
		
		
		#self.PhaseDiagram = ChemicalPotential_Ternary_PhaseDiagramProjected2D(self, main_compound = main_compound, first_element = first_element, second_element = second_element, third_element = third_element)
		self.PhaseDiagram = ChemicalPotential_Ternary_PhaseDiagramProjected2D(main_compound = main_compound, first_element = first_element, second_element = second_element, third_element = third_element)
		
		
		if show_defects_diagram:
			self.DefectsDiagram = Plot_Ternary_DefectsDiagram(main_compound = main_compound, first_element = first_element, second_element = second_element, third_element = third_element)
		
		
		if show_carrier_concentration:
			self.CarrierConcentration = Plot_Ternary_Carrier_Concentration(main_compound = main_compound, first_element = first_element, second_element = second_element, third_element = third_element)
		
		
		super().__init__(	compounds_info = compounds_info, \
							defects_data = defects_data, \
							main_compound_info = main_compound_info, \
							dos_data = dos_data, \
							show_defects_diagram = show_defects_diagram, \
							show_carrier_concentration = show_carrier_concentration, \
							type = "ternary"	)
		
		
		self.Update_PhaseDiagram_Plot_Elements_Function()
	
	
	
	
	
	###############################################################################################
	######################################## Phase Diagram ########################################
	###############################################################################################
	
	def Update_PhaseDiagram_Plot_Elements_Function(self):
		
		# If unique, update each species and the list of elements
		self.first_element = self.mu1_species_selection_box.currentText()
		self.second_element = self.mu2_species_selection_box.currentText()
		self.third_element = self.mu3_species_selection_box.currentText()
		self.elements_list = [self.first_element, self.second_element, self.third_element]
		
		# Reset the mu value displays
		self.deltamu_values[self.first_element]  = 0.0
		self.deltamu_values[self.second_element] = 0.0
		self.deltamu_values[self.third_element]  = (	self.main_compound_enthalpy \
														- self.main_compound_info["dft_"+self.first_element]*self.deltamu_values[self.first_element] \
														- self.main_compound_info["dft_"+self.second_element]*self.deltamu_values[self.second_element] \
													) / self.main_compound_info["dft_"+self.third_element]
		
		self.mu1_display_label.setText(u"\u0394\u03BC<sub>"+self.first_element+"</sub>")
		self.mu2_display_label.setText(u"\u0394\u03BC<sub>"+self.second_element+"</sub>")
		self.mu3_display_label.setText(u"\u0394\u03BC<sub>"+self.third_element+"</sub>")
		self.mu1_display.setText(str(self.deltamu_values[self.first_element]))
		self.mu2_display.setText(str(self.deltamu_values[self.second_element]))
		self.mu3_display.setText(str(self.deltamu_values[self.third_element]))
		
		# Define species
		self.PhaseDiagram.first_element = self.first_element
		self.PhaseDiagram.second_element = self.second_element
		self.PhaseDiagram.third_element = self.third_element
		self.PhaseDiagram.elements_list = self.elements_list
	
	
	
	def Update_DefectsDiagram_Plot_Elements_Function(self):
		
		# Update elements and chemical potentials
		self.DefectsDiagram.first_element	= self.first_element
		self.DefectsDiagram.second_element	= self.second_element
		self.DefectsDiagram.third_element	= self.third_element
		
		self.DefectsDiagram.mu_elements[self.first_element]["mu0"] = self.compounds_info[self.first_element]["mu0"]
		self.DefectsDiagram.mu_elements[self.second_element]["mu0"] = self.compounds_info[self.second_element]["mu0"]
		self.DefectsDiagram.mu_elements[self.third_element]["mu0"] = self.compounds_info[self.third_element]["mu0"]
		
		self.DefectsDiagram.Update_Deltamus(self.deltamu_values)
		"""
		self.DefectsDiagram.mu_elements[self.first_element]["deltamu"] = self.deltamu_values[self.first_element]
		self.DefectsDiagram.mu_elements[self.second_element]["deltamu"] = self.deltamu_values[self.second_element]
		self.DefectsDiagram.mu_elements[self.third_element]["deltamu"] = self.deltamu_values[self.third_element]
		"""
	
	
	
	def Update_CarrierConcentration_Plot_Elements_Function(self):
		
		if len(self.elements_list) > len(set(self.elements_list)):	# Every time someone clicks the species' buttons, the self.elements list gets updated.
																	#	The "set" function checks the self.elements list and omits any that are repeated.
																	#	If any are repeated, then the chosen species are not unique.
			QMessageBox.about(self, "WARNING", "Pick UNIQUE elements!")
			return
		
		self.CarrierConcentration.first_element = self.first_element
		self.CarrierConcentration.second_element = self.second_element
		self.CarrierConcentration.third_element = self.third_element
		
		self.CarrierConcentration.mu_elements[self.first_element]["mu0"] = self.compounds_info[self.first_element]["mu0"]
		self.CarrierConcentration.mu_elements[self.second_element]["mu0"] = self.compounds_info[self.second_element]["mu0"]
		self.CarrierConcentration.mu_elements[self.third_element]["mu0"] = self.compounds_info[self.third_element]["mu0"]
		
		self.CarrierConcentration.mu_elements[self.first_element]["deltamu"] = self.deltamu_values[self.first_element]
		self.CarrierConcentration.mu_elements[self.second_element]["deltamu"] = self.deltamu_values[self.second_element]
		self.CarrierConcentration.mu_elements[self.third_element]["deltamu"] = self.deltamu_values[self.third_element]





