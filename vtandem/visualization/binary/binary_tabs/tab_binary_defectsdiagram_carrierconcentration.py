
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'
__author__ = 'Michael_Lidia_Jiaxing_Elif'

import numpy as np

import PyQt5
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

from vtandem.visualization.binary.binary_plots.plot_binary_defects_diagram import Plot_Binary_DefectsDiagram
from vtandem.visualization.binary.binary_plots.plot_binary_carrier_concentration import Plot_Binary_Carrier_Concentration

from vtandem.visualization.windows.window_defectsdiagram_binary import Window_DefectsDiagram_Binary
from vtandem.visualization.windows.window_carrierconcentration import Window_CarrierConcentration


class Tab_Binary_DefectsDiagram_CarrierConcentration(Window_DefectsDiagram_Binary, Window_CarrierConcentration):
	
	def __init__(self, parent=None, main_compound=None, first_element=None, second_element=None, compounds_info=None, defects_data=None, main_compound_info=None, dos_data=None, show_defects_diagram=True, show_carrier_concentration=True):
		
		QWidget.__init__(self)
		
		# Font description for phase stability diagram plot
		self.font = {'family': 'Arial', 'color':  'black', 'weight': 'normal', 'size': 16 }
		
		# Establish a variable for the main binary compound
		self.main_compound = main_compound
		
		# Label the first and second species of the atoms in the binary compound
		self.first_element = first_element
		self.second_element = second_element
		self.elements_list = [self.first_element, self.second_element]
		
		# Data
		self.compounds_info = compounds_info
		self.defects_data = defects_data
		self.main_compound_info = main_compound_info
		self.dos_data = dos_data
		
		# Display settings
		self.show_defects_diagram = show_defects_diagram
		self.show_carrier_concentration = show_carrier_concentration
		
		# Get enthalpy of main compound
		enthalpy_tracker = main_compound_info["dft_BulkEnergy"]
		for element in self.elements_list:
			enthalpy_tracker -= self.main_compound_info["dft_"+element] * self.compounds_info[element]["mu0"]
		self.main_compound_enthalpy = enthalpy_tracker

		# Keep track of mu values of the species in the binary compound
		self.deltamu_values = {}
		self.deltamu_values[first_element] = 0.0
		self.deltamu_values[second_element] = self.main_compound_enthalpy / self.main_compound_info["dft_"+self.second_element]
		
		# Create delta mu value arrays for first and second elements
		self.deltamu_arrays = {}
		self.deltamu_arrays[first_element] = np.linspace(self.main_compound_enthalpy / self.main_compound_info["dft_"+self.first_element], -0.0, 1001)
		self.deltamu_arrays[second_element] = np.linspace(self.main_compound_enthalpy / self.main_compound_info["dft_"+self.second_element], -0.0, 1001)
		
		# Extrinsic dopant
		self.dopant = "None"
		self.dopant_mu0 = 0.0
		self.dopant_deltamu = 0.0
		self.extrinsic_defects = []
		
		

		###############################################################################################
		###############################################################################################
		#################################### Initialize first tab #####################################
		###############################################################################################
		###############################################################################################
		
		self.tab1 = QWidget()
		self.tab1_layout = QHBoxLayout(self.tab1)
		
		
		if self.show_defects_diagram:
			
			# DEFECTS DIAGRAM
			self.tab1_defectsdiagram_widget = QWidget()												# A layout similar to that of the quabinary phase diagram
																									# 	will be used for the defects diagram. This "sub-main" 
																									#	widget will contain the following widgets:
																									#		- Defects diagram plot
																									#		- y-axis dialog
																									#		- "Generate defects diagram!" button
																									#		- "Save figure!" button
			self.tab1_defectsdiagram_widget_layout = QVBoxLayout(self.tab1_defectsdiagram_widget)	# Again, the above mentioned widgets will be stacked vertically.
			
			# Set up defects diagram object
			self.DefectsDiagram = Plot_Binary_DefectsDiagram(main_compound = main_compound, first_element = first_element, second_element = second_element)
			self.DefectsDiagram.defects_data = defects_data
			self.DefectsDiagram.main_compound_info = main_compound_info
			self.DefectsDiagram.mu_elements[self.first_element]["mu0"] = self.compounds_info[self.first_element]["mu0"]
			self.DefectsDiagram.mu_elements[self.second_element]["mu0"] = self.compounds_info[self.second_element]["mu0"]
			self.DefectsDiagram.mu_elements[self.first_element]["deltamu"] = self.deltamu_values[self.first_element]
			self.DefectsDiagram.mu_elements[self.second_element]["deltamu"] = self.deltamu_values[self.second_element]
			self.DefectsDiagram.EVBM = main_compound_info["VBM"]
			self.DefectsDiagram.ECBM = self.DefectsDiagram.EVBM + main_compound_info["BandGap"]
			self.DefectsDiagram.axis_lims["XMin"] = 0.0
			self.DefectsDiagram.axis_lims["XMax"] = main_compound_info["BandGap"]
			self.DefectsDiagram.fermi_energy_array = np.linspace(self.DefectsDiagram.EVBM, self.DefectsDiagram.ECBM, 2000)
			self.DefectsDiagram.Activate_DefectsDiagram_Plot_Axes()
			
			Window_DefectsDiagram_Binary.__init__(self, show_dopant = True)
			
			# Add the defects diagram widget to Tab 1
			self.tab1_layout.addWidget(self.defectsdiagram_window)
		
		
		
		if self.show_carrier_concentration:
			
			# CARRIER CONCENTRATION
			self.tab1_carrierconcentration_widget = QWidget()													# This "sub-main" widget to contain widgets that 
																												#	are related to the carrier concentration will
																												#	contain:
																												#		- Carrier concentration plot
																												#		- Equilibrium fermi energy display
																												#		- "Generate carrier concentration!" button
																												#		- "Save figure!" button
			self.tab1_carrierconcentration_widget_layout = QVBoxLayout(self.tab1_carrierconcentration_widget)	# The above mentioned widgets will be stacked vertically.
			
			# Set up carrier concentration plot object
			self.CarrierConcentration = Plot_Binary_Carrier_Concentration(main_compound = main_compound, first_element = first_element, second_element = second_element)
			self.CarrierConcentration.defects_data = defects_data
			self.CarrierConcentration.main_compound_info = main_compound_info
			self.CarrierConcentration.dos_data = dos_data[self.main_compound]

			self.CarrierConcentration.mu_elements[self.first_element]["mu0"] = self.compounds_info[self.first_element]["mu0"]
			self.CarrierConcentration.mu_elements[self.second_element]["mu0"] = self.compounds_info[self.second_element]["mu0"]
			
			self.CarrierConcentration.mu_elements[self.first_element]["deltamu"] = self.deltamu_values[self.first_element]
			self.CarrierConcentration.mu_elements[self.second_element]["deltamu"] = self.deltamu_values[self.second_element]
			
			self.CarrierConcentration.vol = main_compound_info["Volume"]
			self.CarrierConcentration.EVBM = self.DefectsDiagram.EVBM
			self.CarrierConcentration.ECBM = self.DefectsDiagram.ECBM
			# Fermi energies should sample outside band gap, in case EFeq is not in gap
			self.CarrierConcentration.fermi_energy_array = np.linspace(self.DefectsDiagram.EVBM-1.0, self.DefectsDiagram.ECBM+1.0, 2000)
			self.CarrierConcentration.Activate_CarrierConcentration_Plot_Axes()
			self.CarrierConcentration.Organize_DOS_Data()
			self.CarrierConcentration.Extract_Relevant_Energies_DOSs()
			self.CarrierConcentration.Calculate_Hole_Electron_Concentration_Matrices()

			Window_CarrierConcentration.__init__(self)

			self.tab1_layout.addWidget(self.carrierconcentration_window)
	
	
	
	
	###############################################################################################
	################################## Chemical Potential Slider ##################################
	###############################################################################################
	
	def Chemical_Potential_Slider(self, element):
		
		# Arguments:
		#	element:	Name of element (str)
		
		if element == self.first_element:
			chemical_potential_first_element_index = self.chemical_potential_first_element_slider.value()
			chemical_potential_second_element_index = 1000 - chemical_potential_first_element_index
			self.chemical_potential_second_element_slider.setValue(chemical_potential_second_element_index)
		
		elif element == self.second_element:
			chemical_potential_second_element_index = self.chemical_potential_second_element_slider.value()
			chemical_potential_first_element_index = 1000 - chemical_potential_second_element_index
			self.chemical_potential_first_element_slider.setValue(chemical_potential_first_element_index)
		
		self.deltamu_values[self.first_element] = self.deltamu_arrays[self.first_element][chemical_potential_first_element_index]
		self.deltamu_values[self.second_element] = self.deltamu_arrays[self.second_element][chemical_potential_second_element_index]
		
		self.chemical_potential_first_element_deltamu.setText('{:.4f}'.format(round(self.deltamu_values[self.first_element], 4)))
		self.chemical_potential_second_element_deltamu.setText('{:.4f}'.format(round(self.deltamu_values[self.second_element], 4)))
		
		# Update deltamu values in DefectsDiagram object
		self.DefectsDiagram.mu_elements[self.first_element]["deltamu"] = self.deltamu_values[self.first_element]
		self.DefectsDiagram.mu_elements[self.second_element]["deltamu"] = self.deltamu_values[self.second_element]
		
		# Calculate defect formation energies
		self.DefectsDiagram.Calculate_DefectFormations()
		
		# Update defects diagram
		try:
			self.DefectsDiagram.Update_Intrinsic_DefectsDiagram_Plot()
		except:
			self.DefectsDiagram.Initialize_Intrinsic_DefectsDiagram_Plot()

		# Update extrinsic defect formation energies
		if self.dopant != "None":
			try:
				self.DefectsDiagram.Update_Extrinsic_DefectsDiagram_Plot()
			except:
				self.DefectsDiagram.Initialize_Extrinsic_DefectsDiagram_Plot()
		
		
		
		if self.show_carrier_concentration:
			
			# Update deltamu values in CarrierConcentration object
			self.CarrierConcentration.mu_elements[self.first_element]["deltamu"] = self.deltamu_values[self.first_element]
			self.CarrierConcentration.mu_elements[self.second_element]["deltamu"] = self.deltamu_values[self.second_element]
			
			# Calculate carrier concentrations
			try:
				self.CarrierConcentration.Update_CarrierConcentration_Plot()
			except:
				self.CarrierConcentration.Initialize_CarrierConcentration_Plot()
			
			# Plot the equilibrium Fermi energy
			self.Update_Equilibrium_Fermi_Energy_Temperature()
	
	
	
	###############################################################################################
	################################# Generate Defects Diagram ####################################
	###############################################################################################
	
	def Generate_DefectsDiagram_Plot_Function_Binary(self, event):
		
		self.Generate_DefectsDiagram_Plot_Function()
		
		if self.show_carrier_concentration:
			
			# Calculate carrier concentrations
			self.CarrierConcentration.Initialize_CarrierConcentration_Plot()

			# Plot the equilibrium Fermi energy
			self.Update_Equilibrium_Fermi_Energy_Temperature()

			# Allow synthesis temperature to be set
			self.defects_synthesis_temperature_box.setEnabled(True)

