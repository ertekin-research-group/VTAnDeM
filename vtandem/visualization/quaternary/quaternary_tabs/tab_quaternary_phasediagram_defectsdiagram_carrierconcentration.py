
__author__ = 'Michael_Lidia_Jiaxing_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'

import numpy as np

import PyQt5
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

from vtandem.visualization.quaternary.quaternary_plots.plot_quaternary_phase_diagram import ChemicalPotential_Quaternary_PhaseDiagramProjected2D
from vtandem.visualization.quaternary.quaternary_plots.plot_quaternary_defects_diagram import Plot_Quaternary_DefectsDiagram
from vtandem.visualization.quaternary.quaternary_plots.plot_quaternary_carrier_concentration import Plot_Quaternary_Carrier_Concentration

from vtandem.visualization.tabs.tab_phasediagram_defectsdiagram_carrierconcentration import Tab_PhaseDiagram_DefectsDiagram_CarrierConcentration



class Tab_PhaseDiagram_DefectsDiagram_CarrierConcentration(Tab_PhaseDiagram_DefectsDiagram_CarrierConcentration):
	
	def __init__(self, parent=None, main_compound=None, first_element=None, second_element=None, third_element=None, fourth_element=None, compounds_info=None, defects_data=None, main_compound_info=None, dos_data=None, show_defects_diagram=True, show_carrier_concentration=True):
		
		###############################################################################################
		########################### Initialize materials-related variables ############################
		###############################################################################################
		
		# Initialize the main quaternary compound
		self.main_compound = main_compound
		
		# Initialize the first, second, third, and fourth species of the atoms in the quaternary compound
		self.first_element = first_element
		self.second_element = second_element
		self.third_element = third_element
		self.fourth_element = fourth_element
		self.elements_list = [self.first_element, self.second_element, self.third_element, self.fourth_element]					# Species list (order MAY change)
		
		self.PhaseDiagram = ChemicalPotential_Quaternary_PhaseDiagramProjected2D(self, main_compound = main_compound, first_element = first_element, second_element = second_element, third_element = third_element, fourth_element = fourth_element)
		
		if show_defects_diagram:
			self.DefectsDiagram = Plot_Quaternary_DefectsDiagram(main_compound = main_compound, first_element = first_element, second_element = second_element, third_element = third_element, fourth_element = fourth_element)
		
		if show_carrier_concentration:
			self.CarrierConcentration = Plot_Quaternary_Carrier_Concentration(main_compound = main_compound, first_element = first_element, second_element = second_element, third_element = third_element, fourth_element = fourth_element)
		
		
		
		super().__init__(	compounds_info = compounds_info, \
							defects_data = defects_data, \
							main_compound_info = main_compound_info, \
							dos_data = dos_data, \
							show_defects_diagram = show_defects_diagram, \
							show_carrier_concentration = show_carrier_concentration, \
							type = "quaternary"	)
		
		
		# Mu4 values (chemical potential of fourth species)
		#self.mu4_value_array = np.linspace(self.main_compound_enthalpy / self.main_compound_elements_count[self.fourth_element], -0.0, 1001)	# Array of mu4 values (for the mu4 slider bar)
		self.mu4_value_array = np.linspace(self.main_compound_enthalpy / self.main_compound_info["dft_"+self.fourth_element], -0.0, 1001)	# Array of mu4 values (for the mu4 slider bar)
		
		
		
		self.Update_PhaseDiagram_Plot_Elements_Function()
	
	
	
	
	
	###############################################################################################
	############################## Mu Values Display Settings #####################################
	###############################################################################################
	
	def Activate_MuValue_FourthElement_Settings(self):
		
		# Create the slider widget for the chemical potential (mu) of the fourth species.
		# *Note: The widget that will contain the slider is already defined in __init__ as muvalue_settings_widget.	This function is 
		#	simply to design and functionalize "sub-main" slider widget.
		
		self.fourth_element_slider_widget = QWidget()										# Create slider widget
		self.fourth_element_slider_layout = QHBoxLayout(self.fourth_element_slider_widget)	# Create (horizontal) layout for the slider widget
		
		self.mu4_species_selection_box = QComboBox()									# Create a drop-down menu so the user can select species
		self.mu4_species_selection_box.addItem(self.first_element)						# Add first species
		self.mu4_species_selection_box.addItem(self.second_element)
		self.mu4_species_selection_box.addItem(self.third_element)
		self.mu4_species_selection_box.addItem(self.fourth_element)
		self.mu4_species_selection_box.setCurrentIndex(3)								# Set initial state of the drop-down menu to the fourth species
		self.fourth_element_slider_layout.addWidget(self.mu4_species_selection_box)		# Add the drop-down menu widget to the slider widget
		
		self.fourth_element_slider_label = QLabel(u"\u0394"+"\u03BC"+"<sub>d</sub>")	# Create a label widget to display the "\Delta\mu_4" text next to the slider (unicode format)
		self.fourth_element_slider_label.setFont(self.mu_display_label_font)			# Set the font for the label
		self.fourth_element_slider_label.setAlignment(Qt.AlignCenter)					# Align the text to center
		self.fourth_element_slider_layout.addWidget(self.fourth_element_slider_label)	# Place the label widget in the main slider widget layout
		
		self.fourth_element_slider = QSlider(Qt.Horizontal)						# QSlider is the actual name of the widget for the slider object
		self.fourth_element_slider.setMinimum(0)								# Set minimum value
		self.fourth_element_slider.setMaximum(1000)								# Set maximum value
		self.fourth_element_slider.setValue(1000)								# Set the initial value
		self.fourth_element_slider.setSingleStep(1)								# Set the step value (how much one "slide" is worth)
		self.fourth_element_slider.setTickInterval(1)							# 
																				# *Note: The slider widget cannot handle floating values, i.e. you can only define
																				#	integers for the maximum, minimum, etc. The way we get around that is by
																				#	saying that the mu4 value is chosen from an array (self.mu4_value_array), and 
																				#	the integer value of the slider is the INDEX of the array.
		self.fourth_element_slider_layout.addWidget(self.fourth_element_slider)	# Add the slider widget to the main slider widget
		
		self.fourth_element_slider.valueChanged.connect(self.Update_Fourth_Species_Slider)	# How to update the slider in the case that someone moves it
		
		self.fourth_element_slider_value_label = QLabel("{0:.4f}".format(-0.0))				# Create a label widget to display the current value of mu4 up to 4 digits
		self.fourth_element_slider_value_label.setAlignment(Qt.AlignCenter)					# Align the text to center
		self.fourth_element_slider_layout.addWidget(self.fourth_element_slider_value_label)	# Add this label widget to the main slider widget layout
		
		
		# Add each mu value display to the mu value settings widget
		self.muvalue_settings_layout.addWidget(self.fourth_element_slider_widget)
	
	
	
	
	
	
	
	
	def Update_Fourth_Species_Slider(self):
		
		# This is how we update the properties of the slider given that the user touches it.
		fourth_element_mu4_value_index = self.fourth_element_slider.value()	# Obtain the (integer) value of the slider when the user uses the slider
		self.deltamu_values[self.fourth_element] = self.mu4_value_array[fourth_element_mu4_value_index]	# Update the mu4 value by 1) using the value of the slider as the index and 
																										#	2) using the index to choose the mu4 value from self.mu4_value_array
		
		#self.deltamu_values[self.third_element] = (self.main_compound_enthalpy - self.main_compound_number_first_specie*self.deltamu_values[self.first_element] - self.main_compound_number_second_specie*self.deltamu_values[self.second_element] - self.main_compound_number_fourth_specie*self.deltamu_values[self.fourth_element]) / self.main_compound_number_third_specie
		#self.deltamu_values[self.third_element] = (self.main_compound_enthalpy - self.main_compound_elements_count[self.first_element]*self.deltamu_values[self.first_element] - self.main_compound_elements_count[self.second_element]*self.deltamu_values[self.second_element] - self.main_compound_elements_count[self.fourth_element]*self.deltamu_values[self.fourth_element]) / self.main_compound_elements_count[self.third_element]
		self.deltamu_values[self.third_element] = (	self.main_compound_enthalpy \
													- self.main_compound_info["dft_"+self.first_element]*self.deltamu_values[self.first_element] \
													- self.main_compound_info["dft_"+self.second_element]*self.deltamu_values[self.second_element] \
													- self.main_compound_info["dft_"+self.fourth_element]*self.deltamu_values[self.fourth_element] \
													) / self.main_compound_info["dft_"+self.third_element]
		self.mu3_display.setText("{0:.4f}".format(self.deltamu_values[self.third_element]))		# Every time the mu4 value changes (and the mu1 and mu2 values are held constant), the mu3 value changes
		
		mu4_rounded = round(self.deltamu_values[self.fourth_element], 4)						# Display the updated mu4 value as a float up to four decimal places
		self.fourth_element_slider_value_label.setText("{0:.4f}".format(mu4_rounded))	# Update the text display
		
		self.PhaseDiagram.deltamu[4] = self.deltamu_values[self.fourth_element]
		
		# Update the quaternary phase diagram
		if (self.PhaseDiagram.main_compound_plot != None) and (self.PhaseDiagram.competing_compound_plots != {}):
			self.PhaseDiagram.Plot_PhaseDiagram()
		
		
		if self.show_defects_diagram:
			
			# Update chemical potentials in defects diagram object
			self.DefectsDiagram.mu_elements[self.third_element]["deltamu"] = self.deltamu_values[self.third_element]
			self.DefectsDiagram.mu_elements[self.fourth_element]["deltamu"] = self.deltamu_values[self.fourth_element]
			
			
			# Recalculate defect formation energies
			self.DefectsDiagram.Calculate_DefectFormations()
			
			# Redraw defects diagram
			if self.DefectsDiagram.intrinsic_defect_plots != {}:
				self.DefectsDiagram.Update_Intrinsic_DefectsDiagram_Plot()
			if self.DefectsDiagram.extrinsic_defect_plots != {}:
				self.DefectsDiagram.Update_Extrinsic_DefectsDiagram_Plot()
		
		
		if self.show_carrier_concentration:
			
			# Update chemical potentials in carrier concentration object
			self.CarrierConcentration.mu_elements[self.third_element]["deltamu"] = self.deltamu_values[self.third_element]
			self.CarrierConcentration.mu_elements[self.fourth_element]["deltamu"] = self.deltamu_values[self.fourth_element]
			
			# Redraw carrier concentration plot
			"""
			if self.CarrierConcentration.carrier_concentration_intrinsic_defect_hole_plot != None:
				self.CarrierConcentration.Update_HoleConcentration_Plot()
			if self.CarrierConcentration.carrier_concentration_intrinsic_defect_electron_plot != None:
				self.CarrierConcentration.Update_ElectronConcentration_Plot()
			"""
			self.CarrierConcentration.Update_CarrierConcentration_Plot()
		
		
		if self.show_defects_diagram and self.show_carrier_concentration:
			
			# Update the equilibrium Fermi energy
			if self.DefectsDiagram.intrinsic_defect_plots != {}:
				self.Update_Equilibrium_Fermi_Energy_Temperature()
	
	
	
	
	
	###############################################################################################
	######################################## Phase Diagram ########################################
	###############################################################################################
	
	def Update_PhaseDiagram_Plot_Elements_Function(self):
		
		# If unique, update each species and the list of elements
		self.first_element = str(self.mu1_species_selection_box.currentText())
		self.second_element = str(self.mu2_species_selection_box.currentText())
		self.third_element = str(self.mu3_species_selection_box.currentText())
		self.fourth_element = str(self.mu4_species_selection_box.currentText())
		self.elements_list = [self.first_element, self.second_element, self.third_element, self.fourth_element]
		
		# Reset the slide bar
		self.fourth_element_slider.setEnabled(True)
		self.fourth_element_slider_label.setText(u"\u0394\u03BC<sub>"+self.fourth_element+"</sub>")
		endpoint_slidebar = self.main_compound_enthalpy / self.main_compound_info["dft_"+self.fourth_element]
		self.mu4_value_array = np.linspace(endpoint_slidebar, -0.0, 1001)
		self.fourth_element_slider.setValue(1000)
		
		# Reset the mu values and their displays
		self.deltamu_values[self.first_element]  = 0.0
		self.deltamu_values[self.second_element] = 0.0
		self.deltamu_values[self.third_element]  = (	self.main_compound_enthalpy \
														- self.main_compound_info["dft_"+self.first_element]*self.deltamu_values[self.first_element] \
														- self.main_compound_info["dft_"+self.second_element]*self.deltamu_values[self.second_element] \
														- self.main_compound_info["dft_"+self.fourth_element]*self.deltamu_values[self.fourth_element] \
														) / self.main_compound_info["dft_"+self.third_element]	# Every time the mu4 value changes (and the mu1 and mu2 values are held constant), the mu3 value changes
		self.deltamu_values[self.fourth_element] = self.mu4_value_array[self.fourth_element_slider.value()]
		
		self.mu1_display_label.setText(u"\u0394\u03BC<sub>"+self.first_element+"</sub> = ")
		self.mu2_display_label.setText(u"\u0394\u03BC<sub>"+self.second_element+"</sub> = ")
		self.mu3_display_label.setText(u"\u0394\u03BC<sub>"+self.third_element+"</sub> = ")
		self.mu1_display.setText("{0:.4f}".format(self.deltamu_values[self.first_element]))
		self.mu2_display.setText("{0:.4f}".format(self.deltamu_values[self.second_element]))
		self.mu3_display.setText("{0:.4f}".format(self.deltamu_values[self.third_element]))
		
		# Define species
		self.PhaseDiagram.first_element = self.first_element
		self.PhaseDiagram.second_element = self.second_element
		self.PhaseDiagram.third_element = self.third_element
		self.PhaseDiagram.fourth_element = self.fourth_element
		self.PhaseDiagram.elements_list = self.elements_list
		
		"""
		# Set enthalpy of main compound
		self.PhaseDiagram.main_compound_enthalpy = self.main_compound_enthalpy
		
		# Set endpoints of phase diagram
		#self.PhaseDiagram.phasediagram_endpoints = min(self.main_compound_enthalpy/self.main_compound_elements_count[self.first_element], self.main_compound_enthalpy/self.main_compound_elements_count[self.second_element], self.main_compound_enthalpy/self.main_compound_elements_count[self.third_element], self.main_compound_enthalpy/self.main_compound_elements_count[self.fourth_element])
		self.PhaseDiagram.phasediagram_endpoints = min(	self.main_compound_enthalpy/self.main_compound_info["dft_"+self.first_element], \
														self.main_compound_enthalpy/self.main_compound_info["dft_"+self.second_element], \
														self.main_compound_enthalpy/self.main_compound_info["dft_"+self.third_element], \
														self.main_compound_enthalpy/self.main_compound_info["dft_"+self.fourth_element])
		"""
	
	
	def Update_DefectsDiagram_Plot_Elements_Function(self):
		
		# Update elements and chemical potentials
		
		self.DefectsDiagram.first_element	= self.first_element
		self.DefectsDiagram.second_element	= self.second_element
		self.DefectsDiagram.third_element	= self.third_element
		self.DefectsDiagram.fourth_element	= self.fourth_element
		
		for element in self.elements_list:
			self.DefectsDiagram.mu_elements[element]["mu0"] = self.compounds_info[element]["mu0"]
			self.DefectsDiagram.mu_elements[element]["deltamu"] = self.deltamu_values[element]
	
	
	
	def Update_CarrierConcentration_Plot_Elements_Function(self):
		
		if len(self.elements_list) > len(set(self.elements_list)):	# Every time someone clicks the species' buttons, the self.elements list gets updated.
																	#	The "set" function checks the self.elements list and omits any that are repeated.
																	#	If any are repeated, then the chosen species are not unique.
			QMessageBox.about(self, "WARNING", "Pick UNIQUE elements!")
			return
		
		self.CarrierConcentration.first_element = self.first_element
		self.CarrierConcentration.second_element = self.second_element
		self.CarrierConcentration.third_element = self.third_element
		self.CarrierConcentration.fourth_element = self.fourth_element
		
		for element in self.elements_list:
			self.CarrierConcentration.mu_elements[element]["mu0"] = self.compounds_info[element]["mu0"]
			self.CarrierConcentration.mu_elements[element]["deltamu"] = self.deltamu_values[element]




