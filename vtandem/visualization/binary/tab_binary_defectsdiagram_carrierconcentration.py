
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'
__author__ = 'Michael_Lidia_Jiaxing_Elif'

import numpy as np
from labellines import labelLine, labelLines

import PyQt5
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

from vtandem.visualization.binary.binary_scripts.plot_binary_defects_diagram import Binary_Defects_Diagram
from vtandem.visualization.binary.binary_scripts.plot_binary_carrier_concentration import Binary_Carrier_Concentration



class Tab_Binary_DefectsDiagram_CarrierConcentration(QWidget):
	
	def __init__(self, parent=None, main_compound=None, first_element=None, second_element=None, compounds_info=None, defects_data=None, dos_data=None, show_defects_diagram=True, show_carrier_concentration=True):
		
		QWidget.__init__(self)
		
		# Font description for phase stability diagram plot
		self.font = {'family': 'Arial', 'color':  'black', 'weight': 'normal', 'size': 16 }
		
		# Display variables
		self.show_defects_diagram = show_defects_diagram
		self.show_carrier_concentration = show_carrier_concentration
		
		# Establish a variable for the main binary compound
		self.main_compound = main_compound
		
		# Label the first and second species of the atoms in the binary compound
		self.first_element = first_element
		self.second_element = second_element
		self.elements_list = [self.first_element, self.second_element]
		
		# Obtain DFT data
		self.compounds_info = compounds_info
		self.defects_data = defects_data
		self.dos_data = dos_data
		
		# Information about main binary compound
		self.main_compound_number_first_specie = self.compounds_info[self.main_compound][self.first_element]	# Number of first species in binary compound
		self.main_compound_number_second_specie = self.compounds_info[self.main_compound][self.second_element]	# Number of second species in binary compound
		self.main_compound_enthalpy = self.compounds_info[self.main_compound]["enthalpy"]						# Enthalpy of binary compound
		
		# Keep track of mu values of the species in the binary compound
		self.deltamu_values = {}
		self.deltamu_values[first_element] = 0.0
		self.deltamu_values[second_element] = self.main_compound_enthalpy / self.main_compound_number_second_specie
		
		# Create delta mu value arrays for first and second elements
		self.deltamu_arrays = {}
		self.deltamu_arrays[first_element] = np.linspace(self.main_compound_enthalpy / self.main_compound_number_first_specie, -0.0, 1001)
		self.deltamu_arrays[second_element] = np.linspace(self.main_compound_enthalpy / self.main_compound_number_second_specie, -0.0, 1001)
		
		# Extrinsic dopant
		self.extrinsic_defect = "None"
		self.extrinsic_defect_mu0 = 0.0
		self.extrinsic_defect_deltamu = 0.0
		
		
		
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
			self.DefectsDiagram = Binary_Defects_Diagram(self, main_compound = main_compound, first_element = first_element, second_element = second_element)
			self.DefectsDiagram.defects_data = self.defects_data[self.main_compound]
			self.DefectsDiagram.main_compound_number_first_specie = self.main_compound_number_first_specie
			self.DefectsDiagram.main_compound_number_second_specie = self.main_compound_number_second_specie
			self.DefectsDiagram.main_compound_total_energy = self.compounds_info[main_compound]["total_energy"]
			self.DefectsDiagram.mu_elements[self.first_element]["mu0"] = self.compounds_info[self.first_element]["mu0"]
			self.DefectsDiagram.mu_elements[self.second_element]["mu0"] = self.compounds_info[self.second_element]["mu0"]
			self.DefectsDiagram.mu_elements[self.first_element]["deltamu"] = self.deltamu_values[self.first_element]
			self.DefectsDiagram.mu_elements[self.second_element]["deltamu"] = self.deltamu_values[self.second_element]
			self.DefectsDiagram.EVBM = self.compounds_info[main_compound]["vbm"]
			self.DefectsDiagram.ECBM = self.DefectsDiagram.EVBM + self.compounds_info[main_compound]["band_gap"]
			self.DefectsDiagram.fermi_energy_array = np.linspace(self.DefectsDiagram.EVBM, self.DefectsDiagram.ECBM, 100)
			self.DefectsDiagram.Activate_DefectsDiagram_Plot_Axes()
			
			# Defects diagram plot
			self.defects_diagram_plot = self.DefectsDiagram.defects_diagram_plot_canvas
			self.tab1_defectsdiagram_widget_layout.addWidget(self.defects_diagram_plot)
			
			# Y-axis limits for defects diagram
			self.defectsdiagram_viewport = QWidget()
			self.defectsdiagram_viewport_layout = QHBoxLayout(self.defectsdiagram_viewport)
			
			# Y-axis limits for defects diagram
			self.defectsdiagram_Ymin_label = QLabel(u"y"+"<sub>min</sub>")
			self.defectsdiagram_Ymin_label.setAlignment(Qt.AlignRight)
			self.defectsdiagram_viewport_layout.addWidget(self.defectsdiagram_Ymin_label)
			self.defectsdiagram_Ymin_box = QLineEdit("-2.0")
			self.defectsdiagram_Ymin_box.editingFinished.connect(lambda: self.Update_WindowSize("DefectsDiagram", "YMin"))
			self.defectsdiagram_viewport_layout.addWidget(self.defectsdiagram_Ymin_box)
			self.defectsdiagram_Ymax_label = QLabel(u"y"+"<sub>max</sub>")
			self.defectsdiagram_Ymax_label.setAlignment(Qt.AlignRight)
			self.defectsdiagram_viewport_layout.addWidget(self.defectsdiagram_Ymax_label)
			self.defectsdiagram_Ymax_box = QLineEdit("2.0")
			self.defectsdiagram_Ymax_box.editingFinished.connect(lambda: self.Update_WindowSize("DefectsDiagram", "YMax"))
			self.defectsdiagram_viewport_layout.addWidget(self.defectsdiagram_Ymax_box)
			self.tab1_defectsdiagram_widget_layout.addWidget(self.defectsdiagram_viewport)
			
			
			
			
			
			
			
			
			
			
			
			
			
			
			# Chemical potential tuners
			self.chemical_potential_displays = QWidget()
			self.chemical_potential_displays_layout = QVBoxLayout(self.chemical_potential_displays)
			self.chemical_potential_first_element = QWidget()
			self.chemical_potential_first_element_layout = QHBoxLayout(self.chemical_potential_first_element)
			self.chemical_potential_first_element_label = QLabel(u"\u0394"+"\u03BC"+"<sub>"+self.first_element+"</sub> =")
			self.chemical_potential_first_element_label.setAlignment(Qt.AlignCenter)
			self.chemical_potential_first_element_layout.addWidget(self.chemical_potential_first_element_label)
			self.chemical_potential_first_element_deltamu = QLabel("-0.0000")
			self.chemical_potential_first_element_layout.addWidget(self.chemical_potential_first_element_deltamu)
			self.chemical_potential_displays_layout.addWidget(self.chemical_potential_first_element)
			self.chemical_potential_second_element = QWidget()
			self.chemical_potential_second_element_layout = QHBoxLayout(self.chemical_potential_second_element)
			self.chemical_potential_second_element_label = QLabel(u"\u0394"+"\u03BC"+"<sub>"+self.second_element+"</sub> =")
			self.chemical_potential_second_element_label.setAlignment(Qt.AlignCenter)
			self.chemical_potential_second_element_layout.addWidget(self.chemical_potential_second_element_label)
			self.chemical_potential_second_element_deltamu = QLabel('{:.4f}'.format(round(self.deltamu_values[second_element], 4)))
			self.chemical_potential_second_element_layout.addWidget(self.chemical_potential_second_element_deltamu)
			self.chemical_potential_displays_layout.addWidget(self.chemical_potential_second_element)
			
			# Chemical potential tuners
			self.chemical_potential_sliders = QWidget()
			self.chemical_potential_sliders_layout = QVBoxLayout(self.chemical_potential_sliders)
			self.chemical_potential_first_element_slider = QSlider(Qt.Horizontal)
			self.chemical_potential_first_element_slider.setMinimum(0)
			self.chemical_potential_first_element_slider.setMaximum(1000)
			self.chemical_potential_first_element_slider.setValue(1000)
			self.chemical_potential_first_element_slider.setSingleStep(1)
			self.chemical_potential_first_element_slider.setTickInterval(1)
			self.chemical_potential_first_element_slider.valueChanged.connect(lambda: self.Chemical_Potential_Slider(self.first_element))
			self.chemical_potential_sliders_layout.addWidget(self.chemical_potential_first_element_slider)
			self.chemical_potential_second_element_slider = QSlider(Qt.Horizontal)
			self.chemical_potential_second_element_slider.setMinimum(0)
			self.chemical_potential_second_element_slider.setMaximum(1000)
			self.chemical_potential_second_element_slider.setValue(0)
			self.chemical_potential_second_element_slider.setSingleStep(1)
			self.chemical_potential_second_element_slider.setTickInterval(1)
			self.chemical_potential_second_element_slider.valueChanged.connect(lambda: self.Chemical_Potential_Slider(self.second_element))
			self.chemical_potential_sliders_layout.addWidget(self.chemical_potential_second_element_slider)
			
			# Chemical potentials section
			self.chemical_potentials_section = QWidget()
			self.chemical_potentials_section_layout = QHBoxLayout(self.chemical_potentials_section)
			self.chemical_potentials_section_layout.addWidget(self.chemical_potential_displays)
			self.chemical_potentials_section_layout.addWidget(self.chemical_potential_sliders)
			self.tab1_defectsdiagram_widget_layout.addWidget(self.chemical_potentials_section)
			
			
			
			
			
			
			
			
			
			
			# Extrinsic defect properties
			self.extrinsic_defect_properties = QWidget()
			self.extrinsic_defect_properties_layout = QHBoxLayout(self.extrinsic_defect_properties)
			
			# Extrinsic defect chemical potential
			self.extrinsic_defect_chemical_potential_label = QLabel(u"\u0394"+"\u03BC"+"<sub>x</sub>")
			self.extrinsic_defect_chemical_potential_label.setAlignment(Qt.AlignCenter)
			self.extrinsic_defect_properties_layout.addWidget(self.extrinsic_defect_chemical_potential_label)
			self.extrinsic_defect_chemical_potential_deltamu = QLineEdit("-0.0000")
			self.extrinsic_defect_chemical_potential_deltamu.setMaxLength(7)
			self.extrinsic_defect_chemical_potential_deltamu.editingFinished.connect(self.Update_ExtrinsicDefect_DeltaMu)
			self.extrinsic_defect_properties_layout.addWidget(self.extrinsic_defect_chemical_potential_deltamu)
			
			# Extrinsic defect selection box
			self.extrinsic_defect_selection_box = QComboBox()
			self.extrinsic_defect_selection_box.addItem("None")
			for defect in self.defects_data[self.main_compound].keys():
				if "_" not in defect:
					continue
				if self.defects_data[self.main_compound][defect]["Extrinsic"] == "Yes":
					self.extrinsic_defect_selection_box.addItem(defect)
			self.extrinsic_defect_selection_box.setCurrentIndex(0)
			self.extrinsic_defect_selection_box.activated.connect(self.Update_ExtrinsicDefect)
			self.extrinsic_defect_properties_layout.addWidget(self.extrinsic_defect_selection_box)
			
			
			if self.show_carrier_concentration:
				
				# Synthesis temperature
				self.defects_synthesis_temperature = QWidget()
				self.defects_synthesis_temperature_layout = QHBoxLayout(self.defects_synthesis_temperature)
				self.defects_synthesis_temperature_label = QLabel(u"T<sub>syn</sub> (K) = ")
				self.defects_synthesis_temperature_label.setAlignment(Qt.AlignRight)
				self.defects_synthesis_temperature_layout.addWidget(self.defects_synthesis_temperature_label)
				self.defects_synthesis_temperature_box = QLineEdit("")
				self.defects_synthesis_temperature_box.editingFinished.connect(self.Update_SynthesisTemperature)
				self.defects_synthesis_temperature_layout.addWidget(self.defects_synthesis_temperature_box)
				self.extrinsic_defect_properties_layout.addWidget(self.defects_synthesis_temperature)
			
			
			self.tab1_defectsdiagram_widget_layout.addWidget(self.extrinsic_defect_properties)
			
			
			
			
			# (WIDGET) Button to generate defects diagram
			self.generate_defects_diagram_plot_button_widget = QPushButton("Generate Defects Diagram")
			self.generate_defects_diagram_plot_button_widget.clicked[bool].connect(self.Generate_DefectsDiagram_Plot_Function)
			self.tab1_defectsdiagram_widget_layout.addWidget(self.generate_defects_diagram_plot_button_widget)
			
			# (WIDGET) Save defects diagram as figure
			self.defects_diagram_savefigure_button = QPushButton("Save Defects Diagram Figure")
			self.defects_diagram_savefigure_button.clicked[bool].connect(lambda: self.SaveFigure("Defects Diagram"))
			self.tab1_defectsdiagram_widget_layout.addWidget(self.defects_diagram_savefigure_button)
			
			# Add the defects diagram widget to Tab 1
			self.tab1_layout.addWidget(self.tab1_defectsdiagram_widget)
		
		
		
		
		
		
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
			self.CarrierConcentration = Binary_Carrier_Concentration(self, main_compound = main_compound, first_element = first_element, second_element = second_element)
			self.CarrierConcentration.binary_defects_data = self.defects_data[self.main_compound]
			self.CarrierConcentration.binary_dos_data = self.dos_data[self.main_compound]
			self.CarrierConcentration.main_compound_number_first_specie = self.main_compound_number_first_specie
			self.CarrierConcentration.main_compound_number_second_specie = self.main_compound_number_second_specie
			self.CarrierConcentration.main_compound_total_energy = self.compounds_info[main_compound]["total_energy"]
			self.CarrierConcentration.mu_elements[self.first_element]["mu0"] = self.compounds_info[self.first_element]["mu0"]
			self.CarrierConcentration.mu_elements[self.second_element]["mu0"] = self.compounds_info[self.second_element]["mu0"]
			self.CarrierConcentration.mu_elements[self.first_element]["deltamu"] = self.deltamu_values[self.first_element]
			self.CarrierConcentration.mu_elements[self.second_element]["deltamu"] = self.deltamu_values[self.second_element]
			self.CarrierConcentration.vol = self.compounds_info[self.main_compound]["volume"]
			self.CarrierConcentration.EVBM = self.DefectsDiagram.EVBM
			self.CarrierConcentration.ECBM = self.DefectsDiagram.ECBM
			self.CarrierConcentration.fermi_energy_array = np.linspace(self.CarrierConcentration.EVBM, self.CarrierConcentration.ECBM, 100)
			self.CarrierConcentration.Activate_CarrierConcentration_Plot_Axes()
			self.CarrierConcentration.Organize_DOS_Data()
			self.CarrierConcentration.Extract_Relevant_Energies_DOSs()
			self.CarrierConcentration.Calculate_Hole_Electron_Concentration_Matrices()
			
			# Carrier concentration plot
			self.carrier_concentration_plot = self.CarrierConcentration.carrier_concentration_plot_canvas
			
			# (WIDGET) Carrier concentration plot
			self.tab1_carrierconcentration_widget_layout.addWidget(self.carrier_concentration_plot)
			
			self.carrierconcentration_viewport = QWidget()
			self.carrierconcentration_viewport_layout = QHBoxLayout(self.carrierconcentration_viewport)
			
			# (WIDGET) Y-axis limits for carrier concentration
			self.carrierconcentration_Ymin_label = QLabel(u"y"+"<sub>min</sub>")
			self.carrierconcentration_Ymin_label.setAlignment(Qt.AlignRight)
			self.carrierconcentration_viewport_layout.addWidget(self.carrierconcentration_Ymin_label)
			self.carrierconcentration_Ymin_box = QLineEdit("1E16")
			self.carrierconcentration_Ymin_box.editingFinished.connect(lambda: self.Update_WindowSize("CarrierConcentration", "YMin"))
			self.carrierconcentration_viewport_layout.addWidget(self.carrierconcentration_Ymin_box)
			self.carrierconcentration_Ymax_label = QLabel(u"y"+"<sub>max</sub>")
			self.carrierconcentration_Ymax_label.setAlignment(Qt.AlignRight)
			self.carrierconcentration_viewport_layout.addWidget(self.carrierconcentration_Ymax_label)
			self.carrierconcentration_Ymax_box = QLineEdit("1E23")
			self.carrierconcentration_Ymax_box.editingFinished.connect(lambda: self.Update_WindowSize("CarrierConcentration", "YMax"))
			self.carrierconcentration_viewport_layout.addWidget(self.carrierconcentration_Ymax_box)
			
			self.carrierconcentration_holes_checkbox = QCheckBox("Holes",self)
			self.carrierconcentration_holes_checkbox.setChecked(True)
			self.carrierconcentration_viewport_layout.addWidget(self.carrierconcentration_holes_checkbox)
			self.carrierconcentration_electrons_checkbox = QCheckBox("Electrons",self)
			self.carrierconcentration_electrons_checkbox.setChecked(True)
			self.carrierconcentration_viewport_layout.addWidget(self.carrierconcentration_electrons_checkbox)
			
			
			self.tab1_carrierconcentration_widget_layout.addWidget(self.carrierconcentration_viewport)
			
			
			
			self.Activate_Equilibrium_Fermi_Energy_Settings()
			
			
			
			# (WIDGET) Button to generate carrier concentration plot
			self.generate_carrier_concentration_plot_button_widget = QPushButton("Generate Carrier Concentration")
			self.generate_carrier_concentration_plot_button_widget.clicked[bool].connect(self.Generate_CarrierConcentration_Plot_Function)
			self.tab1_carrierconcentration_widget_layout.addWidget(self.generate_carrier_concentration_plot_button_widget)
			
			# (WIDGET) Save carrier concentration plot as figure
			self.carrier_concentration_savefigure_button = QPushButton("Save Carrier Concentration Plot")
			self.carrier_concentration_savefigure_button.clicked[bool].connect(lambda: self.SaveFigure("Carrier Concentration"))
			self.tab1_carrierconcentration_widget_layout.addWidget(self.carrier_concentration_savefigure_button)
			
			self.tab1_layout.addWidget(self.tab1_carrierconcentration_widget)
	
	
	
	
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
		self.DefectsDiagram.Update_Intrinsic_DefectsDiagram_Plot()
		
		# Update extrinsic defect formation energies
		if self.extrinsic_defect != "None":
			self.DefectsDiagram.Update_Extrinsic_DefectsDiagram_Plot()
		
		
		
		if self.show_carrier_concentration:
			
			# Update deltamu values in CarrierConcentration object
			self.CarrierConcentration.mu_elements[self.first_element]["deltamu"] = self.deltamu_values[self.first_element]
			self.CarrierConcentration.mu_elements[self.second_element]["deltamu"] = self.deltamu_values[self.second_element]
			
			# Calculate carrier concentrations
			if self.carrierconcentration_holes_checkbox.isChecked():
				self.CarrierConcentration.Update_HoleConcentration_Plot()
		
			if self.carrierconcentration_electrons_checkbox.isChecked():
				self.CarrierConcentration.Update_ElectronConcentration_Plot()
			
			# Plot the equilibrium Fermi energy
			self.Update_Equilibrium_Fermi_Energy_Temperature()
	
	
	
	
	
	###############################################################################################
	################################### Equilibrium Fermi Energy ##################################
	###############################################################################################
	
	def Activate_Equilibrium_Fermi_Energy_Settings(self):
		
		# (WIDGET) Equilibrium Fermi energy widget
		self.equilibrium_fermi_energy_widget = QWidget()
		self.equilibrium_fermi_energy_widget_layout = QHBoxLayout(self.equilibrium_fermi_energy_widget)
		self.equilibrium_fermi_energy_widget_layout.setContentsMargins(0, 0, 0, 0)
		self.equilibrium_fermi_energy_widget_layout.setSpacing(0)
		
		# Label for user-selectable temperature
		temperature_label = QLabel("T (K) = ")
		temperature_label.setMargin(0)
		temperature_label.setAlignment(Qt.AlignCenter)
		self.equilibrium_fermi_energy_widget_layout.addWidget(temperature_label)
		
		# Temperature selection prompt
		self.temperature_selection_box = QComboBox()
		for temperature in self.CarrierConcentration.temperature_array:
			self.temperature_selection_box.addItem(str(int(temperature)))
		self.temperature_selection_box.setCurrentIndex(2)
		self.temperature_selection_box.activated.connect(self.Update_Equilibrium_Fermi_Energy_Temperature)
		self.equilibrium_fermi_energy_widget_layout.addWidget(self.temperature_selection_box)
		
		# Spacer item
		self.equilibrium_fermi_energy_widget_layout.addItem(QSpacerItem(50, 20, QSizePolicy.Expanding, QSizePolicy.Minimum))
		
		# Label for equilibrium Fermi energy
		equilibrium_fermi_energy_label = QLabel(u"E"+"<sub>f</sub>"+"<sup>eq</sup> (eV) = ")
		equilibrium_fermi_energy_label.setAlignment(Qt.AlignCenter)
		self.equilibrium_fermi_energy_widget_layout.addWidget(equilibrium_fermi_energy_label)
		
		# Display for equilibrium Fermi energy
		self.equilibrium_fermi_energy_display = QLineEdit()
		self.equilibrium_fermi_energy_display.setMaxLength(7)
		self.equilibrium_fermi_energy_display.setEnabled(False)
		self.equilibrium_fermi_energy_display.setStyleSheet("""QLineEdit { background-color: white; color: black }""")
		self.equilibrium_fermi_energy_widget_layout.addWidget(self.equilibrium_fermi_energy_display)
		
		self.equilibrium_fermi_energy_widget_layout.addItem(QSpacerItem(50, 40, QSizePolicy.Expanding, QSizePolicy.Minimum))
		
		
		self.equilibrium_fermi_energy_checkbox = QCheckBox("Check Outside\nBand Gap", self)
		self.equilibrium_fermi_energy_checkbox.setChecked(False)
		self.equilibrium_fermi_energy_checkbox.clicked.connect(self.Equilibrium_Fermi_Energy_CheckOutsideBandgap)
		self.equilibrium_fermi_energy_widget_layout.addWidget(self.equilibrium_fermi_energy_checkbox)
		
		
		
		self.tab1_carrierconcentration_widget_layout.addWidget(self.equilibrium_fermi_energy_widget)
	
	
	
	
	def Update_Equilibrium_Fermi_Energy_Temperature(self):
		
		temperature = float(self.temperature_selection_box.currentText())
		intrinsic_equilibrium_fermi_energy = self.CarrierConcentration.intrinsic_equilibrium_fermi_energy[temperature]
		total_equilibrium_fermi_energy = self.CarrierConcentration.total_equilibrium_fermi_energy[temperature]
		
		self.equilibrium_fermi_energy_display.setText(str(total_equilibrium_fermi_energy))
		
		if (str(total_equilibrium_fermi_energy) == "< EVBM") or (str(total_equilibrium_fermi_energy) == "> ECBM"):
			self.equilibrium_fermi_energy_display.setStyleSheet("""QLineEdit { background-color: white; color: red }""")
		else:
			self.equilibrium_fermi_energy_display.setStyleSheet("""QLineEdit { background-color: white; color: black }""")
		
		self.DefectsDiagram.Plot_Equilibrium_Fermi_Energy(temperature=temperature, equilibrium_fermi_energy=total_equilibrium_fermi_energy)
	
	
	
	
	
	def Equilibrium_Fermi_Energy_CheckOutsideBandgap(self):
		if self.equilibrium_fermi_energy_checkbox.isChecked():
			self.CarrierConcentration.check_outside_bandgap = True
		else:
			self.CarrierConcentration.check_outside_bandgap = False
	
	
	
	
	
	
	
	
	
	###############################################################################################
	################################# Generate Defects Diagram ####################################
	###############################################################################################
	
	def Generate_DefectsDiagram_Plot_Function(self, event):
		
		# Update elements and chemical potentials
		self.DefectsDiagram.first_element	= self.first_element
		self.DefectsDiagram.second_element	= self.second_element
		self.DefectsDiagram.mu_elements[self.first_element]["mu0"] = self.compounds_info[self.first_element]["mu0"]
		self.DefectsDiagram.mu_elements[self.second_element]["mu0"] = self.compounds_info[self.second_element]["mu0"]
		self.DefectsDiagram.mu_elements[self.first_element]["deltamu"] = self.deltamu_values[self.first_element]
		self.DefectsDiagram.mu_elements[self.second_element]["deltamu"] = self.deltamu_values[self.second_element]
		
		# Reset defects diagram
		self.DefectsDiagram.defects_diagram_plot_drawing.remove()
		self.DefectsDiagram.defects_diagram_plot_drawing = self.DefectsDiagram.defects_diagram_plot_figure.add_subplot(111)
		self.DefectsDiagram.Activate_DefectsDiagram_Plot_Axes()
		
		# Calculate defect formation energies
		self.DefectsDiagram.Calculate_DefectFormations()
		
		# Plot defect formation energies
		self.DefectsDiagram.intrinsic_defect_plots = {}
		self.DefectsDiagram.extrinsic_defect_plots = {}
		self.DefectsDiagram.Initialize_Intrinsic_DefectsDiagram_Plot()
		
		
		
		
		if self.show_carrier_concentration:
			
			# Update elements and chemical potentials
			self.CarrierConcentration.first_element = self.first_element
			self.CarrierConcentration.second_element = self.second_element
			self.CarrierConcentration.number_species = {self.first_element: self.main_compound_number_first_specie,
														self.second_element: self.main_compound_number_second_specie }
			self.CarrierConcentration.mu_elements[self.first_element]["mu0"] = self.compounds_info[self.first_element]["mu0"]
			self.CarrierConcentration.mu_elements[self.second_element]["mu0"] = self.compounds_info[self.second_element]["mu0"]
			self.CarrierConcentration.mu_elements[self.first_element]["deltamu"] = self.deltamu_values[self.first_element]
			self.CarrierConcentration.mu_elements[self.second_element]["deltamu"] = self.deltamu_values[self.second_element]
			
			# Calculate carrier concentrations
			if self.carrierconcentration_holes_checkbox.isChecked():
				self.CarrierConcentration.Initialize_HoleConcentration_Plot()
		
			if self.carrierconcentration_electrons_checkbox.isChecked():
				self.CarrierConcentration.Initialize_ElectronConcentration_Plot()
			
			# Plot the equilibrium Fermi energy
			self.Update_Equilibrium_Fermi_Energy_Temperature()
	
	
	
	def Update_ExtrinsicDefect_DeltaMu(self):
		
		# Obtain deltamu of dopant
		self.extrinsic_defect_deltamu = float(self.extrinsic_defect_chemical_potential_deltamu.text())
		
		# Recalculate defect formation energies
		self.DefectsDiagram.extrinsic_defect_deltamu = self.extrinsic_defect_deltamu
		self.DefectsDiagram.Calculate_DefectFormations()
		
		# Redraw defects diagram
		if self.DefectsDiagram.intrinsic_defect_plots != {}:
			self.DefectsDiagram.Update_Intrinsic_DefectsDiagram_Plot()
		if self.DefectsDiagram.extrinsic_defect_plots != {}:
			self.DefectsDiagram.Update_Extrinsic_DefectsDiagram_Plot()
		
		
		
		if self.show_carrier_concentration:
			
			# Recalculate carrier concentrations
			self.CarrierConcentration.extrinsic_defect_deltamu = self.extrinsic_defect_deltamu
			
			# Redraw carrier concentration
			if self.CarrierConcentration.carrier_concentration_intrinsic_defect_hole_plot != None:
				self.CarrierConcentration.Update_HoleConcentration_Plot()
			if self.CarrierConcentration.carrier_concentration_intrinsic_defect_electron_plot != None:
				self.CarrierConcentration.Update_ElectronConcentration_Plot()
			
			# Plot the equilibrium Fermi energy
			if (self.DefectsDiagram.intrinsic_defect_plots != {}) and (self.DefectsDiagram.extrinsic_defect_plots != {}):
				self.Update_Equilibrium_Fermi_Energy_Temperature()
	
	
	
	def Update_ExtrinsicDefect(self):
		
		# Obtain selected dopant
		self.extrinsic_defect = self.extrinsic_defect_selection_box.currentText()
		
		# Set intrinsic chemical potential mu0 of dopant
		if self.extrinsic_defect == "None":
			self.extrinsic_defect_chemical_potential_label.setText(u"\u0394"+"\u03BC"+"<sub>x</sub>")
			self.extrinsic_defect_mu0 = 0.0
		else:
			self.extrinsic_defect_chemical_potential_label.setText(u"\u0394"+"\u03BC"+"<sub>"+self.extrinsic_defect.split("_")[0]+"</sub>")
			self.extrinsic_defect_mu0 = self.compounds_info[self.extrinsic_defect.split("_")[0]]["mu0"]
		
		# Reset deltamu of dopant
		self.extrinsic_defect_chemical_potential_deltamu.setText("-0.0000")
		self.extrinsic_defect_deltamu = 0.0
		
		# Recalculate defect formation energies
		self.DefectsDiagram.extrinsic_defect = self.extrinsic_defect
		self.DefectsDiagram.extrinsic_defect_mu0 = self.extrinsic_defect_mu0
		self.DefectsDiagram.extrinsic_defect_deltamu = self.extrinsic_defect_deltamu
		self.DefectsDiagram.Calculate_DefectFormations()
		
		# Draw fresh defects diagram
		self.DefectsDiagram.extrinsic_defect_plots = {}
		self.DefectsDiagram.defects_diagram_plot_drawing.remove()
		self.DefectsDiagram.defects_diagram_plot_drawing = self.DefectsDiagram.defects_diagram_plot_figure.add_subplot(111)
		self.DefectsDiagram.Activate_DefectsDiagram_Plot_Axes()
		self.DefectsDiagram.Initialize_Intrinsic_DefectsDiagram_Plot()
		if self.extrinsic_defect != "None":
			self.DefectsDiagram.Initialize_Extrinsic_DefectsDiagram_Plot()
		
		
		
		if self.show_carrier_concentration:
			
			# Recalculate carrier concentrations
			self.CarrierConcentration.extrinsic_defect = self.extrinsic_defect
			self.CarrierConcentration.extrinsic_defect_mu0 = self.extrinsic_defect_mu0
			self.CarrierConcentration.extrinsic_defect_deltamu = self.extrinsic_defect_deltamu
			
			# Redraw carrier concentration
			if self.CarrierConcentration.carrier_concentration_intrinsic_defect_hole_plot != None:
				self.CarrierConcentration.Initialize_HoleConcentration_Plot()
			if self.CarrierConcentration.carrier_concentration_intrinsic_defect_electron_plot != None:
				self.CarrierConcentration.Initialize_ElectronConcentration_Plot()
			
			# Plot the equilibrium Fermi energy
			if self.DefectsDiagram.intrinsic_defect_plots != {}:
				self.Update_Equilibrium_Fermi_Energy_Temperature()
	
	
	
	def Update_SynthesisTemperature(self):
		
		# Obtain synthesis temperature (K)
		synthesis_temperature = self.defects_synthesis_temperature_box.text()
		
		# Check whether the written synthesis temperature is a possible temperature
		try:
			float(synthesis_temperature)
		except:
			self.defects_synthesis_temperature_box.setText("")
			self.CarrierConcentration.synthesis_temperature = None
			if self.CarrierConcentration.carrier_concentration_intrinsic_defect_hole_plot != None:
				self.CarrierConcentration.Initialize_HoleConcentration_Plot()
			if self.CarrierConcentration.carrier_concentration_intrinsic_defect_electron_plot != None:
				self.CarrierConcentration.Initialize_ElectronConcentration_Plot()
			return
		if float(synthesis_temperature) <= 0.0:
			self.defects_synthesis_temperature_box.setText("")
			self.CarrierConcentration.synthesis_temperature = None
			if self.CarrierConcentration.carrier_concentration_intrinsic_defect_hole_plot != None:
				self.CarrierConcentration.Initialize_HoleConcentration_Plot()
			if self.CarrierConcentration.carrier_concentration_intrinsic_defect_electron_plot != None:
				self.CarrierConcentration.Initialize_ElectronConcentration_Plot()
			return
		
		# Update synthesis temperature in CarrierConcentration object
		self.CarrierConcentration.synthesis_temperature = float(synthesis_temperature)
		
		# Redraw carrier concentration plot with synthesis temperature
		if self.CarrierConcentration.carrier_concentration_intrinsic_defect_hole_plot != None:
			self.CarrierConcentration.Initialize_HoleConcentration_Plot()
		if self.CarrierConcentration.carrier_concentration_intrinsic_defect_electron_plot != None:
			self.CarrierConcentration.Initialize_ElectronConcentration_Plot()
	
	
	
	
	###############################################################################################
	############################# Control Carrier Concentration Plot ##############################
	###############################################################################################
	
	def Generate_CarrierConcentration_Plot_Function(self, event):
		
		self.CarrierConcentration.first_element = self.first_element
		self.CarrierConcentration.second_element = self.second_element
		self.CarrierConcentration.main_compound_number_first_specie = self.main_compound_number_first_specie
		self.CarrierConcentration.main_compound_number_second_specie = self.main_compound_number_second_specie
		self.CarrierConcentration.mu_elements[self.first_element]["mu0"] = self.compounds_info[self.first_element]["mu0"]
		self.CarrierConcentration.mu_elements[self.second_element]["mu0"] = self.compounds_info[self.second_element]["mu0"]
		self.CarrierConcentration.mu_elements[self.first_element]["deltamu"] = self.deltamu_values[self.first_element]
		self.CarrierConcentration.mu_elements[self.second_element]["deltamu"] = self.deltamu_values[self.second_element]
		
		# Reset carrier concentration plot
		self.CarrierConcentration.carrier_concentration_plot_drawing.remove()
		self.CarrierConcentration.carrier_concentration_plot_drawing = self.CarrierConcentration.carrier_concentration_plot_figure.add_subplot(111)
		self.CarrierConcentration.Activate_CarrierConcentration_Plot_Axes()
		
		# Plot the carrier concentration (holes and electrons)
		self.CarrierConcentration.carrier_concentration_intrinsic_defect_hole_plot = None
		self.CarrierConcentration.carrier_concentration_intrinsic_defect_electron_plot = None
		self.CarrierConcentration.carrier_concentration_total_hole_plot = None
		self.CarrierConcentration.carrier_concentration_total_electron_plot = None
		
		if self.carrierconcentration_holes_checkbox.isChecked():
			self.CarrierConcentration.Initialize_HoleConcentration_Plot()
		
		if self.carrierconcentration_electrons_checkbox.isChecked():
			self.CarrierConcentration.Initialize_ElectronConcentration_Plot()
		
		# Plot the equilibrium Fermi energy
		if self.DefectsDiagram.intrinsic_defect_plots != {}:
			self.Update_Equilibrium_Fermi_Energy_Temperature()
	
	
	
	
	
	
	def Update_WindowSize(self, plot_type, ytype):
		
		# Modify defects diagram y-axis
		if plot_type == "DefectsDiagram":
			if ytype == "YMin":
				self.DefectsDiagram.ymin = float(self.defectsdiagram_Ymin_box.text())
			if ytype == "YMax":
				self.DefectsDiagram.ymax = float(self.defectsdiagram_Ymax_box.text())
			self.DefectsDiagram.defects_diagram_plot_drawing.set_ylim(self.DefectsDiagram.ymin, self.DefectsDiagram.ymax)
			self.DefectsDiagram.defects_diagram_plot_canvas.draw()
		
		# Modify carrier concentration y-axis
		elif plot_type == "CarrierConcentration":
			if ytype == "YMin":
				self.CarrierConcentration.ymin = float(self.carrierconcentration_Ymin_box.text())
			if ytype == "YMax":
				self.CarrierConcentration.ymax = float(self.carrierconcentration_Ymax_box.text())
			self.CarrierConcentration.carrier_concentration_plot_drawing.set_ylim(self.CarrierConcentration.ymin, self.CarrierConcentration.ymax)
			self.CarrierConcentration.carrier_concentration_plot_canvas.draw()
	
	
	
	
	###############################################################################################
	###################################### Save Figure ############################################
	###############################################################################################
	
	def SaveFigure(self, figure_type):
		
		options = QFileDialog.Options()
		options |= QFileDialog.DontUseNativeDialog
		filename, extension_type = QFileDialog.getSaveFileName(self, "Save "+figure_type+" Figure", "", "Portable Network Graphics (*.png);;" \
																										+"Portable Document Format (*.pdf);;" \
																										+"Scalable Vector Graphics (*.svg);;" \
																										+"Encapsulated PostScript (*.eps)", options=options)
		if filename:
			extension = extension_type.split(".")[-1].split(")")[0]
			if filename.split(".")[-1] == extension:
				if figure_type == "Phase Diagram":
					self.PhaseDiagram.binary_phase_diagram_plot_figure.savefig(filename, bbox_inches='tight')
				elif figure_type == "Defects Diagram":
					self.DefectsDiagram.defects_diagram_plot_figure.savefig(filename, bbox_inches='tight')
				elif figure_type == "Carrier Concentration":
					self.CarrierConcentration.carrier_concentration_plot_figure.savefig(filename, bbox_inches='tight')
			else:
				if figure_type == "Phase Diagram":
					self.PhaseDiagram.binary_phase_diagram_plot_figure.savefig(filename+"."+extension, bbox_inches='tight')
				elif figure_type == "Defects Diagram":
					self.DefectsDiagram.defects_diagram_plot_figure.savefig(filename+"."+extension, bbox_inches='tight')
				elif figure_type == "Carrier Concentration":
					self.CarrierConcentration.carrier_concentration_plot_figure.savefig(filename+"."+extension, bbox_inches='tight')








