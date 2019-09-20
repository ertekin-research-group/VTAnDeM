
__author__ = 'Michael_Lidia_Jiaxing_Benita_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'







import numpy as np

import PyQt5
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

from vtandem.visualization.quaternary.quaternary_scripts.plot_quaternary_phase_diagram import ChemicalPotential_Quaternary_PhaseDiagramProjected2D
from vtandem.visualization.quaternary.quaternary_scripts.plot_quaternary_defects_diagram import Quaternary_Defects_Diagram
from vtandem.visualization.quaternary.quaternary_scripts.plot_quaternary_carrier_concentration import Quaternary_Carrier_Concentration




class Tab_PhaseDiagram_DefectsDiagram_CarrierConcentration(QWidget):
	
	def __init__(self, parent = None, main_compound = None, first_element = None, second_element = None, third_element = None, fourth_element = None, compounds_info = None, defects_data = None, dos_data = None, show_defects_diagram = True, show_carrier_concentration = True):
		
		QWidget.__init__(self)
		
		
		# Font description for phase stability diagram plot
		self.font = {'family': 'sans-serif', 'color':  'black', 'weight': 'normal', 'size': 16 }
		
		
		# Display variables
		self.show_defects_diagram = show_defects_diagram
		self.show_carrier_concentration = show_carrier_concentration
		
		
		
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
		
		# Obtain DFT data
		self.compounds_info = compounds_info	# Total energies/enthalpies for phase diagram
		self.defects_data = defects_data		# Defect energies for defects diagram
		self.dos_data = dos_data 				# DOS for carrier concentration and equilibrium Fermi energy
		
		# Information about main quaternary compound
		self.main_compound_number_first_specie = self.compounds_info[self.main_compound][self.first_element]	# Number of first species in quaternary compound
		self.main_compound_number_second_specie = self.compounds_info[self.main_compound][self.second_element]	# Number of second species in quaternary compound
		self.main_compound_number_third_specie = self.compounds_info[self.main_compound][self.third_element]	# Number of third species in quaternary compound
		self.main_compound_number_fourth_specie = self.compounds_info[self.main_compound][self.fourth_element]	# Number of fourth species in quaternary compound
		self.main_compound_enthalpy = self.compounds_info[self.main_compound]["enthalpy"]	# Enthalpy of quaternary compound
		
		# Keep track of mu values of the species in the quaternary compound (will be updated as user uses mu4 slider)
		self.deltamu_values = {}
		self.deltamu_values[first_element] = 0.0
		self.deltamu_values[second_element] = 0.0
		self.deltamu_values[third_element] = self.main_compound_enthalpy / self.main_compound_number_third_specie
		self.deltamu_values[fourth_element] = 0.0
		
		# Mu4 values (chemical potential of fourth species)
		self.mu4_value_array = np.linspace(-4.0, -0.0, 1001)	# Array of mu4 values (for the mu4 slider bar)
		
		
		# Extrinsic dopant
		self.extrinsic_defect = "None"
		self.extrinsic_defect_mu0 = 0.0
		self.extrinsic_defect_deltamu = 0.0
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		###############################################################################################
		###############################################################################################
		#################################### Initialize first tab #####################################
		###############################################################################################
		###############################################################################################
		
		
		
		##### Tab 1 (first tab) #####
		self.tab1 = QWidget()
		self.tab1_layout = QHBoxLayout(self.tab1)
		
		
		
		
		
		
		
		
		
		
		
		
		###############################################################################################
		###############################################################################################
		#################### Initialize plot objects (imported from other scripts) ####################
		###############################################################################################
		###############################################################################################
		
		
		
		
		
		# PHASE DIAGRAM
		self.tab1_phasediagram_widget = QWidget()											# Similar to the main widget (defined above), the quaternary
																							#	phase diagram widget will be a "sub-main" widget containing
																							#	all widgets related to the quaternary phase diagram. These 
																							#	include:
																							#		- Quaternary phase diagram plot
																							#		- Mu value displays
																							#		- Slider bar for the fourth mu value
																							#		- "Generate phase diagram!" button
																							#		- "Save figure!" button
		self.tab1_phasediagram_widget_layout = QVBoxLayout(self.tab1_phasediagram_widget)	# The widgets related to the quaternary phase diagram will be 
																							#	stacked vertically.
		
		# Set up quaternary phase diagram object
		self.PhaseDiagram = ChemicalPotential_Quaternary_PhaseDiagramProjected2D(self, main_compound = main_compound, first_element = first_element, second_element = second_element, third_element = third_element, fourth_element = fourth_element)
		self.PhaseDiagram.compounds_info = self.compounds_info
		self.PhaseDiagram.Update_PhaseDiagram_Object()
		self.PhaseDiagram.Update_PhaseDiagram_Plot_Axes()
		
		# Quaternary phase diagram plot
		self.quaternary_phase_diagram_plot = self.PhaseDiagram.quaternary_phase_diagram_plot_canvas
		
		# Plot point that is pressed by user on phase diagram
		self.pressed_point = self.PhaseDiagram.quaternary_phase_diagram_plot_figure.canvas.mpl_connect('button_press_event', self.Pressed_Point)
		self.pressed_point_desc = {'color': 'red', 'marker': 'o'}
		self.pressed_point_plot, = self.PhaseDiagram.quaternary_phase_diagram_plot_drawing.plot([], [], color=self.pressed_point_desc['color'], marker=self.pressed_point_desc['marker'])
		
		# (WIDGET) Quaternary phase diagram plot
		self.tab1_phasediagram_widget_layout.addWidget(self.quaternary_phase_diagram_plot)
		
		# (WIDGET) Mu value settings
		self.muvalue_settings_widget = QWidget()
		self.muvalue_settings_layout = QVBoxLayout(self.muvalue_settings_widget)
		self.tab1_phasediagram_widget_layout.addWidget(self.muvalue_settings_widget)
		self.Activate_MuValue_Settings()
		
		# (WIDGET) Button to generate phase diagram
		self.generate_quaternary_phase_diagram_plot_button_widget = QPushButton("Generate Phase Diagram")
		self.generate_quaternary_phase_diagram_plot_button_widget.clicked[bool].connect(self.Generate_PhaseDiagram_Plot_Function)
		self.tab1_phasediagram_widget_layout.addWidget(self.generate_quaternary_phase_diagram_plot_button_widget)
		
		# (WIDGET) Save phase diagram as figure
		self.phase_diagram_savefigure_button = QPushButton("Save Phase Diagram Figure")
		self.phase_diagram_savefigure_button.clicked[bool].connect(lambda: self.SaveFigure("Phase Diagram"))
		self.tab1_phasediagram_widget_layout.addWidget(self.phase_diagram_savefigure_button)
		
		# Add the phase diagram widget to Tab 1
		self.tab1_layout.addWidget(self.tab1_phasediagram_widget)
		
		
		
		
		
		
		if self.show_defects_diagram:
			
			# DEFECTS DIAGRAM
			self.tab1_defectsdiagram_widget = QWidget()												# A layout similar to that of the quaternary phase diagram
																									# 	will be used for the defects diagram. This "sub-main" 
																									#	widget will contain the following widgets:
																									#		- Defects diagram plot
																									#		- y-axis dialog
																									#		- "Generate defects diagram!" button
																									#		- "Save figure!" button
			self.tab1_defectsdiagram_widget_layout = QVBoxLayout(self.tab1_defectsdiagram_widget)	# Again, the above mentioned widgets will be stacked vertically.
			
			# Set up defects diagram object
			self.DefectsDiagram = Quaternary_Defects_Diagram(self, main_compound = main_compound, first_element = first_element, second_element = second_element, third_element = third_element, fourth_element = fourth_element)
			self.DefectsDiagram.quaternary_defects_data = self.defects_data[self.main_compound]
			self.DefectsDiagram.main_compound_number_first_specie = self.main_compound_number_first_specie
			self.DefectsDiagram.main_compound_number_second_specie = self.main_compound_number_second_specie
			self.DefectsDiagram.main_compound_number_third_specie = self.main_compound_number_third_specie
			self.DefectsDiagram.main_compound_number_fourth_specie = self.main_compound_number_fourth_specie
			self.DefectsDiagram.main_compound_total_energy = self.compounds_info[main_compound]["total_energy"]
			"""
			self.DefectsDiagram.first_element_mu0 = self.compounds_info[first_element]["mu0"]
			self.DefectsDiagram.second_element_mu0 = self.compounds_info[second_element]["mu0"]
			self.DefectsDiagram.third_element_mu0 = self.compounds_info[third_element]["mu0"]
			self.DefectsDiagram.fourth_element_mu0 = self.compounds_info[fourth_element]["mu0"]
			self.DefectsDiagram.mu_values[self.first_element] = self.mu_values[self.first_element]
			self.DefectsDiagram.mu_values[self.second_element] = self.mu_values[self.second_element]
			self.DefectsDiagram.mu_values[self.third_element] = self.mu_values[self.third_element]
			self.DefectsDiagram.mu_values[self.fourth_element] = self.mu_values[self.fourth_element]
			"""
			self.DefectsDiagram.mu_elements[self.first_element]["mu0"] = self.compounds_info[self.first_element]["mu0"]
			self.DefectsDiagram.mu_elements[self.second_element]["mu0"] = self.compounds_info[self.second_element]["mu0"]
			self.DefectsDiagram.mu_elements[self.third_element]["mu0"] = self.compounds_info[self.third_element]["mu0"]
			self.DefectsDiagram.mu_elements[self.fourth_element]["mu0"] = self.compounds_info[self.fourth_element]["mu0"]
			self.DefectsDiagram.mu_elements[self.first_element]["deltamu"] = self.deltamu_values[self.first_element]
			self.DefectsDiagram.mu_elements[self.second_element]["deltamu"] = self.deltamu_values[self.second_element]
			self.DefectsDiagram.mu_elements[self.third_element]["deltamu"] = self.deltamu_values[self.third_element]
			self.DefectsDiagram.mu_elements[self.fourth_element]["deltamu"] = self.deltamu_values[self.fourth_element]
			self.DefectsDiagram.EVBM = self.compounds_info[main_compound]["vbm"]
			self.DefectsDiagram.ECBM = self.DefectsDiagram.EVBM + self.compounds_info[main_compound]["band_gap"]
			self.DefectsDiagram.fermi_energy_array = np.linspace(self.DefectsDiagram.EVBM, self.DefectsDiagram.ECBM, 100)
			self.DefectsDiagram.Activate_DefectsDiagram_Plot_Axes()
			
			# Defects diagram plot
			self.defects_diagram_plot = self.DefectsDiagram.quaternary_defects_diagram_plot_canvas
			
			# (WIDGET) Defects diagram plot
			self.tab1_defectsdiagram_widget_layout.addWidget(self.defects_diagram_plot)
			
			# (WIDGET) Y-axis limits for defects diagram
			self.defectsdiagram_viewport = QWidget()
			self.defectsdiagram_viewport_layout = QHBoxLayout(self.defectsdiagram_viewport)
			
			# (WIDGET) Y-axis limits for defects diagram
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
			
			# (WIDGET) Extrinsic defect properties
			self.extrinsic_defect_properties = QWidget()
			self.extrinsic_defect_properties_layout = QHBoxLayout(self.extrinsic_defect_properties)
			
			# Extrinsic defect chemical potential
			self.extrinsic_defect_chemical_potential_label = QLabel(u"\u0394"+"\u03BC"+"<sub>x</sub> = ")
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
					#self.extrinsic_dopant_selection_box.addItem(u" "+defect.split("_")[0]+"<sub>"+defect.split("_")[-1]+"</sub>")
					#self.extrinsic_dopant_selection_box.addItem(u"\u0394\u03BC"+"<sub>a</sub>"+"<sub>a</sub>")
					#self.extrinsic_dopant_selection_box.addItem(r""+defect.split("_")[0]+"$_{"+defect.split("_")[-1]+"}$")
					#self.extrinsic_dopant_selection_box.addItem(defect.split("_")[0]+'\u2083'+defect.split("_")[-1])
					self.extrinsic_defect_selection_box.addItem(defect)
			self.extrinsic_defect_selection_box.setCurrentIndex(0)
			self.extrinsic_defect_selection_box.activated.connect(self.Update_ExtrinsicDefect)
			self.extrinsic_defect_properties_layout.addWidget(self.extrinsic_defect_selection_box)
			
			self.tab1_defectsdiagram_widget_layout.addWidget(self.extrinsic_defect_properties)
			
			# (WIDGET) Button to generate defects diagram
			self.generate_defects_diagram_plot_button_widget = QPushButton("Generate Defect Diagram")
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
			self.CarrierConcentration = Quaternary_Carrier_Concentration(self, main_compound = main_compound, first_element = first_element, second_element = second_element, third_element = third_element, fourth_element = fourth_element)
			self.CarrierConcentration.quaternary_defects_data = self.defects_data[self.main_compound]
			self.CarrierConcentration.quaternary_dos_data = self.dos_data[self.main_compound]
			self.CarrierConcentration.main_compound_number_first_specie = self.main_compound_number_first_specie
			self.CarrierConcentration.main_compound_number_second_specie = self.main_compound_number_second_specie
			self.CarrierConcentration.main_compound_number_third_specie = self.main_compound_number_third_specie
			self.CarrierConcentration.main_compound_number_fourth_specie = self.main_compound_number_fourth_specie
			self.CarrierConcentration.main_compound_total_energy = self.compounds_info[main_compound]["total_energy"]
			"""
			self.CarrierConcentration.first_element_mu0 = self.compounds_info[first_element]["mu0"]
			self.CarrierConcentration.second_element_mu0 = self.compounds_info[second_element]["mu0"]
			self.CarrierConcentration.third_element_mu0 = self.compounds_info[third_element]["mu0"]
			self.CarrierConcentration.fourth_element_mu0 = self.compounds_info[fourth_element]["mu0"]
			"""
			self.CarrierConcentration.mu_elements[self.first_element]["mu0"] = self.compounds_info[self.first_element]["mu0"]
			self.CarrierConcentration.mu_elements[self.second_element]["mu0"] = self.compounds_info[self.second_element]["mu0"]
			self.CarrierConcentration.mu_elements[self.third_element]["mu0"] = self.compounds_info[self.third_element]["mu0"]
			self.CarrierConcentration.mu_elements[self.fourth_element]["mu0"] = self.compounds_info[self.fourth_element]["mu0"]
			self.CarrierConcentration.mu_elements[self.first_element]["deltamu"] = self.deltamu_values[self.first_element]
			self.CarrierConcentration.mu_elements[self.second_element]["deltamu"] = self.deltamu_values[self.second_element]
			self.CarrierConcentration.mu_elements[self.third_element]["deltamu"] = self.deltamu_values[self.third_element]
			self.CarrierConcentration.mu_elements[self.fourth_element]["deltamu"] = self.deltamu_values[self.fourth_element]
			self.CarrierConcentration.EVBM = self.DefectsDiagram.EVBM
			self.CarrierConcentration.ECBM = self.DefectsDiagram.ECBM
			self.CarrierConcentration.fermi_energy_array = np.linspace(self.CarrierConcentration.EVBM, self.CarrierConcentration.ECBM, 100)
			self.CarrierConcentration.Activate_CarrierConcentration_Plot_Axes()
			self.CarrierConcentration.Organize_DOS_Data()
			self.CarrierConcentration.Extract_Relevant_Energies_DOSs()
			self.CarrierConcentration.Calculate_Hole_Electron_Concentration_Matrices()
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
	############################## Mu Values Display Settings #####################################
	###############################################################################################
	
	def Activate_MuValue_Settings(self):
		
		# Create user-editable input line widgets for the chemical potential (mu) values. These mu values will change based on
		#	where in the quaternary phase diagram plot the user clicks.
		# *Note: The "sub-main" widget that will contain these mu display ports for each species is already created in __init__ as
		#	the variable self.muvalue_display_widget. It's layout is self.muvalue_settings_layout, and this is where each display
		#	port will be placed.
		
		self.mu_display_label_font = QFont("sans-serif", 12)							# Define the font
		
		# First species
		self.mu1_display_widget = QWidget()										# Create a widget for the display port
		self.mu1_display_layout = QHBoxLayout(self.mu1_display_widget)			# Create a vertical layout for the widget (should only contain the label and the mu4 number, stacked on top of each other)
		self.mu1_species_selection_box = QComboBox()							# Create a drop-down menu so the user can select species
		self.mu1_species_selection_box.addItem(self.first_element)				# Add first species to drop-down menu
		self.mu1_species_selection_box.addItem(self.second_element)
		self.mu1_species_selection_box.addItem(self.third_element)
		self.mu1_species_selection_box.addItem(self.fourth_element)
		self.mu1_species_selection_box.setCurrentIndex(0)						# Set drop-down menu option to the first species initially
		#self.mu1_species_selection_box.activated.connect(lambda: self.Update_Species(1))	# When drop-down menu is used, update species list internally
		self.mu1_display_layout.addWidget(self.mu1_species_selection_box)		# Add the drop-down menu to the mu1 display widget
		self.mu1_display_label = QLabel(u"\u0394"+"\u03BC"+"<sub>a</sub> = ")	# Create a label widget to display the "\Delta\mu_1" text above the actual mu1 value (unicode format)
		self.mu1_display_label.setFont(self.mu_display_label_font)				# Set the font
		self.mu1_display_label.setAlignment(Qt.AlignCenter)						# Align the label to center
		self.mu1_display_layout.addWidget(self.mu1_display_label)				# Add the label	to the mu1 display widget
		self.mu1_display = QLineEdit("-0.0000")									# Create the user-editable line that displays the current mu1 value
		self.mu1_display.setMaxLength(7)										# Set the maximum character count of the entry to 7 (up to 4 decimal places)
		self.mu1_display.editingFinished.connect(lambda: self.Update_MuValue_Displays(1))	# When value is edited by user, update mu values internally
		self.mu1_display_layout.addWidget(self.mu1_display)						# Add the mu1 display line
		
		# Second species
		self.mu2_display_widget = QWidget()
		self.mu2_display_layout = QHBoxLayout(self.mu2_display_widget)
		self.mu2_species_selection_box = QComboBox()
		self.mu2_species_selection_box.addItem(self.first_element)
		self.mu2_species_selection_box.addItem(self.second_element)
		self.mu2_species_selection_box.addItem(self.third_element)
		self.mu2_species_selection_box.addItem(self.fourth_element)
		self.mu2_species_selection_box.setCurrentIndex(1)
		#self.mu2_species_selection_box.activated.connect(lambda: self.Update_Species(2))
		self.mu2_display_layout.addWidget(self.mu2_species_selection_box)
		self.mu2_display_label = QLabel(u"\u0394"+"\u03BC"+"<sub>b</sub> = ")
		self.mu2_display_label.setFont(self.mu_display_label_font)
		self.mu2_display_label.setAlignment(Qt.AlignCenter)
		self.mu2_display_layout.addWidget(self.mu2_display_label)
		self.mu2_display = QLineEdit("-0.0000")
		self.mu2_display.setMaxLength(7)
		self.mu2_display.editingFinished.connect(lambda: self.Update_MuValue_Displays(2))
		self.mu2_display_layout.addWidget(self.mu2_display)
		
		# Third species
		self.mu3_display_widget = QWidget()
		self.mu3_display_layout = QHBoxLayout(self.mu3_display_widget)
		self.mu3_species_selection_box = QComboBox()
		self.mu3_species_selection_box.addItem(self.first_element)
		self.mu3_species_selection_box.addItem(self.second_element)
		self.mu3_species_selection_box.addItem(self.third_element)
		self.mu3_species_selection_box.addItem(self.fourth_element)
		self.mu3_species_selection_box.setCurrentIndex(2)
		#self.mu3_species_selection_box.activated.connect(lambda: self.Update_Species(3))
		self.mu3_display_layout.addWidget(self.mu3_species_selection_box)
		self.mu3_display_label = QLabel(u"\u0394"+"\u03BC"+"<sub>c</sub> = ")
		self.mu3_display_label.setFont(self.mu_display_label_font)
		self.mu3_display_label.setAlignment(Qt.AlignCenter)
		self.mu3_display_layout.addWidget(self.mu3_display_label)
		self.mu3_display = QLineEdit("-0.0000")
		self.mu3_display.setMaxLength(7)
		self.mu3_display.setEnabled(False)										# The mu3 value should NOT be editable by the user
		self.mu3_display.setStyleSheet("""QLineEdit { background-color: white; color: black }""")
		self.mu3_display_layout.addWidget(self.mu3_display)
		
		
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
		#self.mu4_species_selection_box.activated.connect(lambda: self.Update_Species(4))# When the user selects a species to be the fourth species, update the species list internally
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
		self.muvalue_settings_layout.addWidget(self.mu1_display_widget)
		self.muvalue_settings_layout.addWidget(self.mu2_display_widget)
		self.muvalue_settings_layout.addWidget(self.mu3_display_widget)
		self.muvalue_settings_layout.addWidget(self.fourth_element_slider_widget)
	
	
	"""
	def Update_Species(self, species_number):
		
		if species_number == 1:
			self.first_element = str(self.mu1_species_selection_box.currentText())
		elif species_number == 2:
			self.second_element = str(self.mu2_species_selection_box.currentText())
		elif species_number == 3:
			self.third_element = str(self.mu3_species_selection_box.currentText())
		elif species_number == 4:
			self.fourth_element = str(self.mu4_species_selection_box.currentText())
		self.elements_list = [self.first_element, self.second_element, self.third_element, self.fourth_element]
	"""
	
	
	def Update_MuValue_Displays(self, display_number):
		
		# Update mu values of first and second elements
		if display_number == 1:
			try:
				float(self.mu1_display.text())
			except:
				self.mu1_display.setText("-0.0000")
				pass
			if float(self.mu1_display.text()) > 0.0:
				self.mu1_display.setText("-0.0000")
			self.deltamu_values[self.first_element] = float(self.mu1_display.text())
		elif display_number == 2:
			try:
				float(self.mu2_display.text())
			except:
				self.mu2_display.setText("-0.0000")
				pass
			if float(self.mu2_display.text()) > 0.0:
				self.mu2_display.setText("-0.0000")
			self.deltamu_values[self.second_element] = float(self.mu2_display.text())
		
		# Update mu value of third element
		self.deltamu_values[self.third_element] = (self.main_compound_enthalpy - self.main_compound_number_first_specie*self.deltamu_values[self.first_element] - self.main_compound_number_second_specie*self.deltamu_values[self.second_element] - self.main_compound_number_fourth_specie*self.deltamu_values[self.fourth_element]) / self.main_compound_number_third_specie
		self.mu3_display.setText(str(self.deltamu_values[self.third_element]))
		
		# Redraw red dot on phase diagram
		if (self.PhaseDiagram.main_compound_plot != None) and (self.PhaseDiagram.competing_compound_plots != {}):
			self.pressed_point_plot.set_data([self.deltamu_values[self.first_element]], [self.deltamu_values[self.second_element]])
			self.PhaseDiagram.quaternary_phase_diagram_plot_canvas.draw()
		
		
		if self.show_defects_diagram:
			
			# Update chemical potentials in defects diagram object
			"""
			self.DefectsDiagram.mu_values[self.first_element] = self.mu_values[self.first_element]
			self.DefectsDiagram.mu_values[self.second_element] = self.mu_values[self.second_element]
			self.DefectsDiagram.mu_values[self.third_element] = self.mu_values[self.third_element]
			"""
			self.DefectsDiagram.mu_elements[self.first_element]["deltamu"] = self.deltamu_values[self.first_element]
			self.DefectsDiagram.mu_elements[self.second_element]["deltamu"] = self.deltamu_values[self.second_element]
			self.DefectsDiagram.mu_elements[self.third_element]["deltamu"] = self.deltamu_values[self.third_element]
			
			# Recalculate defect formation energies
			self.DefectsDiagram.Calculate_DefectFormations()
			
			# Redraw defects diagram
			if self.DefectsDiagram.intrinsic_defect_plots != {}:
				self.DefectsDiagram.Update_Intrinsic_DefectsDiagram_Plot()
			if self.DefectsDiagram.extrinsic_defect_plots != {}:
				self.DefectsDiagram.Update_Extrinsic_DefectsDiagram_Plot()
		
		
		
		if self.show_carrier_concentration:
			
			# Update chemical potentials in carrier concentration object
			"""
			self.CarrierConcentration.mu_values[self.first_element] = self.mu_values[self.first_element]
			self.CarrierConcentration.mu_values[self.second_element] = self.mu_values[self.second_element]
			self.CarrierConcentration.mu_values[self.third_element] = self.mu_values[self.third_element]
			"""
			self.CarrierConcentration.mu_elements[self.first_element]["deltamu"] = self.deltamu_values[self.first_element]
			self.CarrierConcentration.mu_elements[self.second_element]["deltamu"] = self.deltamu_values[self.second_element]
			self.CarrierConcentration.mu_elements[self.third_element]["deltamu"] = self.deltamu_values[self.third_element]
			
			# Recalculate carrier concentrations
			self.CarrierConcentration.Calculate_CarrierConcentration()
			
			# Redraw carrier concentration
			if self.CarrierConcentration.carrier_concentration_intrinsic_defect_hole_plot != None:
				self.CarrierConcentration.Update_HoleConcentration_Plot()
			if self.CarrierConcentration.carrier_concentration_intrinsic_defect_electron_plot != None:
				self.CarrierConcentration.Update_ElectronConcentration_Plot()
		
		
		
		if self.show_defects_diagram and self.show_carrier_concentration:
			
			# Plot the equilibrium Fermi energy
			if self.DefectsDiagram.intrinsic_defect_plots != {}:
				self.Update_Equilibrium_Fermi_Energy_Temperature()
	
	
	
	def Update_Fourth_Species_Slider(self):
		
		# This is how we update the properties of the slider given that the user touches it.
		fourth_element_mu4_value_index = self.fourth_element_slider.value()	# Obtain the (integer) value of the slider when the user uses the slider
		self.deltamu_values[self.fourth_element] = self.mu4_value_array[fourth_element_mu4_value_index]		# Update the mu4 value by 1) using the value of the slider as the index and 
																										#	2) using the index to choose the mu4 value from self.mu4_value_array
		
		self.deltamu_values[self.third_element] = (self.main_compound_enthalpy - self.main_compound_number_first_specie*self.deltamu_values[self.first_element] - self.main_compound_number_second_specie*self.deltamu_values[self.second_element] - self.main_compound_number_fourth_specie*self.deltamu_values[self.fourth_element]) / self.main_compound_number_third_specie
		self.mu3_display.setText("{0:.4f}".format(self.deltamu_values[self.third_element]))		# Every time the mu4 value changes (and the mu1 and mu2 values are held constant), the mu3 value changes
		
		mu4_rounded = round(self.deltamu_values[self.fourth_element], 4)						# Display the updated mu4 value as a float up to four decimal places
		self.fourth_element_slider_value_label.setText("{0:.4f}".format(mu4_rounded))	# Update the text display
		
		
		self.PhaseDiagram.mu4 = self.deltamu_values[self.fourth_element]
		
		
		
		# Update the quaternary phase diagram
		if (self.PhaseDiagram.main_compound_plot != None) and (self.PhaseDiagram.competing_compound_plots != {}):
			self.PhaseDiagram.Plot_PhaseDiagram()
		
		
		
		if self.show_defects_diagram:
			
			# Update chemical potentials in defects diagram object
			"""
			self.DefectsDiagram.mu_values[self.third_element] = self.mu_values[self.third_element]
			self.DefectsDiagram.mu_values[self.fourth_element] = self.mu_values[self.fourth_element]
			"""
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
			"""
			self.CarrierConcentration.mu_values[self.third_element] = self.mu_values[self.third_element]
			self.CarrierConcentration.mu_values[self.fourth_element] = self.mu_values[self.fourth_element]
			"""
			self.CarrierConcentration.mu_elements[self.third_element]["deltamu"] = self.deltamu_values[self.third_element]
			self.CarrierConcentration.mu_elements[self.fourth_element]["deltamu"] = self.deltamu_values[self.fourth_element]
			
			
			# Recalculate carrier concentrations
			self.CarrierConcentration.Calculate_CarrierConcentration()
			
			# Redraw carrier concentration plot
			if self.CarrierConcentration.carrier_concentration_intrinsic_defect_hole_plot != None:
				self.CarrierConcentration.Update_HoleConcentration_Plot()
			if self.CarrierConcentration.carrier_concentration_intrinsic_defect_electron_plot != None:
				self.CarrierConcentration.Update_ElectronConcentration_Plot()
		
		
		
		if self.show_defects_diagram and self.show_carrier_concentration:
			
			# Update the equilibrium Fermi energy
			if self.DefectsDiagram.intrinsic_defect_plots != {}:
				self.Update_Equilibrium_Fermi_Energy_Temperature()
		
		print(self.PhaseDiagram.PSR_vertices)
	
	
	
	
	
	
	
	
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
		
		
		self.equilibrium_fermi_energy_widget_layout.addItem(QSpacerItem(50, 40, QSizePolicy.Expanding, QSizePolicy.Minimum))
		
		
		
		# Label for equilibrium Fermi energy
		equilibrium_fermi_energy_label = QLabel(u"E"+"<sub>f</sub>"+"<sup>eq</sup> = ")
		equilibrium_fermi_energy_label.setAlignment(Qt.AlignCenter)
		self.equilibrium_fermi_energy_widget_layout.addWidget(equilibrium_fermi_energy_label)
		
		# Display for equilibrium Fermi energy
		"""
		self.equilibrium_fermi_energy_display = QLCDNumber()
		self.equilibrium_fermi_energy_display.setStyleSheet("QLCDNumber {color: black}")
		self.equilibrium_fermi_energy_widget_layout.addWidget(self.equilibrium_fermi_energy_display)
		"""
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
		
		#self.CarrierConcentration.Calculate_CarrierConcentration()
		
		temperature = float(self.temperature_selection_box.currentText())
		intrinsic_equilibrium_fermi_energy = self.CarrierConcentration.intrinsic_equilibrium_fermi_energy[temperature]
		total_equilibrium_fermi_energy = self.CarrierConcentration.total_equilibrium_fermi_energy[temperature]
		
		#self.equilibrium_fermi_energy_display.setText(str(round(total_equilibrium_fermi_energy, 5)))
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
	################################## Generate Phase Diagram #####################################
	###############################################################################################
	
	def Generate_PhaseDiagram_Plot_Function(self, event):
		
		# This function specifies what happens when the user clicks the "Generate Plot" button. The 
		#	only thing it needs to check is whether the elements chosen as the first, second, third, 
		#	and fourth species are unique.
		
		if len(self.elements_list) > len(set(self.elements_list)):	# Every time someone clicks the species' buttons, the self.elements_list gets updated.
																	#	The "set" function checks the self.elements_list and omits any that are repeated.
																	#	If any are repeated, then the chosen species are not unique.
			QMessageBox.about(self, "WARNING", "Pick UNIQUE elements!")
			self.fourth_element_slider.setEnabled(False)	# Disable the slider temporarily, until the user chooses unique elements to plot
			return
		
		# If unique, update each species and the list of elements
		self.first_element = str(self.mu1_species_selection_box.currentText())
		self.second_element = str(self.mu2_species_selection_box.currentText())
		self.third_element = str(self.mu3_species_selection_box.currentText())
		self.fourth_element = str(self.mu4_species_selection_box.currentText())
		self.elements_list = [self.first_element, self.second_element, self.third_element, self.fourth_element]
		
		
		
		# Update the number of each species in the main compound based on the change in species
		self.main_compound_number_first_specie = self.compounds_info[self.main_compound][self.first_element]
		self.main_compound_number_second_specie = self.compounds_info[self.main_compound][self.second_element]
		self.main_compound_number_third_specie = self.compounds_info[self.main_compound][self.third_element]
		self.main_compound_number_fourth_specie = self.compounds_info[self.main_compound][self.fourth_element]
		
		# Reset the slide bar
		self.fourth_element_slider.setEnabled(True)
		self.fourth_element_slider_label.setText(u"\u0394\u03BC<sub>"+self.fourth_element+"</sub>")
		endpoint_slidebar = self.main_compound_enthalpy / self.main_compound_number_fourth_specie
		self.mu4_value_array = np.linspace(endpoint_slidebar, -0.0, 1001)
		self.fourth_element_slider.setValue(1000)
		
		# Reset the mu values and their displays
		self.deltamu_values[self.first_element]  = 0.0
		self.deltamu_values[self.second_element] = 0.0
		self.deltamu_values[self.third_element]  = (self.main_compound_enthalpy - self.main_compound_number_first_specie*self.deltamu_values[self.first_element] - self.main_compound_number_second_specie*self.deltamu_values[self.second_element] - self.main_compound_number_fourth_specie*self.deltamu_values[self.fourth_element]) / self.main_compound_number_third_specie	# Every time the mu4 value changes (and the mu1 and mu2 values are held constant), the mu3 value changes
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
		self.PhaseDiagram.Update_PhaseDiagram_Object()
		"""
		self.PhaseDiagram.main_compound_number_first_specie = self.main_compound_number_first_specie
		self.PhaseDiagram.main_compound_number_second_specie = self.main_compound_number_second_specie
		self.PhaseDiagram.main_compound_number_third_specie = self.main_compound_number_third_specie
		self.PhaseDiagram.main_compound_number_fourth_specie = self.main_compound_number_fourth_specie
		"""
		
		# Set up the new plot
		self.PhaseDiagram.quaternary_phase_diagram_plot_drawing.remove()
		self.PhaseDiagram.quaternary_phase_diagram_plot_drawing = self.PhaseDiagram.quaternary_phase_diagram_plot_figure.add_subplot(111)
		
		# Reset the plots of the main and competing compounds
		self.PhaseDiagram.main_compound_plot = None
		self.PhaseDiagram.competing_compound_plots = {}
		
		# Reset the axes of the plot
		self.PhaseDiagram.Update_PhaseDiagram_Plot_Axes()
		
		# Reset clicked point
		self.pressed_point_plot, = self.PhaseDiagram.quaternary_phase_diagram_plot_drawing.plot([], [], color=self.pressed_point_desc['color'], marker=self.pressed_point_desc['marker'])
		
		# Plot the phase stability diagram of the quaternary compound using the new settings
		self.PhaseDiagram.Plot_PhaseDiagram()
		
		
		if self.show_defects_diagram:
			
			# Update global variables in defects diagram object
			self.DefectsDiagram.first_element = self.first_element
			self.DefectsDiagram.second_element = self.second_element
			self.DefectsDiagram.third_element = self.third_element
			self.DefectsDiagram.fourth_element = self.fourth_element
			"""
			self.DefectsDiagram.first_element_mu0 = self.compounds_info[self.first_element]["mu0"]
			self.DefectsDiagram.second_element_mu0 = self.compounds_info[self.second_element]["mu0"]
			self.DefectsDiagram.third_element_mu0 = self.compounds_info[self.third_element]["mu0"]
			self.DefectsDiagram.fourth_element_mu0 = self.compounds_info[self.fourth_element]["mu0"]
			self.DefectsDiagram.mu_values[self.first_element] = self.mu_values[self.first_element]
			self.DefectsDiagram.mu_values[self.second_element] = self.mu_values[self.second_element]
			self.DefectsDiagram.mu_values[self.third_element] = self.mu_values[self.third_element]
			self.DefectsDiagram.mu_values[self.fourth_element] = self.mu_values[self.fourth_element]
			"""
			self.DefectsDiagram.mu_elements[self.first_element]["mu0"] = self.compounds_info[self.first_element]["mu0"]
			self.DefectsDiagram.mu_elements[self.second_element]["mu0"] = self.compounds_info[self.second_element]["mu0"]
			self.DefectsDiagram.mu_elements[self.third_element]["mu0"] = self.compounds_info[self.third_element]["mu0"]
			self.DefectsDiagram.mu_elements[self.fourth_element]["mu0"] = self.compounds_info[self.fourth_element]["mu0"]
			self.DefectsDiagram.mu_elements[self.first_element]["deltamu"] = self.deltamu_values[self.first_element]
			self.DefectsDiagram.mu_elements[self.second_element]["deltamu"] = self.deltamu_values[self.second_element]
			self.DefectsDiagram.mu_elements[self.third_element]["deltamu"] = self.deltamu_values[self.third_element]
			self.DefectsDiagram.mu_elements[self.fourth_element]["deltamu"] = self.deltamu_values[self.fourth_element]
			
			
			# Reset defects diagram
			self.DefectsDiagram.quaternary_defects_diagram_plot_drawing.remove()
			self.DefectsDiagram.quaternary_defects_diagram_plot_drawing = self.DefectsDiagram.quaternary_defects_diagram_plot_figure.add_subplot(111)
			self.DefectsDiagram.Activate_DefectsDiagram_Plot_Axes()
			self.DefectsDiagram.quaternary_defects_diagram_plot_canvas.draw()
		
		
		
		
		if self.show_carrier_concentration:
			
			# Update global variables in carrier concentration object
			self.CarrierConcentration.first_element = self.first_element
			self.CarrierConcentration.second_element = self.second_element
			self.CarrierConcentration.third_element = self.third_element
			self.CarrierConcentration.fourth_element = self.fourth_element
			self.CarrierConcentration.main_compound_number_first_specie = self.main_compound_number_first_specie
			self.CarrierConcentration.main_compound_number_second_specie = self.main_compound_number_second_specie
			self.CarrierConcentration.main_compound_number_third_specie = self.main_compound_number_third_specie
			self.CarrierConcentration.main_compound_number_fourth_specie = self.main_compound_number_fourth_specie
			"""
			self.CarrierConcentration.first_element_mu0 = self.compounds_info[self.first_element]["mu0"]
			self.CarrierConcentration.second_element_mu0 = self.compounds_info[self.second_element]["mu0"]
			self.CarrierConcentration.third_element_mu0 = self.compounds_info[self.third_element]["mu0"]
			self.CarrierConcentration.fourth_element_mu0 = self.compounds_info[self.fourth_element]["mu0"]
			self.CarrierConcentration.mu_values[self.first_element] = self.mu_values[self.first_element]
			self.CarrierConcentration.mu_values[self.second_element] = self.mu_values[self.second_element]
			self.CarrierConcentration.mu_values[self.third_element] = self.mu_values[self.third_element]
			self.CarrierConcentration.mu_values[self.fourth_element] = self.mu_values[self.fourth_element]
			"""
			self.CarrierConcentration.mu_elements[self.first_element]["mu0"] = self.compounds_info[self.first_element]["mu0"]
			self.CarrierConcentration.mu_elements[self.second_element]["mu0"] = self.compounds_info[self.second_element]["mu0"]
			self.CarrierConcentration.mu_elements[self.third_element]["mu0"] = self.compounds_info[self.third_element]["mu0"]
			self.CarrierConcentration.mu_elements[self.fourth_element]["mu0"] = self.compounds_info[self.fourth_element]["mu0"]
			self.CarrierConcentration.mu_elements[self.first_element]["deltamu"] = self.deltamu_values[self.first_element]
			self.CarrierConcentration.mu_elements[self.second_element]["deltamu"] = self.deltamu_values[self.second_element]
			self.CarrierConcentration.mu_elements[self.third_element]["deltamu"] = self.deltamu_values[self.third_element]
			self.CarrierConcentration.mu_elements[self.fourth_element]["deltamu"] = self.deltamu_values[self.fourth_element]
			
			# Reset carrier concentration plot
			self.CarrierConcentration.carrier_concentration_hole_plot = None
			self.CarrierConcentration.carrier_concentration_electron_plot = None
			self.CarrierConcentration.carrier_concentration_plot_drawing.remove()
			self.CarrierConcentration.carrier_concentration_plot_drawing = self.CarrierConcentration.carrier_concentration_plot_figure.add_subplot(111)
			self.CarrierConcentration.Activate_CarrierConcentration_Plot_Axes()
			self.CarrierConcentration.carrier_concentration_plot_canvas.draw()
	
	
	
	
	###############################################################################################
	################################# Generate Defects Diagram ####################################
	###############################################################################################
	
	def Generate_DefectsDiagram_Plot_Function(self, event):
		
		if self.PhaseDiagram.main_compound_plot == None:
			QMessageBox.about(self, "WARNING", "Generate the phase diagram first!")
			return
		
		if len(self.elements_list) > len(set(self.elements_list)):	# Every time someone clicks the species' buttons, the self.elements list gets updated.
																	#	The "set" function checks the self.elements list and omits any that are repeated.
																	#	If any are repeated, then the chosen species are not unique.
			QMessageBox.about(self, "WARNING", "Pick UNIQUE elements!")
			return
		
		# Update elements and chemical potentials
		self.DefectsDiagram.first_element = self.first_element
		self.DefectsDiagram.second_element = self.second_element
		self.DefectsDiagram.third_element = self.third_element
		self.DefectsDiagram.fourth_element = self.fourth_element
		"""
		self.DefectsDiagram.mu_values[self.first_element] = self.mu_values[self.first_element]
		self.DefectsDiagram.mu_values[self.second_element] = self.mu_values[self.second_element]
		self.DefectsDiagram.mu_values[self.third_element] = self.mu_values[self.third_element]
		self.DefectsDiagram.mu_values[self.fourth_element] = self.mu_values[self.fourth_element]
		"""
		self.DefectsDiagram.mu_elements[self.first_element]["mu0"] = self.compounds_info[self.first_element]["mu0"]
		self.DefectsDiagram.mu_elements[self.second_element]["mu0"] = self.compounds_info[self.second_element]["mu0"]
		self.DefectsDiagram.mu_elements[self.third_element]["mu0"] = self.compounds_info[self.third_element]["mu0"]
		self.DefectsDiagram.mu_elements[self.fourth_element]["mu0"] = self.compounds_info[self.fourth_element]["mu0"]
		self.DefectsDiagram.mu_elements[self.first_element]["deltamu"] = self.deltamu_values[self.first_element]
		self.DefectsDiagram.mu_elements[self.second_element]["deltamu"] = self.deltamu_values[self.second_element]
		self.DefectsDiagram.mu_elements[self.third_element]["deltamu"] = self.deltamu_values[self.third_element]
		self.DefectsDiagram.mu_elements[self.fourth_element]["deltamu"] = self.deltamu_values[self.fourth_element]
		
		# Reset defects diagram
		self.DefectsDiagram.quaternary_defects_diagram_plot_drawing.remove()
		self.DefectsDiagram.quaternary_defects_diagram_plot_drawing = self.DefectsDiagram.quaternary_defects_diagram_plot_figure.add_subplot(111)
		self.DefectsDiagram.Activate_DefectsDiagram_Plot_Axes()
		
		# Calculate defect formation energies
		self.DefectsDiagram.Calculate_DefectFormations()
		
		# Plot defect formation energies
		self.DefectsDiagram.intrinsic_defect_plots = {}
		self.DefectsDiagram.extrinsic_defect_plots = {}
		self.DefectsDiagram.Initialize_Intrinsic_DefectsDiagram_Plot()
		
		"""
		# Plot defect formation energies
		self.DefectsDiagram.defect_plots = {}
		if len(self.defects_diagram_Yaxis_box.text()) == 0:
			self.DefectsDiagram.Plot_DefectsDiagram('HSE06')
		else:
			self.DefectsDiagram.ymin = float(self.defects_diagram_Yaxis_box.text().split(',')[0])
			self.DefectsDiagram.ymax = float(self.defects_diagram_Yaxis_box.text().split(',')[1])
			self.DefectsDiagram.Plot_DefectsDiagram('HSE06')
			self.DefectsDiagram.defects_diagram_plot_drawing.set_ylim(self.DefectsDiagram.ymin, self.DefectsDiagram.ymax)
			self.DefectsDiagram.defects_diagram_plot_canvas.draw()
		"""
		
		if self.show_carrier_concentration:
			
			# Update elements and chemical potentials
			self.CarrierConcentration.first_element = self.first_element
			self.CarrierConcentration.second_element = self.second_element
			self.CarrierConcentration.third_element = self.third_element
			self.CarrierConcentration.fourth_element = self.fourth_element
			self.CarrierConcentration.mu_elements[self.first_element]["mu0"] = self.compounds_info[self.first_element]["mu0"]
			self.CarrierConcentration.mu_elements[self.second_element]["mu0"] = self.compounds_info[self.second_element]["mu0"]
			self.CarrierConcentration.mu_elements[self.third_element]["mu0"] = self.compounds_info[self.third_element]["mu0"]
			self.CarrierConcentration.mu_elements[self.fourth_element]["mu0"] = self.compounds_info[self.fourth_element]["mu0"]
			self.CarrierConcentration.mu_elements[self.first_element]["deltamu"] = self.deltamu_values[self.first_element]
			self.CarrierConcentration.mu_elements[self.second_element]["deltamu"] = self.deltamu_values[self.second_element]
			self.CarrierConcentration.mu_elements[self.third_element]["deltamu"] = self.deltamu_values[self.third_element]
			self.CarrierConcentration.mu_elements[self.fourth_element]["deltamu"] = self.deltamu_values[self.fourth_element]
			
			# Calculate carrier concentrations
			self.CarrierConcentration.Calculate_CarrierConcentration()
			
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
			self.CarrierConcentration.Calculate_CarrierConcentration()
			
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
		self.DefectsDiagram.quaternary_defects_diagram_plot_drawing.remove()
		self.DefectsDiagram.quaternary_defects_diagram_plot_drawing = self.DefectsDiagram.quaternary_defects_diagram_plot_figure.add_subplot(111)
		self.DefectsDiagram.Activate_DefectsDiagram_Plot_Axes()
		self.DefectsDiagram.Initialize_Intrinsic_DefectsDiagram_Plot()
		if self.extrinsic_defect != "None":
			self.DefectsDiagram.Initialize_Extrinsic_DefectsDiagram_Plot()
		
		
		if self.show_carrier_concentration:
			
			# Recalculate carrier concentrations
			self.CarrierConcentration.extrinsic_defect = self.extrinsic_defect
			self.CarrierConcentration.extrinsic_defect_mu0 = self.extrinsic_defect_mu0
			self.CarrierConcentration.extrinsic_defect_deltamu = self.extrinsic_defect_deltamu
			self.CarrierConcentration.Calculate_CarrierConcentration()
			
			# Redraw carrier concentration
			if self.CarrierConcentration.carrier_concentration_intrinsic_defect_hole_plot != None:
				self.CarrierConcentration.Initialize_HoleConcentration_Plot()
			if self.CarrierConcentration.carrier_concentration_intrinsic_defect_electron_plot != None:
				self.CarrierConcentration.Initialize_ElectronConcentration_Plot()
			
			# Plot the equilibrium Fermi energy
			if self.DefectsDiagram.intrinsic_defect_plots != {}:
				self.Update_Equilibrium_Fermi_Energy_Temperature()
	
	
	
	
	
	
	###############################################################################################
	############################## Control Carrier Concentration Plot #############################
	###############################################################################################
	
	def Generate_CarrierConcentration_Plot_Function(self, event):
		
		if self.PhaseDiagram.main_compound_plot == None:
			QMessageBox.about(self, "WARNING", "Generate the phase diagram first!")
		
		if len(self.elements_list) > len(set(self.elements_list)):	# Every time someone clicks the species' buttons, the self.elements list gets updated.
																	#	The "set" function checks the self.elements list and omits any that are repeated.
																	#	If any are repeated, then the chosen species are not unique.
			QMessageBox.about(self, "WARNING", "Pick UNIQUE elements!")
			return
		
		self.CarrierConcentration.first_element = self.first_element
		self.CarrierConcentration.second_element = self.second_element
		self.CarrierConcentration.third_element = self.third_element
		self.CarrierConcentration.fourth_element = self.fourth_element
		self.CarrierConcentration.main_compound_number_first_specie = self.main_compound_number_first_specie
		self.CarrierConcentration.main_compound_number_second_specie = self.main_compound_number_second_specie
		self.CarrierConcentration.main_compound_number_third_specie = self.main_compound_number_third_specie
		self.CarrierConcentration.main_compound_number_fourth_specie = self.main_compound_number_fourth_specie
		"""
		self.CarrierConcentration.mu_values[self.first_element] = self.mu_values[self.first_element]
		self.CarrierConcentration.mu_values[self.second_element] = self.mu_values[self.second_element]
		self.CarrierConcentration.mu_values[self.third_element] = self.mu_values[self.third_element]
		self.CarrierConcentration.mu_values[self.fourth_element] = self.mu_values[self.fourth_element]
		"""
		self.CarrierConcentration.mu_elements[self.first_element]["mu0"] = self.compounds_info[self.first_element]["mu0"]
		self.CarrierConcentration.mu_elements[self.second_element]["mu0"] = self.compounds_info[self.second_element]["mu0"]
		self.CarrierConcentration.mu_elements[self.third_element]["mu0"] = self.compounds_info[self.third_element]["mu0"]
		self.CarrierConcentration.mu_elements[self.fourth_element]["mu0"] = self.compounds_info[self.fourth_element]["mu0"]
		self.CarrierConcentration.mu_elements[self.first_element]["deltamu"] = self.deltamu_values[self.first_element]
		self.CarrierConcentration.mu_elements[self.second_element]["deltamu"] = self.deltamu_values[self.second_element]
		self.CarrierConcentration.mu_elements[self.third_element]["deltamu"] = self.deltamu_values[self.third_element]
		self.CarrierConcentration.mu_elements[self.fourth_element]["deltamu"] = self.deltamu_values[self.fourth_element]
		
		
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
	
	
	
	
	###############################################################################################
	################################ Clicking on Phase Diagram ####################################
	###############################################################################################
	
	def Pressed_Point(self, event):
		
		# This function reads in the coordinates of the point in the phase stability diagram plot that the user clicks.

		point_x = event.xdata	# x-coordinate
		point_y = event.ydata	# y-coordinate
		
		if (not isinstance(point_x, float)) or (not isinstance(point_y, float)):	# Check that the user presses somewhere on the plot (and not anywhere else)
			return
		
		elif self.PhaseDiagram.main_compound_plot == None:
			return
		
		
		# Update the clicked mu values
		self.deltamu_values[self.first_element] = point_x
		self.deltamu_values[self.second_element] = point_y
		self.deltamu_values[self.third_element] = (self.main_compound_enthalpy - self.main_compound_number_first_specie*self.deltamu_values[self.first_element] - self.main_compound_number_second_specie*self.deltamu_values[self.second_element] - self.main_compound_number_fourth_specie*self.deltamu_values[self.fourth_element]) / self.main_compound_number_third_specie
		
		# Update the display ports for the mu values
		self.mu1_display.setText("{0:.4f}".format(self.deltamu_values[self.first_element]))
		self.mu2_display.setText("{0:.4f}".format(self.deltamu_values[self.second_element]))
		self.mu3_display.setText("{0:.4f}".format(self.deltamu_values[self.third_element]))
		
		# Update the red dot where the user clicked on the phase diagram
		self.pressed_point_plot.set_data([self.deltamu_values[self.first_element]], [self.deltamu_values[self.second_element]])
		self.PhaseDiagram.quaternary_phase_diagram_plot_figure.canvas.draw_idle()
		
		
		
		if self.show_defects_diagram:
			
			# Update chemical potentials in defects diagram
			"""
			self.DefectsDiagram.mu_values[self.first_element] = self.mu_values[self.first_element]
			self.DefectsDiagram.mu_values[self.second_element] = self.mu_values[self.second_element]
			self.DefectsDiagram.mu_values[self.third_element] = self.mu_values[self.third_element]
			"""
			self.DefectsDiagram.mu_elements[self.first_element]["deltamu"] = self.deltamu_values[self.first_element]
			self.DefectsDiagram.mu_elements[self.second_element]["deltamu"] = self.deltamu_values[self.second_element]
			self.DefectsDiagram.mu_elements[self.third_element]["deltamu"] = self.deltamu_values[self.third_element]
			
			
			# Recalculate defect formation energies
			self.DefectsDiagram.Calculate_DefectFormations()
			
			# Update defects diagram plot
			if self.DefectsDiagram.intrinsic_defect_plots != {}:
				self.DefectsDiagram.Update_Intrinsic_DefectsDiagram_Plot()
			if self.DefectsDiagram.extrinsic_defect_plots != {}:
				self.DefectsDiagram.Update_Extrinsic_DefectsDiagram_Plot()
		
		
		
		if self.show_carrier_concentration:
			
			# Update chemical potentials in carrier concentration
			"""
			self.CarrierConcentration.mu_values[self.first_element] = self.mu_values[self.first_element]
			self.CarrierConcentration.mu_values[self.second_element] = self.mu_values[self.second_element]
			self.CarrierConcentration.mu_values[self.third_element] = self.mu_values[self.third_element]
			"""
			self.CarrierConcentration.mu_elements[self.first_element]["deltamu"] = self.deltamu_values[self.first_element]
			self.CarrierConcentration.mu_elements[self.second_element]["deltamu"] = self.deltamu_values[self.second_element]
			self.CarrierConcentration.mu_elements[self.third_element]["deltamu"] = self.deltamu_values[self.third_element]
			
			# Recalculate carrier concentrations
			self.CarrierConcentration.Calculate_CarrierConcentration()
			
			# Update carrier concentration plot
			if self.CarrierConcentration.carrier_concentration_intrinsic_defect_hole_plot != None:
				self.CarrierConcentration.Update_HoleConcentration_Plot()
			if self.CarrierConcentration.carrier_concentration_intrinsic_defect_electron_plot != None:
				self.CarrierConcentration.Update_ElectronConcentration_Plot()
		
		
		if self.show_defects_diagram and self.show_carrier_concentration:
			
			# Update the equilibrium Fermi energy
			if self.DefectsDiagram.intrinsic_defect_plots != {}:
				self.Update_Equilibrium_Fermi_Energy_Temperature()
	
	
	
	
	
	def Update_WindowSize(self, plot_type, ytype):
		
		# Modify defects diagram y-axis
		if plot_type == "DefectsDiagram":
			if ytype == "YMin":
				self.DefectsDiagram.ymin = float(self.defectsdiagram_Ymin_box.text())
			if ytype == "YMax":
				self.DefectsDiagram.ymax = float(self.defectsdiagram_Ymax_box.text())
			self.DefectsDiagram.quaternary_defects_diagram_plot_drawing.set_ylim(self.DefectsDiagram.ymin, self.DefectsDiagram.ymax)
			self.DefectsDiagram.quaternary_defects_diagram_plot_canvas.draw()
		
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
					self.PhaseDiagram.quaternary_phase_diagram_plot_figure.savefig(filename, bbox_inches='tight')
				elif figure_type == "Defects Diagram":
					self.DefectsDiagram.quaternary_defects_diagram_plot_figure.savefig(filename, bbox_inches='tight')
				elif figure_type == "Carrier Concentration":
					self.CarrierConcentration.carrier_concentration_plot_figure.savefig(filename, bbox_inches='tight')
				elif figure_type == "Phase Diagram 3D":
					self.PhaseDiagram3D.quaternary_phasediagram_3d_plot_figure.savefig(filename, bbox_inches='tight')
			else:
				if figure_type == "Phase Diagram":
					self.PhaseDiagram.quaternary_phase_diagram_plot_figure.savefig(filename+"."+extension, bbox_inches='tight')
				elif figure_type == "Defects Diagram":
					self.DefectsDiagram.quaternary_defects_diagram_plot_figure.savefig(filename+"."+extension, bbox_inches='tight')
				elif figure_type == "Carrier Concentration":
					self.CarrierConcentration.carrier_concentration_plot_figure.savefig(filename+"."+extension, bbox_inches='tight')
				elif figure_type == "Phase Diagram 3D":
					self.PhaseDiagram3D.quaternary_phasediagram_3d_plot_figure.savefig(filename+"."+extension, bbox_inches='tight')




















