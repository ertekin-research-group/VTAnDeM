
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'
__author__ = 'Michael_Lidia_Jiaxing_Elif'

import numpy as np

import PyQt5
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

from vtandem.visualization.windows.window_defectsdiagram import Window_DefectsDiagram
from vtandem.visualization.windows.window_carrierconcentration import Window_CarrierConcentration



#class Tab_PhaseDiagram_DefectsDiagram_CarrierConcentration(QWidget, Window_DefectsDiagram):
class Tab_PhaseDiagram_DefectsDiagram_CarrierConcentration(Window_DefectsDiagram, Window_CarrierConcentration):
	
	def __init__(self, 	compounds_info = None, \
						defects_data = None, \
						main_compound_info = None, \
						dos_data = None, \
						show_defects_diagram = None, \
						show_carrier_concentration = None, \
						type = None	):
		
		
		# Make this tab object into widget so all button connections work!!
		QWidget.__init__(self)
		
		
		
		# Check legitimacy of 'type' argument
		if (type != "ternary") and (type != "quaternary"):
			raise Exception("'type' variable must be either 'ternary' or 'quaternary'. Exiting...")
		self.type = type
		
		
		
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
		
		
		# Find endpoints of phase diagram
		endpoint_candidates = []
		for element in self.elements_list:
			endpoint_candidates.append( self.main_compound_enthalpy/self.main_compound_info["dft_"+element] )
		self.phasediagram_endpoints = min( endpoint_candidates )
		
		
		
		# Keep track of mu values of the species in the quaternary compound (will be updated as user uses mu4 slider)
		self.deltamu_values = {}
		for index, element in enumerate(self.elements_list):
			if index == 2:	# Give nonzero deltamu to third element only
				#self.deltamu_values[element] = self.main_compound_enthalpy / self.main_compound_elements_count[self.third_element]
				self.deltamu_values[element] = self.main_compound_enthalpy / self.main_compound_info["dft_"+self.third_element]
			else:
				self.deltamu_values[element] = 0.0
		
		
		
		
		
		"""
		# Extrinsic dopant
		self.dopant = "None"
		self.dopant_mu0 = 0.0
		self.dopant_deltamu = 0.0
		self.extrinsic_defects = []
		"""
		
		
		
		###############################################################################################
		###############################################################################################
		#################################### Initialize first tab #####################################
		###############################################################################################
		###############################################################################################
		
		self.tab1 = QWidget()
		self.tab1_layout = QHBoxLayout(self.tab1)
		
		
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
		
		
		# Set up ternary phase diagram object
		# Note: We create the PhaseDiagram object in tab_quaternary.../tab_ternary...
		self.PhaseDiagram.main_compound_enthalpy = self.main_compound_enthalpy
		self.PhaseDiagram.phasediagram_endpoints = self.phasediagram_endpoints
		self.PhaseDiagram.compounds_info = self.compounds_info
		self.PhaseDiagram.main_compound_info = self.main_compound_info
		self.PhaseDiagram.Update_PhaseDiagram_Plot_Axes()
		
		
		# (WIDGET) Phase diagram plot
		self.phase_diagram_plot = self.PhaseDiagram.phase_diagram_plot_canvas
		self.tab1_phasediagram_widget_layout.addWidget(self.phase_diagram_plot)
		
		
		# Plot point that is pressed by user on phase diagram
		self.pressed_point = self.PhaseDiagram.phase_diagram_plot_figure.canvas.mpl_connect('button_press_event', self.Pressed_Point)
		self.pressed_point_desc = {'color': 'red', 'marker': 'o'}
		self.pressed_point_plot, = self.PhaseDiagram.phase_diagram_plot_drawing.plot([], [], color=self.pressed_point_desc['color'], marker=self.pressed_point_desc['marker'])
		
		
		# (WIDGET) Mu value settings
		self.muvalue_settings_widget = QWidget()
		self.muvalue_settings_layout = QVBoxLayout(self.muvalue_settings_widget)
		self.tab1_phasediagram_widget_layout.addWidget(self.muvalue_settings_widget)
		self.Activate_MuValue_Settings()
		
		# Activate fourth element slider if main compound is quaternary
		if self.type == "quaternary":
			self.Activate_MuValue_FourthElement_Settings()
		
		
		# (WIDGET) Button to generate phase diagram
		self.generate_phase_diagram_plot_button_widget = QPushButton("Generate Phase Diagram")
		self.generate_phase_diagram_plot_button_widget.clicked[bool].connect(self.Generate_PhaseDiagram_Plot_Function)
		self.tab1_phasediagram_widget_layout.addWidget(self.generate_phase_diagram_plot_button_widget)
		
		# (WIDGET) Save phase diagram as figure
		self.phase_diagram_savefigure_button = QPushButton("Save Phase Diagram Figure")
		self.phase_diagram_savefigure_button.clicked[bool].connect(lambda: self.PhaseDiagram.SaveFigure())
		self.tab1_phasediagram_widget_layout.addWidget(self.phase_diagram_savefigure_button)
		
		# Add the phase diagram widget to Tab 1
		self.tab1_layout.addWidget(self.tab1_phasediagram_widget)
		
		
		
		
		
		if self.show_defects_diagram:
			
			# DEFECTS DIAGRAM
			
			# Set up defects diagram object
			# Note: We create the DefectsDiagram object in tab_quaternary.../tab_ternary...
			self.DefectsDiagram.defects_data = defects_data
			self.DefectsDiagram.main_compound_info = main_compound_info
			
			self.DefectsDiagram.mu_elements[self.first_element]["mu0"] = self.compounds_info[self.first_element]["mu0"]
			self.DefectsDiagram.mu_elements[self.second_element]["mu0"] = self.compounds_info[self.second_element]["mu0"]
			self.DefectsDiagram.mu_elements[self.third_element]["mu0"] = self.compounds_info[self.third_element]["mu0"]
			
			self.DefectsDiagram.mu_elements[self.first_element]["deltamu"] = self.deltamu_values[self.first_element]
			self.DefectsDiagram.mu_elements[self.second_element]["deltamu"] = self.deltamu_values[self.second_element]
			self.DefectsDiagram.mu_elements[self.third_element]["deltamu"] = self.deltamu_values[self.third_element]
			
			if self.type == "quaternary":
				self.DefectsDiagram.mu_elements[self.fourth_element]["mu0"] = self.compounds_info[self.fourth_element]["mu0"]
				self.DefectsDiagram.mu_elements[self.fourth_element]["deltamu"] = self.deltamu_values[self.fourth_element]
			
			self.DefectsDiagram.EVBM = main_compound_info["VBM"]
			self.DefectsDiagram.ECBM = self.DefectsDiagram.EVBM + main_compound_info["BandGap"]
			self.DefectsDiagram.fermi_energy_array = np.linspace(self.DefectsDiagram.EVBM, self.DefectsDiagram.ECBM, 100)
			self.DefectsDiagram.Activate_DefectsDiagram_Plot_Axes()
			
			
			Window_DefectsDiagram.__init__(self, show_dopant = True)
			
			
			# Add the defects diagram widget to Tab 1
			#self.tab1_layout.addWidget(self.defectsdiagram_window.defectsdiagram_window)
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
			# Note: We create the CarrierConcentration object in tab_quaternary.../tab_ternary...
			#self.defectsdiagram_window.CarrierConcentration_Plot = self.CarrierConcentration
			
			self.CarrierConcentration.defects_data = defects_data
			self.CarrierConcentration.main_compound_info = main_compound_info
			self.CarrierConcentration.dos_data = dos_data[self.main_compound]
			
			self.CarrierConcentration.mu_elements[self.first_element]["mu0"] = self.compounds_info[self.first_element]["mu0"]
			self.CarrierConcentration.mu_elements[self.second_element]["mu0"] = self.compounds_info[self.second_element]["mu0"]
			self.CarrierConcentration.mu_elements[self.third_element]["mu0"] = self.compounds_info[self.third_element]["mu0"]
			
			self.CarrierConcentration.mu_elements[self.first_element]["deltamu"] = self.deltamu_values[self.first_element]
			self.CarrierConcentration.mu_elements[self.second_element]["deltamu"] = self.deltamu_values[self.second_element]
			self.CarrierConcentration.mu_elements[self.third_element]["deltamu"] = self.deltamu_values[self.third_element]
			
			if self.type == "quaternary":
				#self.CarrierConcentration.main_compound_number_fourth_specie = self.main_compound_number_fourth_specie
				#self.CarrierConcentration.number_species[self.fourth_element] = self.main_compound_number_fourth_specie
				self.CarrierConcentration.mu_elements[self.fourth_element]["mu0"] = self.compounds_info[self.fourth_element]["mu0"]
				self.CarrierConcentration.mu_elements[self.fourth_element]["deltamu"] = self.deltamu_values[self.fourth_element]
			
			
			self.CarrierConcentration.vol = main_compound_info["Volume"]
			self.CarrierConcentration.EVBM = self.DefectsDiagram.EVBM
			self.CarrierConcentration.ECBM = self.DefectsDiagram.ECBM
			self.CarrierConcentration.fermi_energy_array = np.linspace(self.CarrierConcentration.EVBM, self.CarrierConcentration.ECBM, 100)
			self.CarrierConcentration.Activate_CarrierConcentration_Plot_Axes()
			self.CarrierConcentration.Organize_DOS_Data()
			self.CarrierConcentration.Extract_Relevant_Energies_DOSs()
			self.CarrierConcentration.Calculate_Hole_Electron_Concentration_Matrices()
			
			
			"""
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
			
			
			# (WIDGET) Save carrier concentration plot as figure
			self.carrier_concentration_savefigure_button = QPushButton("Save Carrier Concentration Plot")
			self.carrier_concentration_savefigure_button.clicked[bool].connect(lambda: self.CarrierConcentration.SaveFigure())
			self.tab1_carrierconcentration_widget_layout.addWidget(self.carrier_concentration_savefigure_button)
			"""
			Window_CarrierConcentration.__init__(self)
			
			
			#self.tab1_layout.addWidget(self.tab1_carrierconcentration_widget)
			self.tab1_layout.addWidget(self.carrierconcentration_window)
	
	
	
	
	###############################################################################################
	################################## Mu Values Displays #########################################
	###############################################################################################
	
	def Activate_MuValue_Settings(self):
		
		# Create widgets for the display ports of the mu values. These mu values will change based on where on the ternary phase
		#	diagram plot the user clicks.
		# *Note: The "sub-main" widget that will contain these mu display ports for each species is already created in __init__ as
		#	the variable self.muvalue_display_widget. It's layout is self.muvalue_display_layout, and this is where each display
		#	port will be placed.
		
		self.mu_display_label_font = QFont("Arial", 12)						# Define the font
		
		self.mu1_display_widget = QWidget()									# Create a widget for the display port
		self.mu1_display_layout = QHBoxLayout(self.mu1_display_widget)		# Create a vertical layout for the widget (should only contain the label and the mu4 number, stacked on top of each other)
		self.mu1_species_selection_box = QComboBox()						# Create a drop-down menu so the user can select species
		self.mu1_species_selection_box.addItem(self.first_element)			# Add first species to drop-down menu
		self.mu1_species_selection_box.addItem(self.second_element)
		self.mu1_species_selection_box.addItem(self.third_element)
		self.mu1_species_selection_box.setCurrentIndex(0)					# Set drop-down menu option to the first species initially
		self.mu1_display_layout.addWidget(self.mu1_species_selection_box)	# Add the drop-down menu to the mu1 display widget
		self.mu1_display_label = QLabel(u"\u0394"+"\u03BC"+"<sub>a</sub>")	# Create a label widget to display the "\Delta\mu_1" text above the actual mu1 value (unicode format)
		self.mu1_display_label.setFont(self.mu_display_label_font)			# Set the font
		self.mu1_display_label.setAlignment(Qt.AlignCenter)					# Align the label to center
		self.mu1_display_layout.addWidget(self.mu1_display_label)			# Add the label
		self.mu1_display = QLineEdit("-0.0000")								# Create the actual port that will display the current mu1 value
		self.mu1_display.setMaxLength(7)									# Set the maximum character count of the entry to 7 (up to 4 decimal places)
		self.mu1_display.editingFinished.connect(lambda: self.Update_MuValue_Displays(1))
		self.mu1_display_layout.addWidget(self.mu1_display)					# Add the display port
		
		self.mu2_display_widget = QWidget()
		self.mu2_display_layout = QHBoxLayout(self.mu2_display_widget)
		self.mu2_species_selection_box = QComboBox()
		self.mu2_species_selection_box.addItem(self.first_element)
		self.mu2_species_selection_box.addItem(self.second_element)
		self.mu2_species_selection_box.addItem(self.third_element)
		self.mu2_species_selection_box.setCurrentIndex(1)
		self.mu2_display_layout.addWidget(self.mu2_species_selection_box)
		self.mu2_display_label = QLabel(u"\u0394"+"\u03BC"+"<sub>b</sub>")
		self.mu2_display_label.setFont(self.mu_display_label_font)
		self.mu2_display_label.setAlignment(Qt.AlignCenter)
		self.mu2_display_layout.addWidget(self.mu2_display_label)
		self.mu2_display = QLineEdit("-0.0000")
		self.mu2_display.setMaxLength(7)
		self.mu2_display.editingFinished.connect(lambda: self.Update_MuValue_Displays(2))
		self.mu2_display_layout.addWidget(self.mu2_display)
		
		self.mu3_display_widget = QWidget()
		self.mu3_display_layout = QHBoxLayout(self.mu3_display_widget)
		self.mu3_species_selection_box = QComboBox()
		self.mu3_species_selection_box.addItem(self.first_element)
		self.mu3_species_selection_box.addItem(self.second_element)
		self.mu3_species_selection_box.addItem(self.third_element)
		self.mu3_species_selection_box.setCurrentIndex(2)
		self.mu3_display_layout.addWidget(self.mu3_species_selection_box)
		self.mu3_display_label = QLabel(u"\u0394"+"\u03BC"+"<sub>c</sub>")
		self.mu3_display_label.setFont(self.mu_display_label_font)
		self.mu3_display_label.setAlignment(Qt.AlignCenter)
		self.mu3_display_layout.addWidget(self.mu3_display_label)
		self.mu3_display = QLineEdit("-0.0000")
		self.mu3_display.setMaxLength(7)
		self.mu3_display.setEnabled(False)
		self.mu3_display.setStyleSheet("""QLineEdit { background-color: white; color: black }""")
		self.mu3_display_layout.addWidget(self.mu3_display)
		
		if self.type == "quaternary":
			self.mu1_species_selection_box.addItem(self.fourth_element)
			self.mu2_species_selection_box.addItem(self.fourth_element)
			self.mu3_species_selection_box.addItem(self.fourth_element)
		
		self.muvalue_settings_layout.addWidget(self.mu1_display_widget)
		self.muvalue_settings_layout.addWidget(self.mu2_display_widget)
		self.muvalue_settings_layout.addWidget(self.mu3_display_widget)
	
	
	
	
	
	
	
	
	
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
		if self.type == "ternary":
			#self.deltamu_values[self.third_element] = (self.main_compound_enthalpy - self.main_compound_number_first_specie*self.deltamu_values[self.first_element] - self.main_compound_number_second_specie*self.deltamu_values[self.second_element]) / self.main_compound_number_third_specie
			self.deltamu_values[self.third_element] = (	self.main_compound_enthalpy \
														- self.main_compound_info["dft_"+self.first_element]*self.deltamu_values[self.first_element] \
														- self.main_compound_info["dft_"+self.second_element]*self.deltamu_values[self.second_element] \
														) / self.main_compound_info["dft_"+self.third_element]
		elif self.type == "quaternary":
			#self.deltamu_values[self.third_element] = (self.main_compound_enthalpy - self.main_compound_number_first_specie*self.deltamu_values[self.first_element] - self.main_compound_number_second_specie*self.deltamu_values[self.second_element] - self.main_compound_number_fourth_specie*self.deltamu_values[self.fourth_element]) / self.main_compound_number_third_specie
			self.deltamu_values[self.third_element] = (	self.main_compound_enthalpy \
														- self.main_compound_info["dft_"+self.first_element]*self.deltamu_values[self.first_element] \
														- self.main_compound_info["dft_"+self.second_element]*self.deltamu_values[self.second_element] \
														- self.main_compound_info["dft_"+self.fourth_element]*self.deltamu_values[self.fourth_element] \
														) / self.main_compound_info["dft_"+self.third_element]
		self.mu3_display.setText(str(self.deltamu_values[self.third_element]))
		
		
		
		
		# Redraw red dot on phase diagram
		if (self.PhaseDiagram.main_compound_plot != None) and (self.PhaseDiagram.competing_compound_plots != {}):
			self.pressed_point_plot.set_data([self.deltamu_values[self.first_element]], [self.deltamu_values[self.second_element]])
			self.PhaseDiagram.phase_diagram_plot_canvas.draw()
		
		
		
		if self.show_defects_diagram:
			
			# Update chemical potentials in defects diagram object
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
			self.CarrierConcentration.mu_elements[self.first_element]["deltamu"] = self.deltamu_values[self.first_element]
			self.CarrierConcentration.mu_elements[self.second_element]["deltamu"] = self.deltamu_values[self.second_element]
			self.CarrierConcentration.mu_elements[self.third_element]["deltamu"] = self.deltamu_values[self.third_element]
			
			# Redraw carrier concentration plot
			if self.CarrierConcentration.carrier_concentration_intrinsic_defect_hole_plot != None:
				self.CarrierConcentration.Update_HoleConcentration_Plot()
			if self.CarrierConcentration.carrier_concentration_intrinsic_defect_electron_plot != None:
				self.CarrierConcentration.Update_ElectronConcentration_Plot()
		
		
		
		if self.show_defects_diagram and self.show_carrier_concentration:
			
			# Update the equilibrium Fermi energy
			if self.DefectsDiagram.intrinsic_defect_plots != {}:
				self.Update_Equilibrium_Fermi_Energy_Temperature()
	
	
	
	
	
	
	
	
	
	###############################################################################################
	################################## Generate Phase Diagram #####################################
	###############################################################################################
	
	def Generate_PhaseDiagram_Plot_Function(self, event):
		
		# This function specifies what happens when the user clicks the "Generate Phase Diagram" button.
		
		self.Update_PhaseDiagram_Plot_Elements_Function()	# Found in tab_quaternary.../tab_ternary... script
		
		# Check whether the chosen elements are unique initially
		if len(self.elements_list) > len(set(self.elements_list)):	# Every time someone clicks the species' buttons, the self.elements_list gets updated.
																	#	The "set" function checks the self.elements_list and omits any that are repeated.
																	#	If any are repeated, then the chosen species are not unique.
			QMessageBox.about(self, "WARNING", "Pick UNIQUE elements!")
			self.fourth_element_slider.setEnabled(False)	# Disable the slider temporarily, until the user chooses unique elements to plot
			return
		
		
		# Set up the new plot
		self.PhaseDiagram.phase_diagram_plot_drawing.remove()
		self.PhaseDiagram.phase_diagram_plot_drawing = self.PhaseDiagram.phase_diagram_plot_figure.add_subplot(111)
		
		# Reset the plots of the main and competing compounds
		self.PhaseDiagram.main_compound_plot = None
		self.PhaseDiagram.competing_compound_plots = {}
		
		# Reset the axes of the plot
		self.PhaseDiagram.Update_PhaseDiagram_Plot_Axes()
		
		# Reset clicked point
		self.pressed_point_plot, = self.PhaseDiagram.phase_diagram_plot_drawing.plot([0.0], [0.0], color=self.pressed_point_desc['color'], marker=self.pressed_point_desc['marker'])
		
		# Plot the phase stability diagram of the ternary compound using the new settings
		self.PhaseDiagram.Plot_PhaseDiagram()
		
		
		if self.show_defects_diagram:
			
			#self.defectsdiagram_window.dopant_selection_box.setEnabled(True)
			self.dopant_selection_box.setEnabled(True)
			
			self.Update_DefectsDiagram_Plot_Elements_Function()	# Found in tab_quaternary.../tab_ternary... script
			
			self.Generate_DefectsDiagram_Plot_Function()
		
		
		
		if self.show_carrier_concentration:
			
			#self.defectsdiagram_window.defects_synthesis_temperature_box.setEnabled(True)
			self.defects_synthesis_temperature_box.setEnabled(True)
			self.Generate_CarrierConcentration_Plot_Function()
	
	
	
	
	
	
	###############################################################################################
	################################ Clicking on Phase Diagram ####################################
	###############################################################################################
	
	def Pressed_Point(self, event):
		
		# This function reads in the coordinates of the point in the phase stability diagram plot that the user clicks.

		point_x = event.xdata	# x-coordinate
		point_y = event.ydata	# y-coordinate
		
		if (not isinstance(point_x, float)) or (not isinstance(point_y, float)):	# Check that the user presses somewhere on the plot (and not anywhere else)
			return ""
		elif self.PhaseDiagram.main_compound_plot == None:
			return
		
		
		# Update the clicked mu values
		self.deltamu_values[self.first_element] = point_x
		self.deltamu_values[self.second_element] = point_y
		if self.type == "ternary":
			#self.deltamu_values[self.third_element] = (self.main_compound_enthalpy - self.main_compound_number_first_specie*self.deltamu_values[self.first_element] - self.main_compound_number_second_specie*self.deltamu_values[self.second_element]) / self.main_compound_number_third_specie
			self.deltamu_values[self.third_element] = (	self.main_compound_enthalpy \
														- self.main_compound_info["dft_"+self.first_element]*self.deltamu_values[self.first_element] \
														- self.main_compound_info["dft_"+self.second_element]*self.deltamu_values[self.second_element] \
														) / self.main_compound_info["dft_"+self.third_element]
		elif self.type == "quaternary":
			#self.deltamu_values[self.third_element] = (self.main_compound_enthalpy - self.main_compound_number_first_specie*self.deltamu_values[self.first_element] - self.main_compound_number_second_specie*self.deltamu_values[self.second_element] - self.main_compound_number_fourth_specie*self.deltamu_values[self.fourth_element]) / self.main_compound_number_third_specie
			self.deltamu_values[self.third_element] = (	self.main_compound_enthalpy \
														- self.main_compound_info["dft_"+self.first_element]*self.deltamu_values[self.first_element] \
														- self.main_compound_info["dft_"+self.second_element]*self.deltamu_values[self.second_element] \
														- self.main_compound_info["dft_"+self.fourth_element]*self.deltamu_values[self.fourth_element] \
														) / self.main_compound_info["dft_"+self.third_element]
		
		
		# Update the display ports for the mu values
		self.mu1_display.setText("{0:.4f}".format(self.deltamu_values[self.first_element]))
		self.mu2_display.setText("{0:.4f}".format(self.deltamu_values[self.second_element]))
		self.mu3_display.setText("{0:.4f}".format(self.deltamu_values[self.third_element]))
		
		
		# Update the red dot where the user clicked on the phase diagram
		self.pressed_point_plot.set_data([self.deltamu_values[self.first_element]], [self.deltamu_values[self.second_element]])
		self.PhaseDiagram.phase_diagram_plot_canvas.draw()
		
		
		
		if self.show_defects_diagram:
			
			# Update chemical potentials in defects diagram
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
			self.CarrierConcentration.mu_elements[self.first_element]["deltamu"] = self.deltamu_values[self.first_element]
			self.CarrierConcentration.mu_elements[self.second_element]["deltamu"] = self.deltamu_values[self.second_element]
			self.CarrierConcentration.mu_elements[self.third_element]["deltamu"] = self.deltamu_values[self.third_element]
			
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






