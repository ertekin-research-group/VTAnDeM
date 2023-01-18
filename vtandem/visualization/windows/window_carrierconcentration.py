
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'
__author__ = 'Michael_Lidia_Jiaxing_Elif'


import numpy as np
import os

import PyQt5
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

script_path = os.path.dirname(__file__)
vtandem_source_path = "/".join(script_path.split("/")[:-3])

title_font = 16

class Window_CarrierConcentration(QWidget):
	
	def __init__(self):
		
		QWidget.__init__(self)
		
		###### Main carrier concentration window widget
		self.carrierconcentration_window = QWidget()										# One of the main sub-widgets is where the user defines the settings of the plots.
		self.carrierconcentration_window_layout = QVBoxLayout(self.carrierconcentration_window)	# The settings should be placed on top of one another, i.e. vertically.
		
		# Title
		self.carrierconcentration_title = QLabel("Carrier \n Concentration")				# QLabel is a widget that displays text
		self.carrierconcentration_title.setAlignment(Qt.AlignCenter)					# Align the text to center
		self.carrierconcentration_title_font = QFont("sans-serif", title_font, QFont.Bold) 	# Declare font
		self.carrierconcentration_title.setFont(self.carrierconcentration_title_font)		# Set the font for the QLabel text
		self.carrierconcentration_window_layout.addWidget(self.carrierconcentration_title)

		# Carrier concentration plot
		self.carrierconcentration_plot = self.CarrierConcentration.carrier_concentration_plot_canvas
		self.carrierconcentration_window_layout.addWidget(self.carrierconcentration_plot)
		
		# Add y-axis limits boxes
		self.Activate_CarrierConcentration_YLimits()

		# Add temperature and equilibrium Fermi energy
		self.Activate_Equilibrium_Fermi_Energy_Settings()

		# Add effective mass boxes
		self.Activate_EffectiveMass_Settings()

		# (WIDGET) Save carrier concentration plot as figure
		self.carrier_concentration_savefigure_button = QPushButton("Save Carrier Concentration Plot")
		self.carrier_concentration_savefigure_button.clicked[bool].connect(lambda: self.CarrierConcentration.SaveFigure())
		self.carrierconcentration_window_layout.addWidget(self.carrier_concentration_savefigure_button)
	
	
	
	###############################################################################################
	############################# Generate Carrier Concentration Plot #############################
	###############################################################################################
	
	def Generate_CarrierConcentration_Plot_Function(self):
		
		self.Update_CarrierConcentration_Plot_Elements_Function()	# Found in tab_quaternary.../tab_ternary... script
		
		# Reset carrier concentration plot
		self.CarrierConcentration.carrier_concentration_plot_drawing.remove()
		self.CarrierConcentration.carrier_concentration_plot_drawing = self.CarrierConcentration.carrier_concentration_plot_figure.add_subplot(111)
		self.CarrierConcentration.Activate_CarrierConcentration_Plot_Axes()
		
		# Plot the carrier concentration (holes and electrons)
		self.CarrierConcentration.carrier_concentration_intrinsic_defect_hole_plot = None
		self.CarrierConcentration.carrier_concentration_intrinsic_defect_electron_plot = None
		self.CarrierConcentration.carrier_concentration_total_hole_plot = None
		self.CarrierConcentration.carrier_concentration_total_electron_plot = None
		
		self.CarrierConcentration.Initialize_CarrierConcentration_Plot()

		# Plot the equilibrium Fermi energy
		if self.DefectsDiagram.intrinsic_defect_plots != {}:
			self.Update_Equilibrium_Fermi_Energy_Temperature()
	


	###############################################################################################
	################################ Carrier Concentration Limits #################################
	###############################################################################################
	
	def Activate_CarrierConcentration_YLimits(self):

		# (WIDGET) Y-axis limits for carrier concentration
		self.carrierconcentration_viewport = QWidget()
		self.carrierconcentration_viewport_layout = QHBoxLayout(self.carrierconcentration_viewport)
		
		# Y min label
		self.carrierconcentration_Ymin_label = QLabel(u"y"+"<sub>min</sub>")
		self.carrierconcentration_Ymin_label.setAlignment(Qt.AlignCenter)
		self.carrierconcentration_viewport_layout.addWidget(self.carrierconcentration_Ymin_label)
		
		# Y min value
		self.carrierconcentration_Ymin_box = QLineEdit("1E16")
		self.carrierconcentration_Ymin_box.editingFinished.connect(lambda: self.CarrierConcentration.Update_WindowSize("YMin", self.carrierconcentration_Ymin_box))
		self.carrierconcentration_viewport_layout.addWidget(self.carrierconcentration_Ymin_box)

		# Y max label
		self.carrierconcentration_Ymax_label = QLabel(u"y"+"<sub>max</sub>")
		self.carrierconcentration_Ymax_label.setAlignment(Qt.AlignCenter)
		self.carrierconcentration_viewport_layout.addWidget(self.carrierconcentration_Ymax_label)
		
		# Y max value
		self.carrierconcentration_Ymax_box = QLineEdit("1E23")
		self.carrierconcentration_Ymax_box.editingFinished.connect(lambda: self.CarrierConcentration.Update_WindowSize("YMax", self.carrierconcentration_Ymax_box))
		self.carrierconcentration_viewport_layout.addWidget(self.carrierconcentration_Ymax_box)
		
		# Add widget to window
		self.carrierconcentration_window_layout.addWidget(self.carrierconcentration_viewport)

	

	###############################################################################################
	########################## Temperature and Equilibrium Fermi Energy ###########################
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
		self.temperature_selection_box.setEnabled(False)
		self.equilibrium_fermi_energy_widget_layout.addWidget(self.temperature_selection_box)
		
		# Spacer item
		self.equilibrium_fermi_energy_widget_layout.addItem(QSpacerItem(50, 20, QSizePolicy.Expanding, QSizePolicy.Minimum))
		
		# Label for equilibrium Fermi energy
		equilibrium_fermi_energy_label = QLabel(u"E"+"<sub>f</sub>"+"<sup>eq</sup> (eV) = ")
		equilibrium_fermi_energy_label.setAlignment(Qt.AlignCenter)
		self.equilibrium_fermi_energy_widget_layout.addWidget(equilibrium_fermi_energy_label)
		
		# Display for equilibrium Fermi energy
		self.equilibrium_fermi_energy_display = QLineEdit("0.00000")
		self.equilibrium_fermi_energy_display.setMaxLength(7)
		self.equilibrium_fermi_energy_display.setEnabled(False)
		self.equilibrium_fermi_energy_display.setStyleSheet("""QLineEdit { background-color: white; color: black }""")
		self.equilibrium_fermi_energy_widget_layout.addWidget(self.equilibrium_fermi_energy_display)
		
		# Spacer item
		self.equilibrium_fermi_energy_widget_layout.addItem(QSpacerItem(50, 20, QSizePolicy.Expanding, QSizePolicy.Minimum))
		
		# Add widget to window
		self.carrierconcentration_window_layout.addWidget(self.equilibrium_fermi_energy_widget)
	
	
	
	def Update_Equilibrium_Fermi_Energy_Temperature(self):
		
		temperature = float(self.temperature_selection_box.currentText())
		intrinsic_equilibrium_fermi_energy = self.CarrierConcentration.intrinsic_equilibrium_fermi_energy[temperature]
		total_equilibrium_fermi_energy = self.CarrierConcentration.total_equilibrium_fermi_energy[temperature]
		
		self.equilibrium_fermi_energy_display.setText(str(total_equilibrium_fermi_energy))
		
		self.DefectsDiagram.Plot_Equilibrium_Fermi_Energy(equilibrium_fermi_energy=total_equilibrium_fermi_energy)

		if (self.defectconc_stat_box.currentText() == "Maxwell-Boltzmann") and self.DefectsDiagram.negative_formation_energy:
			QMessageBox.about(self, "WARNING", "Negative defect formation energy detected.\nChanging to Fermi-Dirac statistics.")
			self.defectconc_stat_box.setCurrentIndex(1) # Switch to Fermi-Dirac
			self.CarrierConcentration.Update_DefectConcentration_DistributionStatistics(self.defectconc_stat_box)
			self.Update_Equilibrium_Fermi_Energy_Temperature()



	###############################################################################################
	############################ Effective Mass, Parabolic Band Model #############################
	###############################################################################################
	
	def Activate_EffectiveMass_Settings(self):

		# (WIDGET) Effective masses of DOS for carrier concentrations
		self.effective_masses_widget = QWidget()
		self.effective_masses_widget_layout = QVBoxLayout(self.effective_masses_widget)
		
		# Label for "Calculated from:"
		dos_calculated_label = QLabel("Carrier concentrations calculated from:")
		dos_calculated_label.setAlignment(Qt.AlignCenter)
		self.effective_masses_widget_layout.addWidget(dos_calculated_label)

		# Widgets for inputs
		self.effective_mass_inputs_widget = QWidget()
		self.effective_mass_inputs_widget_layout = QHBoxLayout(self.effective_mass_inputs_widget)

		# Options bar for DFT vs. parabolic band
		self.dos_option_box = QComboBox()
		if self.CarrierConcentration.dos_data["DOS"] != {}:
			self.dos_option_box.addItem("Real DOS")
		self.dos_option_box.addItem("Parabolic Approx.")
		self.dos_option_box.setCurrentIndex(0)
		self.dos_option_box.setEnabled(False)
		self.effective_mass_inputs_widget_layout.addWidget(self.dos_option_box)

		# Help button
		self.dos_type_question_button = QPushButton()
		self.dos_type_question_button.setIcon(QIcon(vtandem_source_path+"/icon/QuestionIcon.png"))
		self.dos_type_question_button.clicked[bool].connect(self.DOS_Type_Help_Function)
		self.effective_mass_inputs_widget_layout.addWidget(self.dos_type_question_button)

		# Label for hole effective mass
		hole_effective_mass_label = QLabel(u"m"+"<sub>h</sub>"+"<sup>*</sup> (m"+"<sub>0</sub>"+") = ")
		hole_effective_mass_label.setAlignment(Qt.AlignCenter)
		self.effective_mass_inputs_widget_layout.addWidget(hole_effective_mass_label)

		# Input for hole effective mass
		self.hole_effective_mass_input = QLineEdit("1.00")
		self.hole_effective_mass_input.setEnabled(False)
		self.effective_mass_inputs_widget_layout.addWidget(self.hole_effective_mass_input)

		# Label for electron effective mass
		electron_effective_mass_label = QLabel(u"m"+"<sub>e</sub>"+"<sup>*</sup> (m"+"<sub>0</sub>"+") = ")
		electron_effective_mass_label.setAlignment(Qt.AlignCenter)
		self.effective_mass_inputs_widget_layout.addWidget(electron_effective_mass_label)

		# Input for electron effective mass
		self.electron_effective_mass_input = QLineEdit("1.00")
		self.electron_effective_mass_input.setEnabled(False)
		self.effective_mass_inputs_widget_layout.addWidget(self.electron_effective_mass_input)

		# Add functionalities to effective mass inputs
		self.dos_option_box.activated.connect(lambda: self.CarrierConcentration.Update_DOS_Data(self.dos_option_box, self.hole_effective_mass_input, self.electron_effective_mass_input))
		self.hole_effective_mass_input.editingFinished.connect(lambda: self.CarrierConcentration.Update_EffMass(self.hole_effective_mass_input, self.electron_effective_mass_input))
		self.hole_effective_mass_input.editingFinished.connect(self.Update_Equilibrium_Fermi_Energy_Temperature)
		self.electron_effective_mass_input.editingFinished.connect(lambda: self.CarrierConcentration.Update_EffMass(self.hole_effective_mass_input, self.electron_effective_mass_input))
		self.electron_effective_mass_input.editingFinished.connect(self.Update_Equilibrium_Fermi_Energy_Temperature)

		# Add inputs widget to overarching widget
		self.effective_masses_widget_layout.addWidget(self.effective_mass_inputs_widget)

		# Add widget to window
		self.carrierconcentration_window_layout.addWidget(self.effective_masses_widget)

		# If DOS hasn't been imported, update and calculate DOS using parabolic band approximation
		if self.CarrierConcentration.dos_data["DOS"] == {}:
			self.CarrierConcentration.energies_ValenceBand = np.linspace(-20, self.CarrierConcentration.EVBM, 2000)
			self.CarrierConcentration.energies_ConductionBand = np.linspace(self.CarrierConcentration.ECBM, 20, 2000)
			self.hole_effective_mass_input.setEnabled(True)
			self.electron_effective_mass_input.setEnabled(True)
			self.CarrierConcentration.Update_CarrierConcentrations_EffMasses(self.hole_effective_mass_input, self.electron_effective_mass_input)



	def DOS_Type_Help_Function(self):
		dialog_instructions = 	"""
The free electron and hole concentrations can be calculated using either:
	- Real DOS (i.e. calculated using DFT), or
	- Parabolic Approximation (i.e. using the DOS of a single parabolic band).\n
For more info, see:
M.Y. Toriyama, et al., "How to analyse a density of states," Mater. Today. Elec. 1, 100002 (2022).
M.Y. Toriyama, et al., "VTAnDeM: A Python Toolkit for Simultaneously Visualizing Phase Stability, Defect Energetics, and Carrier Concentrations of Materials," Submitted.
								"""

		self.dos_type_message_window = QMainWindow()
		self.dos_type_message_window.setWindowTitle("Help")
		self.dos_type_message_window.setWindowIcon(QIcon(vtandem_source_path+"/logo/LogoSmall.png"))
		self.dos_type_message_widget = QWidget()
		self.dos_type_message_widget_layout = QVBoxLayout(self.dos_type_message_widget)
		dos_type_dialog_instructions = QLabel(dialog_instructions)
		dos_type_dialog_instructions.setTextInteractionFlags(Qt.TextSelectableByMouse)
		self.dos_type_message_widget_layout.addWidget(dos_type_dialog_instructions)
		self.dos_type_message_window.setCentralWidget(self.dos_type_message_widget)
		self.dos_type_message_window.show()
	





