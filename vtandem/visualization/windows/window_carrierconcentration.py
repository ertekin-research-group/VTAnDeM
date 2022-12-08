
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'
__author__ = 'Michael_Lidia_Jiaxing_Elif'


import numpy as np

import PyQt5
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

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
		
		self.carrierconcentration_viewport = QWidget()
		self.carrierconcentration_viewport_layout = QHBoxLayout(self.carrierconcentration_viewport)
		
		# (WIDGET) Y-axis limits for carrier concentration
		self.carrierconcentration_Ymin_label = QLabel(u"y"+"<sub>min</sub>")
		self.carrierconcentration_Ymin_label.setAlignment(Qt.AlignCenter)
		self.carrierconcentration_viewport_layout.addWidget(self.carrierconcentration_Ymin_label)
		self.carrierconcentration_Ymin_box = QLineEdit("1E16")
		self.carrierconcentration_Ymin_box.editingFinished.connect(lambda: self.CarrierConcentration.Update_WindowSize("YMin", self.carrierconcentration_Ymin_box))
		self.carrierconcentration_viewport_layout.addWidget(self.carrierconcentration_Ymin_box)
		self.carrierconcentration_Ymax_label = QLabel(u"y"+"<sub>max</sub>")
		self.carrierconcentration_Ymax_label.setAlignment(Qt.AlignCenter)
		self.carrierconcentration_viewport_layout.addWidget(self.carrierconcentration_Ymax_label)
		self.carrierconcentration_Ymax_box = QLineEdit("1E23")
		self.carrierconcentration_Ymax_box.editingFinished.connect(lambda: self.CarrierConcentration.Update_WindowSize("YMax", self.carrierconcentration_Ymax_box))
		self.carrierconcentration_viewport_layout.addWidget(self.carrierconcentration_Ymax_box)
				
		self.carrierconcentration_window_layout.addWidget(self.carrierconcentration_viewport)
		
		self.Activate_Equilibrium_Fermi_Energy_Settings()
		
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
		self.equilibrium_fermi_energy_display = QLineEdit("0.00000")
		self.equilibrium_fermi_energy_display.setMaxLength(7)
		self.equilibrium_fermi_energy_display.setEnabled(False)
		self.equilibrium_fermi_energy_display.setStyleSheet("""QLineEdit { background-color: white; color: black }""")
		self.equilibrium_fermi_energy_widget_layout.addWidget(self.equilibrium_fermi_energy_display)
		
		self.equilibrium_fermi_energy_widget_layout.addItem(QSpacerItem(50, 20, QSizePolicy.Expanding, QSizePolicy.Minimum))
		
		self.carrierconcentration_window_layout.addWidget(self.equilibrium_fermi_energy_widget)
	
	
	def Update_Equilibrium_Fermi_Energy_Temperature(self):
		
		temperature = float(self.temperature_selection_box.currentText())
		intrinsic_equilibrium_fermi_energy = self.CarrierConcentration.intrinsic_equilibrium_fermi_energy[temperature]
		total_equilibrium_fermi_energy = self.CarrierConcentration.total_equilibrium_fermi_energy[temperature]
		
		self.equilibrium_fermi_energy_display.setText(str(total_equilibrium_fermi_energy))
		
		if (str(total_equilibrium_fermi_energy) == "< EVBM") or (str(total_equilibrium_fermi_energy) == "> ECBM"):
			self.equilibrium_fermi_energy_display.setStyleSheet("""QLineEdit { background-color: white; color: red }""")
		else:
			self.equilibrium_fermi_energy_display.setStyleSheet("""QLineEdit { background-color: white; color: black }""")

		self.DefectsDiagram.Plot_Equilibrium_Fermi_Energy(equilibrium_fermi_energy=total_equilibrium_fermi_energy)

