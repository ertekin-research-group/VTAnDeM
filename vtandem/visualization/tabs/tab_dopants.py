
__author__ = 'Michael_Lidia_Jiaxing_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.path as mpltPath
from copy import deepcopy

import PyQt5
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

# Import compositional phase diagram objects
from vtandem.visualization.ternary.ternary_plots.plot_composition_ternary_phase_diagram import Plot_Composition_Ternary_PhaseDiagram
from vtandem.visualization.quaternary.quaternary_plots.plot_composition_quaternary_phase_diagram import Plot_Composition_Quaternary_PhaseDiagram

# Import defect diagram objects
from vtandem.visualization.binary.binary_plots.plot_binary_defects_diagram import Plot_Binary_DefectsDiagram
from vtandem.visualization.ternary.ternary_plots.plot_ternary_defects_diagram import Plot_Ternary_DefectsDiagram

# Import carrier concentration objects
from vtandem.visualization.binary.binary_plots.plot_binary_carrier_concentration import Plot_Binary_Carrier_Concentration
from vtandem.visualization.ternary.ternary_plots.plot_ternary_carrier_concentration import Plot_Ternary_Carrier_Concentration

from vtandem.visualization.windows.window_phasediagram_composition import Window_Compositional_PhaseDiagram
from vtandem.visualization.windows.window_defectsdiagram import Window_DefectsDiagram
from vtandem.visualization.windows.window_carrierconcentration import Window_CarrierConcentration

from vtandem.dft.obtain_dft import *


class Tab_Dopants(Window_DefectsDiagram, Window_CarrierConcentration):
	
	def __init__(self, 	type: str, \
						main_compound = None, \
						compounds_info = None, \
						defects_data = None, \
						main_compound_info = None, \
						dos_data = None, \
						show_defects_diagram = True, \
						show_carrier_concentration = True
						):	# User specifies the main compound and its constituents
		

		# Check that the 'type' argument is legitimate
		if (type != "binary") and (type != "ternary"):
			raise Exception("The argument 'type' must be either 'binary' or 'ternary'. Exiting...")
		self.type = type


		###############################################################################################
		########################### Initialize materials-related variables ############################
		###############################################################################################
		
		# Initialize the main ternary compound
		self.main_compound = main_compound
		
		# Get extrinsic defects
		self.extrinsic_defects = []
		for defect in defects_data.keys():
			if "_" not in defect:
				continue
			if defects_data[defect]["Extrinsic"] == "Yes":
				self.extrinsic_defects.append(defect)
		
		# Track selected dopant
		try:
			self.dopant = self.extrinsic_defects[0].split("_")[0]
		except:
			self.dopant = None
		
		# Obtain DFT compounds data
		elements_list_withdopant = deepcopy(self.elements_list)
		elements_list_withdopant.append(self.dopant)
		self.compounds_info = Obtain_Compounds_Data( elements_list_withdopant )
		self.main_compound_info = main_compound_info
		
		# Plot settings
		self.show_defects_diagram = show_defects_diagram
		self.show_carrier_concentration = show_carrier_concentration


		###############################################################################################
		################################# Compositional phase diagram #################################
		###############################################################################################
		
		if self.type == "binary":
			self.Compositional_PhaseDiagram = Plot_Composition_Ternary_PhaseDiagram(main_compound = self.main_compound, first_element = self.first_element, second_element = self.second_element, third_element = self.dopant, compounds_info = self.compounds_info, main_compound_info = self.main_compound_info)
		elif self.type == "ternary":
			self.Compositional_PhaseDiagram = Plot_Composition_Quaternary_PhaseDiagram(main_compound = self.main_compound, first_element = self.first_element, second_element = self.second_element, third_element = self.third_element, fourth_element = self.dopant, compounds_info = self.compounds_info, main_compound_info = self.main_compound_info)
		self.Compositional_PhaseDiagram.composition_phasediagram_plot_figure.canvas.mpl_connect('button_press_event', self.Generate_DefectsDiagram_Plot_Function)
		self.Compositional_PhaseDiagram.composition_phasediagram_plot_figure.canvas.mpl_connect('button_press_event', self.Generate_CarrierConcentration_Plot_Function)
		

		if self.show_defects_diagram:

			# Defects diagram
			if self.type == "binary":
				self.DefectsDiagram = Plot_Binary_DefectsDiagram(main_compound = self.main_compound, first_element = self.first_element, second_element = self.second_element)
			elif self.type == "ternary":
				self.DefectsDiagram = Plot_Ternary_DefectsDiagram(main_compound = self.main_compound, first_element = self.first_element, second_element = self.second_element, third_element = self.third_element)
			self.DefectsDiagram.defects_data = defects_data
			self.DefectsDiagram.main_compound_info = main_compound_info

			for element in self.elements_list:
				self.DefectsDiagram.mu_elements[element]["mu0"] = self.compounds_info[element]["mu0"]
				self.DefectsDiagram.mu_elements[element]["deltamu"] = self.Compositional_PhaseDiagram.deltamu_values[element]

			self.DefectsDiagram.EVBM = main_compound_info["VBM"]
			self.DefectsDiagram.ECBM = self.DefectsDiagram.EVBM + main_compound_info["BandGap"]
			self.DefectsDiagram.axis_lims["XMin"] = 0.0
			self.DefectsDiagram.axis_lims["XMax"] = main_compound_info["BandGap"]
			self.DefectsDiagram.fermi_energy_array = np.linspace(self.DefectsDiagram.EVBM, self.DefectsDiagram.ECBM, 2000)
			self.DefectsDiagram.Activate_DefectsDiagram_Plot_Axes()

			self.DefectsDiagram.dopant = self.dopant
			self.DefectsDiagram.dopant_mu0 = self.compounds_info[self.dopant]["mu0"]
			self.DefectsDiagram.extrinsic_defects = self.extrinsic_defects
		


		if self.show_carrier_concentration:

			# Carrier concentration
			if self.type == "binary":
				self.CarrierConcentration = Plot_Binary_Carrier_Concentration(main_compound = self.main_compound, first_element = self.first_element, second_element = self.second_element)
			elif self.type == "ternary":
				self.CarrierConcentration = Plot_Ternary_Carrier_Concentration(main_compound = self.main_compound, first_element = self.first_element, second_element = self.second_element, third_element = self.third_element)
			self.CarrierConcentration.defects_data = defects_data
			self.CarrierConcentration.main_compound_info = main_compound_info

			try:
				self.CarrierConcentration.dos_data = dos_data[self.main_compound]
			except:
				self.CarrierConcentration.dos_data = {"DOS": {}, "Volume": 0.0}

			for element in self.elements_list:
				self.CarrierConcentration.mu_elements[element]["mu0"] = self.compounds_info[element]["mu0"]
				self.CarrierConcentration.mu_elements[element]["deltamu"] = self.Compositional_PhaseDiagram.deltamu_values[element]

			self.CarrierConcentration.vol = main_compound_info["Volume"]
			self.CarrierConcentration.EVBM = self.DefectsDiagram.EVBM
			self.CarrierConcentration.ECBM = self.DefectsDiagram.ECBM
			# Fermi energies should sample outside band gap, in case EFeq is not in gap
			self.CarrierConcentration.fermi_energy_array = np.linspace(self.DefectsDiagram.EVBM-1.0, self.DefectsDiagram.ECBM+1.0, 2000)
			self.CarrierConcentration.Activate_CarrierConcentration_Plot_Axes()
			if self.CarrierConcentration.dos_data["DOS"] != {}:
				self.CarrierConcentration.Organize_DOS_Data()
				self.CarrierConcentration.Extract_Relevant_Energies_DOSs()
				self.CarrierConcentration.Calculate_Hole_Electron_Concentration_Matrices()
			
			
		
		
		###############################################################################################
		###############################################################################################
		#################################### Initialize fourth tab ####################################
		###############################################################################################
		###############################################################################################
		
		self.tab4 = QWidget()
		self.tab4_layout = QHBoxLayout(self.tab4)
		
		# Add compositional phase diagram window widget to tab4
		self.Compositional_PhaseDiagram_Window = Window_Compositional_PhaseDiagram(self.main_compound, self.Compositional_PhaseDiagram)
		self.Compositional_PhaseDiagram_Window.compositional_phasediagram_window_layout.addWidget(self.Compositional_PhaseDiagram_Window.composition_phase_diagram_title)
		self.Compositional_PhaseDiagram_Window.compositional_phasediagram_window_layout.addWidget(self.Compositional_PhaseDiagram_Window.composition_phase_diagram_plot)
		
		# Add dopant selection box
		self.dopant_selection_widget = QWidget()
		self.dopant_selection_widget_layout = QHBoxLayout(self.dopant_selection_widget)
		self.dopant_selection_label = QLabel("Dopant: ")
		self.dopant_selection_label.setAlignment(Qt.AlignCenter)
		self.dopant_selection_widget_layout.addWidget(self.dopant_selection_label)
		self.dopant_selection_box = QComboBox()
		for defect in self.extrinsic_defects:
			dopant = defect.split("_")[0]
			if dopant not in [self.dopant_selection_box.itemText(i) for i in range(self.dopant_selection_box.count())]:
				self.dopant_selection_box.addItem(dopant)
		self.dopant_selection_box.setCurrentIndex(0)
		self.dopant_selection_box.activated.connect(self.Update_Dopant)
		self.dopant_selection_widget_layout.addWidget(self.dopant_selection_box)
		self.Compositional_PhaseDiagram_Window.compositional_phasediagram_window_layout.addWidget(self.dopant_selection_widget)
		
		# Add save figures button
		self.Compositional_PhaseDiagram_Window.compositional_phasediagram_window_layout.addWidget(self.Compositional_PhaseDiagram_Window.phasediagram_savefigure_button)
		
		# Add compositional phase diagram window object to tab4
		self.tab4_layout.addWidget(self.Compositional_PhaseDiagram_Window.compositional_phasediagram_window)
		
		# Add defect formation energy diagram window widget to tab4
		if self.show_defects_diagram:
			Window_DefectsDiagram.__init__(self, show_dopant=False)
			self.tab4_layout.addWidget(self.defectsdiagram_window)

		# Add carrier concentration window widget to tab4
		if self.show_carrier_concentration:
			Window_CarrierConcentration.__init__(self)
			self.tab4_layout.addWidget(self.carrierconcentration_window)



	###############################################################################################
	################################ Update Dopant in Phase Diagram ###############################
	###############################################################################################
	
	def Update_Dopant(self, event):
		
		self.dopant = self.dopant_selection_box.currentText()

		self.Compositional_PhaseDiagram.dopant = self.dopant
		if self.type == "binary":
			self.Compositional_PhaseDiagram.elements_list = [self.first_element, self.second_element, self.dopant]
		elif self.type == "ternary":
			self.Compositional_PhaseDiagram.elements_list = [self.first_element, self.second_element, self.third_element, self.dopant]

		self.compounds_info = Obtain_Compounds_Data( self.Compositional_PhaseDiagram.elements_list )
		self.Compositional_PhaseDiagram.compounds_info = self.compounds_info

		self.Compositional_PhaseDiagram.Generate_Compositional_PhaseDiagram(self.Compositional_PhaseDiagram.compounds_info, self.Compositional_PhaseDiagram.elements_list)

		self.Compositional_PhaseDiagram.Find_All_PhaseRegions()
		self.Compositional_PhaseDiagram.Plot_Centroids()
		
		self.DefectsDiagram.dopant = self.dopant
		self.DefectsDiagram.dopant_mu0 = self.compounds_info[self.dopant]["mu0"]

		self.CarrierConcentration.dopant = self.dopant
		self.CarrierConcentration.dopant_mu0 = self.compounds_info[self.dopant]["mu0"]
	
	
	
	###############################################################################################
	################################# Generate Defects Diagram ####################################
	###############################################################################################
	
	def Generate_DefectsDiagram_Plot_Function(self, event):
		
		# Check whether a region is actually clicked
		if not self.Check_Click_Function(event):
			return

		# Update elements and chemical potentials
		self.DefectsDiagram.dopant_deltamu = self.Compositional_PhaseDiagram.deltamu_values[self.dopant]
		self.DefectsDiagram.Update_Deltamus(self.Compositional_PhaseDiagram.deltamu_values)

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
		self.DefectsDiagram.Initialize_Extrinsic_DefectsDiagram_Plot()



	###############################################################################################
	################################ Generate Carrier Concentration ###############################
	###############################################################################################
	
	def Generate_CarrierConcentration_Plot_Function(self, event):
		
		# Check whether a region is actually clicked
		if not self.Check_Click_Function(event):
			return

		# Update elements and chemical potentials
		self.CarrierConcentration.dopant_deltamu = self.Compositional_PhaseDiagram.deltamu_values[self.dopant]
		self.CarrierConcentration.Update_Deltamus(self.Compositional_PhaseDiagram.deltamu_values)

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
	################################ Generate Carrier Concentration ###############################
	###############################################################################################
	
	def Check_Click_Function(self, event):
		
		# This function checks whether a region in the compositional phase diagram was
		# 	actually clicked by the user. If not, then nothing will happen (as shown by 
		# 	the usage of this function in Generate_DefectsDiagram_Plot_Function and 
		# 	Generate_CarrierConcentration_Plot_Function).

		if self.type == "binary":
			point_x = event.xdata
			point_y = event.ydata
			is_in_any_threephaseregion = False
			for three_phase_region in self.Compositional_PhaseDiagram.phase_region_objects:
				is_in_threephaseregion = mpltPath.Path(three_phase_region.vertices).contains_point([point_x, point_y])
				if is_in_threephaseregion:
					is_in_any_threephaseregion = True
			if not is_in_any_threephaseregion:
				return False
		elif self.type == "ternary":
			# Check that a centroid of one of the four-phase regions has been clicked on
			contains, index = self.Compositional_PhaseDiagram.centroids_plot.contains(event)
			if not contains:
				return False

		return True


