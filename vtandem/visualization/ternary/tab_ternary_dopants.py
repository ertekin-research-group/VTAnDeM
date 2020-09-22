
__author__ = 'Michael_Lidia_Jiaxing_Benita_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'


import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

import PyQt5
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

from vtandem.visualization.quaternary.quaternary_scripts.plot_composition_quaternary_phase_diagram import Composition_Quaternary_PhaseDiagram
from vtandem.visualization.ternary.ternary_scripts.plot_ternary_defects_diagram import Ternary_Defects_Diagram

from vtandem.visualization.widgets.phasediagram_window_composition import Compositional_PhaseDiagram_Window



class Tab_Ternary_Dopants(object):
	
	def __init__(self, main_compound=None, first_element=None, second_element=None, third_element=None, compounds_info=None, defects_data=None):	# User specifies the main compound and its constituents
		
		###############################################################################################
		########################### Initialize materials-related variables ############################
		###############################################################################################
		
		# Initialize the main ternary compound
		self.main_compound = main_compound
		
		# Label the first, second, and third species of the atoms in the ternary compound
		self.first_element = first_element
		self.second_element = second_element
		self.third_element = third_element
		self.elements_list = [self.first_element, self.second_element, self.third_element]					# Species list (order MAY change)
		
		# Obtain DFT data
		self.compounds_info = compounds_info
		self.defects_data = defects_data
		
		# Get extrinsic defects
		self.extrinsic_defects = []
		for defect in self.defects_data[self.main_compound].keys():
			if "_" not in defect:
				continue
			if self.defects_data[self.main_compound][defect]["Extrinsic"] == "Yes":
				self.extrinsic_defects.append(defect)
		
		# Track selected dopant
		self.dopant = self.extrinsic_defects[0].split("_")[0]
		
		
		###############################################################################################
		################################# Compositional phase diagram #################################
		###############################################################################################
		
		self.Compositional_PhaseDiagram = Composition_Quaternary_PhaseDiagram(main_compound = self.main_compound, first_element = self.dopant, second_element = self.first_element, third_element = self.second_element, fourth_element = self.third_element, compounds_info = self.compounds_info)
		self.Compositional_PhaseDiagram.composition_phasediagram_plot_figure.canvas.mpl_connect('button_press_event', self.Generate_DefectsDiagram_Plot_Function)
		
		
		
		# Defects diagram
		self.DefectsDiagram = Ternary_Defects_Diagram(main_compound = self.main_compound, first_element = self.first_element, second_element = self.second_element, third_element = self.third_element)
		self.DefectsDiagram.defects_data = self.defects_data[self.main_compound]
		self.DefectsDiagram.main_compound_number_first_specie = self.compounds_info[main_compound][self.first_element]
		self.DefectsDiagram.main_compound_number_second_specie = self.compounds_info[main_compound][self.second_element]
		self.DefectsDiagram.main_compound_number_third_specie = self.compounds_info[main_compound][self.third_element]
		self.DefectsDiagram.main_compound_total_energy = self.compounds_info[main_compound]["total_energy"]
		self.DefectsDiagram.mu_elements[self.first_element]["mu0"] = self.compounds_info[self.first_element]["mu0"]
		self.DefectsDiagram.mu_elements[self.second_element]["mu0"] = self.compounds_info[self.second_element]["mu0"]
		self.DefectsDiagram.mu_elements[self.third_element]["mu0"] = self.compounds_info[self.third_element]["mu0"]
		self.DefectsDiagram.mu_elements[self.first_element]["deltamu"] = self.Compositional_PhaseDiagram.deltamu_values[self.first_element]
		self.DefectsDiagram.mu_elements[self.second_element]["deltamu"] = self.Compositional_PhaseDiagram.deltamu_values[self.second_element]
		self.DefectsDiagram.mu_elements[self.third_element]["deltamu"] = self.Compositional_PhaseDiagram.deltamu_values[self.third_element]
		self.DefectsDiagram.EVBM = self.compounds_info[main_compound]["vbm"]
		self.DefectsDiagram.ECBM = self.DefectsDiagram.EVBM + self.compounds_info[main_compound]["band_gap"]
		self.DefectsDiagram.fermi_energy_array = np.linspace(self.DefectsDiagram.EVBM, self.DefectsDiagram.ECBM, 100)
		self.DefectsDiagram.Activate_DefectsDiagram_Plot_Axes()
		
		
		
		###############################################################################################
		###############################################################################################
		#################################### Initialize fourth tab ####################################
		###############################################################################################
		###############################################################################################
		
		self.tab4 = QWidget()
		self.tab4_layout = QHBoxLayout(self.tab4)
		
		### Compositional phase diagram window
		# Initialize compositional phase diagram window object
		self.Compositional_PhaseDiagram_Window = Compositional_PhaseDiagram_Window(self.main_compound, self.Compositional_PhaseDiagram)
		# Add compound title to window
		self.Compositional_PhaseDiagram_Window.compositional_phasediagram_window_layout.addWidget(self.Compositional_PhaseDiagram_Window.compound_title)
		# Add compositional phase diagram plot
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
		
		
		
		# Defect formation energy diagram window widget
		self.defectsdiagram_window = QWidget()									# One of the main sub-widgets is where the user defines the settings of the plots.
		self.defectsdiagram_window_layout = QVBoxLayout(self.defectsdiagram_window)		# The settings should be placed on top of one another, i.e. vertically.
		
		# Defects diagram plot
		self.defects_diagram_plot = self.DefectsDiagram.defects_diagram_plot_canvas
		self.defectsdiagram_window_layout.addWidget(self.defects_diagram_plot)
		
		# Y-axis limits for defects diagram
		self.defectsdiagram_viewport = QWidget()
		self.defectsdiagram_viewport_layout = QHBoxLayout(self.defectsdiagram_viewport)
		
		# Y-axis limits for defects diagram
		self.defectsdiagram_Ymin_label = QLabel(u"y"+"<sub>min</sub>")
		self.defectsdiagram_Ymin_label.setAlignment(Qt.AlignRight)
		self.defectsdiagram_viewport_layout.addWidget(self.defectsdiagram_Ymin_label)
		self.defectsdiagram_Ymin_box = QLineEdit("-2.0")
		self.defectsdiagram_Ymin_box.editingFinished.connect(lambda: self.DefectsDiagram.Update_WindowSize("YMin", self.defectsdiagram_Ymin_box))
		self.defectsdiagram_viewport_layout.addWidget(self.defectsdiagram_Ymin_box)
		self.defectsdiagram_Ymax_label = QLabel(u"y"+"<sub>max</sub>")
		self.defectsdiagram_Ymax_label.setAlignment(Qt.AlignRight)
		self.defectsdiagram_viewport_layout.addWidget(self.defectsdiagram_Ymax_label)
		self.defectsdiagram_Ymax_box = QLineEdit("2.0")
		self.defectsdiagram_Ymax_box.editingFinished.connect(lambda: self.DefectsDiagram.Update_WindowSize("YMax", self.defectsdiagram_Ymax_box))
		self.defectsdiagram_viewport_layout.addWidget(self.defectsdiagram_Ymax_box)
		self.defectsdiagram_window_layout.addWidget(self.defectsdiagram_viewport)
		
		# (WIDGET) Save defects diagram as figure
		self.defects_diagram_savefigure_button = QPushButton("Save Defects Diagram Figure")
		self.defects_diagram_savefigure_button.clicked[bool].connect(lambda: self.SaveFigure_Function())
		self.defectsdiagram_window_layout.addWidget(self.defects_diagram_savefigure_button)
		
		self.tab4_layout.addWidget(self.defectsdiagram_window)
	
	
	
	###############################################################################################
	################################ Update Dopant in Phase Diagram ###############################
	###############################################################################################
	
	def Update_Dopant(self, event):
		
		self.dopant = self.dopant_selection_box.currentText()
		
		self.Compositional_PhaseDiagram.first_element = self.dopant
		self.Compositional_PhaseDiagram.elements_list = [self.dopant, self.first_element, self.second_element, self.third_element]
		self.Compositional_PhaseDiagram.composition_phasediagram_legend.remove()
		self.Compositional_PhaseDiagram.composition_phasediagram_plot_drawing.remove()
		self.Compositional_PhaseDiagram.composition_phasediagram_plot_drawing = Axes3D(self.Compositional_PhaseDiagram.composition_phasediagram_plot_figure)
		self.Compositional_PhaseDiagram.Create_Compositional_PhaseDiagram()
		self.Compositional_PhaseDiagram.Plot_Compositional_PhaseDiagram()
		
		self.Compositional_PhaseDiagram.fourphaseregions = []
		self.Compositional_PhaseDiagram.fourphaseregion_names = []
		self.Compositional_PhaseDiagram.fourphaseregion_centroids = []
		self.Compositional_PhaseDiagram.Find_All_FourPhaseRegions()
		self.Compositional_PhaseDiagram.Plot_Centroids()
		
		self.Compositional_PhaseDiagram.fourphaseregion_selected = None
		self.Compositional_PhaseDiagram.fourphaseregion_shade = []
		
		self.Compositional_PhaseDiagram.fourphaseregion_annotation = self.Compositional_PhaseDiagram.composition_phasediagram_plot_drawing.annotate("", xy=(0,0), xytext=(20,20), textcoords="offset points",
																																					bbox = dict(boxstyle="round", fc="w"),
																																					arrowprops = dict(arrowstyle = "->") )
		
		self.DefectsDiagram.extrinsic_defect = self.dopant
		self.DefectsDiagram.extrinsic_defect_mu0 = self.compounds_info[self.dopant]["enthalpy"]
	
	
	
	###############################################################################################
	################################# Generate Defects Diagram ####################################
	###############################################################################################
	
	def Generate_DefectsDiagram_Plot_Function(self, event):
		
		# Check that a centroid of one of the four-phase regions has been clicked on
		contains, index = self.Compositional_PhaseDiagram.centroids_plot.contains(event)
		if not contains:
			return
		
		# Update elements and chemical potentials
		self.DefectsDiagram.extrinsic_defect_deltamu = self.Compositional_PhaseDiagram.deltamu_values[self.dopant]
		self.DefectsDiagram.mu_elements[self.first_element]["deltamu"] = self.Compositional_PhaseDiagram.deltamu_values[self.first_element]
		self.DefectsDiagram.mu_elements[self.second_element]["deltamu"] = self.Compositional_PhaseDiagram.deltamu_values[self.second_element]
		self.DefectsDiagram.mu_elements[self.third_element]["deltamu"] = self.Compositional_PhaseDiagram.deltamu_values[self.third_element]
		
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
	
	
	
	
	
	###############################################################################################
	###################################### Save Figure ############################################
	###############################################################################################
	
	def SaveFigure_Function(self):
		
		options = QFileDialog.Options()
		options |= QFileDialog.DontUseNativeDialog
		filename, extension_type = QFileDialog.getSaveFileName(filter = "Portable Network Graphics (*.png);;" \
																+"Portable Document Format (*.pdf);;" \
																+"Scalable Vector Graphics (*.svg);;" \
																+"Encapsulated PostScript (*.eps)", options=options)
		self.DefectsDiagram.SaveFigure(filename, extension_type)





