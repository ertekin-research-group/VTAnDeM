
__author__ = 'Michael_Lidia_Jiaxing_Benita_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'


import numpy as np

import PyQt5
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

from vtandem.visualization.quaternary.quaternary_plots.plot_composition_quaternary_phase_diagram import Plot_Composition_Quaternary_PhaseDiagram
from vtandem.visualization.quaternary.quaternary_plots.plot_quaternary_defects_diagram import Plot_Quaternary_DefectsDiagram

from vtandem.visualization.windows.window_phasediagram_composition import Window_Compositional_PhaseDiagram
from vtandem.visualization.windows.window_defectsdiagram import Window_DefectsDiagram


class Tab_Quaternary_Compositional_PhaseDiagram3D(object):
	
	def __init__(self, main_compound=None, first_element=None, second_element=None, third_element=None, fourth_element=None, compounds_info=None, defects_data=None, show_defects_diagram=True):	# User specifies the main compound and its constituents
		
		###############################################################################################
		########################### Initialize materials-related variables ############################
		###############################################################################################
		
		# Initialize the main quaternary compound
		self.main_compound = main_compound
		
		# Label the first, second, third, and fourth species of the atoms in the quaternary compound
		self.first_element = first_element
		self.second_element = second_element
		self.third_element = third_element
		self.fourth_element = fourth_element
		self.elements_list = [self.first_element, self.second_element, self.third_element, self.fourth_element]					# Species list (order MAY change)
		
		# Obtain DFT data
		self.compounds_info = compounds_info
		self.defects_data = defects_data
		
		###############################################################################################
		################################# Compositional phase diagram #################################
		###############################################################################################
		
		self.Compositional_PhaseDiagram3D = Plot_Composition_Quaternary_PhaseDiagram(main_compound = self.main_compound, first_element = self.first_element, second_element = self.second_element, third_element = self.third_element, fourth_element = self.fourth_element, compounds_info = self.compounds_info)
		
		# Defects diagram
		if show_defects_diagram:
			self.DefectsDiagram = Plot_Quaternary_DefectsDiagram(main_compound = self.main_compound, first_element = self.first_element, second_element = self.second_element, third_element = self.third_element, fourth_element = self.fourth_element)
			self.DefectsDiagram.defects_data = self.defects_data[self.main_compound]
			self.DefectsDiagram.main_compound_number_first_specie = self.compounds_info[main_compound][self.first_element]
			self.DefectsDiagram.main_compound_number_second_specie = self.compounds_info[main_compound][self.second_element]
			self.DefectsDiagram.main_compound_number_third_specie = self.compounds_info[main_compound][self.third_element]
			self.DefectsDiagram.main_compound_number_fourth_specie = self.compounds_info[main_compound][self.fourth_element]
			self.DefectsDiagram.main_compound_total_energy = self.compounds_info[main_compound]["total_energy"]
			self.DefectsDiagram.mu_elements[self.first_element]["mu0"] = self.compounds_info[self.first_element]["mu0"]
			self.DefectsDiagram.mu_elements[self.second_element]["mu0"] = self.compounds_info[self.second_element]["mu0"]
			self.DefectsDiagram.mu_elements[self.third_element]["mu0"] = self.compounds_info[self.third_element]["mu0"]
			self.DefectsDiagram.mu_elements[self.fourth_element]["mu0"] = self.compounds_info[self.fourth_element]["mu0"]
			self.DefectsDiagram.mu_elements[self.first_element]["deltamu"] = self.Compositional_PhaseDiagram3D.deltamu_values[self.first_element]
			self.DefectsDiagram.mu_elements[self.second_element]["deltamu"] = self.Compositional_PhaseDiagram3D.deltamu_values[self.second_element]
			self.DefectsDiagram.mu_elements[self.third_element]["deltamu"] = self.Compositional_PhaseDiagram3D.deltamu_values[self.third_element]
			self.DefectsDiagram.mu_elements[self.fourth_element]["deltamu"] = self.Compositional_PhaseDiagram3D.deltamu_values[self.fourth_element]
			self.DefectsDiagram.EVBM = self.compounds_info[main_compound]["vbm"]
			self.DefectsDiagram.ECBM = self.DefectsDiagram.EVBM + self.compounds_info[main_compound]["band_gap"]
			self.DefectsDiagram.fermi_energy_array = np.linspace(self.DefectsDiagram.EVBM, self.DefectsDiagram.ECBM, 100)
			self.DefectsDiagram.Activate_DefectsDiagram_Plot_Axes()
			# Update defects diagram in response to clicking on phase diagram
			self.Compositional_PhaseDiagram3D.composition_phasediagram_plot_figure.canvas.mpl_connect('button_press_event', self.Generate_DefectsDiagram_Plot_Function)
		
		
		
		
		
		###############################################################################################
		###############################################################################################
		#################################### Initialize third tab #####################################
		###############################################################################################
		###############################################################################################
		
		self.tab3 = QWidget()
		self.tab3_layout = QHBoxLayout(self.tab3)
		
		# Add compositional phase diagram window widget to tab3
		self.Compositional_PhaseDiagram_Window = Window_Compositional_PhaseDiagram(self.main_compound, self.Compositional_PhaseDiagram3D)
		#self.Compositional_PhaseDiagram_Window.compositional_phasediagram_window_layout.addWidget(self.Compositional_PhaseDiagram_Window.compound_title)
		self.Compositional_PhaseDiagram_Window.compositional_phasediagram_window_layout.addWidget(self.Compositional_PhaseDiagram_Window.composition_phase_diagram_title)
		self.Compositional_PhaseDiagram_Window.compositional_phasediagram_window_layout.addWidget(self.Compositional_PhaseDiagram_Window.composition_phase_diagram_plot)
		self.Compositional_PhaseDiagram_Window.compositional_phasediagram_window_layout.addWidget(self.Compositional_PhaseDiagram_Window.phasediagram_savefigure_button)
		self.tab3_layout.addWidget(self.Compositional_PhaseDiagram_Window.compositional_phasediagram_window)
		
		
		# Add defect formation energy diagram window widget to tab3
		if show_defects_diagram:
			self.DefectsDiagram_Window = Window_DefectsDiagram(self.main_compound, self.DefectsDiagram, show_carrier_concentration=False, show_dopant=False)
			self.tab3_layout.addWidget(self.DefectsDiagram_Window.defectsdiagram_window)
	
	
	
	
	###############################################################################################
	################################# Generate Defects Diagram ####################################
	###############################################################################################
	
	def Generate_DefectsDiagram_Plot_Function(self, event):
		
		# Check that a centroid of one of the four-phase regions has been clicked on
		contains, index = self.Compositional_PhaseDiagram3D.centroids_plot.contains(event)
		if not contains:
			return
		
		# Update elements and chemical potentials
		self.DefectsDiagram.mu_elements[self.first_element]["deltamu"] = self.Compositional_PhaseDiagram3D.deltamu_values[self.first_element]
		self.DefectsDiagram.mu_elements[self.second_element]["deltamu"] = self.Compositional_PhaseDiagram3D.deltamu_values[self.second_element]
		self.DefectsDiagram.mu_elements[self.third_element]["deltamu"] = self.Compositional_PhaseDiagram3D.deltamu_values[self.third_element]
		self.DefectsDiagram.mu_elements[self.fourth_element]["deltamu"] = self.Compositional_PhaseDiagram3D.deltamu_values[self.fourth_element]
		
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





