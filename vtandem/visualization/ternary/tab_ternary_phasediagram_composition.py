
__author__ = 'Michael_Lidia_Jiaxing_Benita_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'


import numpy as np
import matplotlib.pyplot as plt

import PyQt5
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

from vtandem.visualization.ternary.ternary_scripts.plot_composition_ternary_phase_diagram import Composition_Ternary_PhaseDiagram
from vtandem.visualization.ternary.ternary_scripts.plot_ternary_defects_diagram import Ternary_Defects_Diagram

from vtandem.visualization.widgets.phasediagram_window_composition import Compositional_PhaseDiagram_Window



class Tab_Ternary_Compositional_PhaseDiagram(object):
	
	def __init__(self, main_compound=None, first_element=None, second_element=None, third_element=None, compounds_info=None, defects_data=None, show_defects_diagram=True):	# User specifies the main compound and its constituents
		
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
		
		"""
		# Keep track of mu values of the species in the ternary compound
		self.deltamu_values = {}
		self.deltamu_values[first_element] = 0.0
		self.deltamu_values[second_element] = 0.0
		self.deltamu_values[third_element] = self.compounds_info[main_compound]["enthalpy"] / self.compounds_info[main_compound][self.third_element]
		"""
		
		###############################################################################################
		################################# Compositional phase diagram #################################
		###############################################################################################
		
		self.Compositional_PhaseDiagram = Composition_Ternary_PhaseDiagram(main_compound = self.main_compound, first_element = self.first_element, second_element = self.second_element, third_element = self.third_element, compounds_info = self.compounds_info)
		self.Compositional_PhaseDiagram.composition_phasediagram_plot_figure.canvas.mpl_connect('button_press_event', self.Generate_DefectsDiagram_Plot_Function)
		"""
		self.Compositional_PhaseDiagram.compounds_info = self.compounds_info
		self.Compositional_PhaseDiagram.Create_Compositional_PhaseDiagram()
		self.Compositional_PhaseDiagram.Plot_Compositional_PhaseDiagram()	# Draws the compositional phase diagram
		"""
		
		
		"""
		# Find all three-phase regions after compositional phase diagram is drawn.
		#	Store lists of three points constituting the three-phase regions into self.threephaseregions.
		#	Note that the function is called AFTER the phase diagram is drawn.
		self.threephaseregions = []
		self.Find_All_ThreePhaseRegions()
		
		# Plot selected three-phase region.
		#	Store list of the three points constituting the three-phase region in self.threephaseregion_selected.
		#	Store the polygon plot object in self.threephaseregion_shade.
		self.threephaseregion_selected = None
		self.threephaseregion_shade = None
		self.Compositional_PhaseDiagram.composition_phasediagram_plot_figure.canvas.mpl_connect('button_press_event', self.Shade_ThreePhaseRegion)
		self.Compositional_PhaseDiagram.composition_phasediagram_plot_figure.canvas.mpl_connect('button_press_event', self.Calculate_ChemicalPotentials_ThreePhaseRegion)
		self.Compositional_PhaseDiagram.composition_phasediagram_plot_figure.canvas.mpl_connect('button_press_event', self.Generate_DefectsDiagram_Plot_Function)
		"""
		
		
		# Defects diagram
		if show_defects_diagram:
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
		#################################### Initialize third tab #####################################
		###############################################################################################
		###############################################################################################
		
		self.tab3 = QWidget()
		self.tab3_layout = QHBoxLayout(self.tab3)
		
		# Add compositional phase diagram window widget to tab3
		self.Compositional_PhaseDiagram_Window = Compositional_PhaseDiagram_Window(self.main_compound, self.Compositional_PhaseDiagram)
		self.tab3_layout.addWidget(self.Compositional_PhaseDiagram_Window.compositional_phasediagram_window)
		
		
		# Add defect formation energy diagram window widget to tab3
		if show_defects_diagram:
			
			###### Main defects diagram window widget
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
			
			self.tab3_layout.addWidget(self.defectsdiagram_window)
	
	
	
	"""
	###############################################################################################
	################# Find All Three-Phase Regions in Compositional Phase Diagram #################
	###############################################################################################
	
	def Find_All_ThreePhaseRegions(self):
		
		# Obtain all three-phase regions.
		#	The data structure is as follows:
		#		lines --> List of arrays, each array is 2x2 for ternary (3x2 for quaternary, etc.), column vector represents point on phase diagram.
		#					ex: array([ [0.3, 0.5], [1.0, 0.0] ]) is a line that goes from point [x=0.3, y=1.0] to point [x=0.5, y=0.0]
		#		labels --> Dictionary with point-PDEntry pairs.
		#	Algorithm works as follows:
		#		1) Loop through all lines and get the endpoints.
		#		2) Loop through all points.
		#		3) If the endpoints both draw lines to the third point, then a three-phase region exists between the three points.
		for line in self.Compositional_PhaseDiagram.lines:
			endpoint1, endpoint2 = np.transpose(line)
			for point in self.Compositional_PhaseDiagram.labels.keys():
				point = np.asarray(point)
				if 		( any(np.array_equal(np.transpose(np.vstack((point, endpoint1))), i) for i in self.Compositional_PhaseDiagram.lines) or any(np.array_equal(np.transpose(np.vstack((endpoint1, point))), i) for i in self.Compositional_PhaseDiagram.lines) ) \
					and ( any(np.array_equal(np.transpose(np.vstack((point, endpoint2))), i) for i in self.Compositional_PhaseDiagram.lines) or any(np.array_equal(np.transpose(np.vstack((endpoint2, point))), i) for i in self.Compositional_PhaseDiagram.lines) ):
					
					# Check if three phase region has already been recorded
					if any( all( list(x) in threephaseregion for x in [point, endpoint1, endpoint2] ) for threephaseregion in self.threephaseregions):
						continue
					
					# Record three phase region
					three_phase_region = [list(point), list(endpoint1), list(endpoint2)]
					self.threephaseregions.append(three_phase_region)
	
	
	
	###############################################################################################
	############ Shade Three-Phase Region in Compositional Phase Diagram when Clicked #############
	###############################################################################################
	
	def Shade_ThreePhaseRegion(self, event):
		
		# Check to see that all three-phase regions have been calculated
		if self.threephaseregions is None:
			return
		
		# Read coordinates of clicked point
		point_x = event.xdata
		point_y = event.ydata
		
		# Check that the user presses somewhere on the plot (and not anywhere else)
		if (not isinstance(point_x, float)) or (not isinstance(point_y, float)):
			return
		
		# Shade in the selected three-phase region
		for threephaseregion in self.threephaseregions:
			is_in_threephaseregion = mpltPath.Path(threephaseregion).contains_point([point_x, point_y])
			if is_in_threephaseregion:
				self.threephaseregion_selected = threephaseregion
				if self.threephaseregion_shade is not None:
					self.threephaseregion_shade.remove()
				self.threephaseregion_shade = plt.Polygon(threephaseregion, color='gray', alpha=0.5)
				self.Compositional_PhaseDiagram.composition_phasediagram_plot_drawing.add_patch(self.threephaseregion_shade)
				self.Compositional_PhaseDiagram.composition_phasediagram_plot_canvas.draw()
	
	
	
	###############################################################################################
	############### Calculate Chemical Potentials of Elements in Three-Phase Region ###############
	###############################################################################################
	
	def Calculate_ChemicalPotentials_ThreePhaseRegion(self, event):
		
		# Check to see that a three-phase region has been selected
		if self.threephaseregion_selected is None:
			return
		
		# Get the names of compounds constituting the three-phase region
		threephaseregion_compounds = []
		for point in self.threephaseregion_selected:
			entry = self.Compositional_PhaseDiagram.labels[tuple(point)]	# PDEntry object
			threephaseregion_compounds.append( entry.name )
		
		# Get the matrix encoding the stoichiometry of each compound constituting the three-phase region
		composition_matrix = np.zeros((3, 3))
		for compound_index, compound in enumerate(threephaseregion_compounds):
			for element_index, element in enumerate(self.elements_list):
				try:
					composition_matrix[compound_index][element_index] = self.compounds_info[compound][element]
				except:
					composition_matrix[compound_index][element_index] = 0.0
		
		# Get the enthalpies of formation of each compound in the three-phase region
		enthalpies_array = np.zeros((3, 1))
		for compound_index, compound in enumerate(threephaseregion_compounds):
			enthalpies_array[compound_index] = self.compounds_info[compound]["enthalpy"]
		
		# Solve for the delta mu values
		deltamus = np.dot( np.linalg.inv(composition_matrix), enthalpies_array )
		print(deltamus)
		
		
		# Record delta mu values
		for element, deltamu in zip(self.elements_list, deltamus):
			self.deltamu_values[element] = deltamu
		self.DefectsDiagram.mu_elements[self.first_element]["deltamu"] = self.deltamu_values[self.first_element]
		self.DefectsDiagram.mu_elements[self.second_element]["deltamu"] = self.deltamu_values[self.second_element]
		self.DefectsDiagram.mu_elements[self.third_element]["deltamu"] = self.deltamu_values[self.third_element]
		print(self.deltamu_values)
		print(self.DefectsDiagram.mu_elements)
	"""
	
	
	
	###############################################################################################
	################################# Generate Defects Diagram ####################################
	###############################################################################################
	
	def Generate_DefectsDiagram_Plot_Function(self, event):
		
		# This function specifies what happens when the user clicks the "Generate Defects Diagram" button.
		
		# Update elements and chemical potentials
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
	
	
	
	
	
	
	"""
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
	"""
	
	
	
	
	
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
		"""
		if filename:
			extension = extension_type.split(".")[-1].split(")")[0]
			if filename.split(".")[-1] == extension:
				self.DefectsDiagram.defects_diagram_plot_figure.savefig(filename, bbox_inches='tight')
			else:
				self.DefectsDiagram.defects_diagram_plot_figure.savefig(filename+"."+extension, bbox_inches='tight')
		"""







