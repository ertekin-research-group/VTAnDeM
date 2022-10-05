
__author__ = 'Michael_Lidia_Jiaxing_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from pymatgen.core.composition import Composition

import PyQt5
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

from vtandem.visualization.plots.save_plot import SaveFigure

from vtandem.visualization.utils.compositional_phasediagram import *
from vtandem.visualization.utils.compound_name import Compound_Name_Formal


class Plot_Composition_PhaseDiagram(SaveFigure):
	
	def __init__(self, type: str, main_compound_info: dict, compounds_info: dict):
		
		# Record type of phase diagram (either ternary or quaternary)
		if type not in ["ternary", "quaternary"]:
			raise ValueError("Argument 'type' can be either 'ternary' or 'quaternary'.")
		self.type = type
		
		# Font description for phase stability diagram plot
		fontsize = 12
		"""
		self.font_labels = {'family': 'sans-serif', 'color': 'black', 'weight': 'bold', 'size': fontsize }
		self.font_legend = {'family': 'sans-serif', 'size': fontsize }
		self.font_annotation = {'family': 'sans-serif', 'color': 'black', 'weight': 'bold', 'size': fontsize }
		"""
		self.font_labels = {'color': 'black', 'weight': 'bold', 'size': fontsize }
		self.font_legend = {'size': fontsize }
		self.font_annotation = {'color': 'black', 'weight': 'bold', 'size': fontsize }
		
		# Store all extracted DFT data
		self.main_compound_info = main_compound_info
		self.compounds_info = compounds_info
		print(main_compound_info)

		# Phase diagram (in composition space) object
		self.composition_phasediagram_plot_figure = plt.figure()
		self.composition_phasediagram_plot_canvas = FigureCanvas(self.composition_phasediagram_plot_figure)


		# Initialize plot drawing
		if self.type == "ternary":
			self.composition_phasediagram_plot_drawing = self.composition_phasediagram_plot_figure.add_subplot(111)
		if self.type == "quaternary":
			self.composition_phasediagram_plot_drawing = self.composition_phasediagram_plot_figure.add_subplot(111, projection='3d')
		

		# Generate compositional phase diagram
		self.Generate_Compositional_PhaseDiagram(self.compounds_info, self.elements_list)
		
		# Save figure feature
		SaveFigure.__init__(self, self.composition_phasediagram_plot_figure)
	
	

	def Generate_Compositional_PhaseDiagram(self, compounds_info, elements_list):

		try:
			self.composition_phasediagram_legend.remove()
			self.composition_phasediagram_plot_drawing.clear()  # Doesn't remove the drawing, it simply "clears" it ;)
		except:
			pass
		
		# Generate compositional phase diagram
		#self.pmg_phasediagram, lines, self.labels = Create_Compositional_PhaseDiagram(compounds_info, elements_list)
		self.pmg_phasediagram, lines, self.labels = Create_Compositional_PhaseDiagram(compounds_info, elements_list, self.main_compound_info, self.main_compound)
		newlabels = Plot_Compositional_PhaseDiagram(self.composition_phasediagram_plot_drawing, 
													self.type,
													lines,
													self.labels,
													self.font_labels
													)
		self.composition_phasediagram_legend = self.composition_phasediagram_plot_figure.text(0.01, 0.01, "\n".join(newlabels), fontdict=self.font_legend)

		self.composition_phasediagram_plot_canvas.draw()


		### Find all phase regions after compositional phase diagram is drawn.
		self.phase_region_objects = []
		self.phaseregion_selected = None
		self.phaseregion_shade = []
		
		# Trigger for shading (MUST come before chemical potentials calculation)
		if self.type == "ternary":
			self.composition_phasediagram_plot_figure.canvas.mpl_connect('button_press_event', self.Shade_ThreePhaseRegion)
		if self.type == "quaternary":
			self.composition_phasediagram_plot_figure.canvas.mpl_connect('button_press_event', self.Shade_FourPhaseRegion)
		
		# Other triggers
		self.composition_phasediagram_plot_figure.canvas.mpl_connect('button_press_event', self.Calculate_ChemicalPotentials_Region)
		self.composition_phasediagram_plot_figure.canvas.mpl_connect('motion_notify_event', self.Hover)

		self.text_x_offset = 20
		self.text_y_offset = 20
		self.phaseregion_annotation = self.composition_phasediagram_plot_drawing.annotate("", xy=(0,0), xytext=(self.text_x_offset, self.text_y_offset), 
																								textcoords="offset pixels",
																								bbox = dict(boxstyle="round", fc="w"),
																								arrowprops = dict(arrowstyle = "->"),
																								fontsize = self.font_annotation['size'] )
		self.phaseregion_annotation.get_bbox_patch().set_facecolor("w")
		self.phaseregion_annotation.get_bbox_patch().set_alpha(1.0)




	###############################################################################################
	#################### Find All Phase Regions in Compositional Phase Diagram ####################
	###############################################################################################
	
	def Find_All_PhaseRegions(self):
		
		# Get all phase region names
		all_phase_regions = self.pmg_phasediagram.get_all_chempots(Composition(self.main_compound)).keys()
		
		# Get the coordinates of each stable phase in the phase diagram
		pd_coordinates = {}
		for coordinate in self.labels.keys():
			pd_coordinates[ self.labels[coordinate].composition.reduced_composition ] = coordinate

		# Create and store phase region objects
		for phase_region in all_phase_regions:
			
			phase_region_compounds = []
			for compound in phase_region.split("-"):
				phase_region_compounds.append(compound)
			
			phase_region_vertices = []
			for compound in phase_region_compounds:
				#phase_region_vertices.append( pd_coordinates[Composition(compound)] )
				phase_region_vertices.append( pd_coordinates[Composition(compound).reduced_composition] )
			
			centroid = np.sum( np.asarray(phase_region_vertices), axis=0 ) / len(phase_region_vertices)
			
			phase_region_obj = PhaseRegion()
			phase_region_obj.name = ', '.join(phase_region_compounds)
			phase_region_obj.vertices = phase_region_vertices
			phase_region_obj.centroid = centroid
			
			self.phase_region_objects.append(phase_region_obj)
	
	
	
	
	
	
	###############################################################################################
	########################## Plot Centroids of All Four-Phase Regions ###########################
	###############################################################################################
	
	def Plot_Centroids(self):
		
		# Check that all four-phase regions have been found
		if self.phase_region_objects == []:
			return
		
		# Set color of scatter plot points
		scatterplot_color = 'k'
		scatterplot_marker = '*'
		
		# Plot all centroids individually in scatter plot
		centroids = []
		for phase_region in self.phase_region_objects:
			centroid = phase_region.centroid
			centroids.append(centroid)
			phase_region.centroid_plot = self.composition_phasediagram_plot_drawing.scatter(	*zip(*centroids),
																								color = scatterplot_color,
																								marker = scatterplot_marker	)
		
		# Plot all centroids together in scatter plot
		centroids = np.asarray(centroids)
		self.centroids_plot = self.composition_phasediagram_plot_drawing.scatter(	*zip(*centroids),
																					color = scatterplot_color,
																					marker = scatterplot_marker )
	
	
	
	def Update_Annotation(self, index, event):
		
		# Get xy position of cursor on screen (in data coordinates)
		screen_position_xy = self.centroids_plot.get_offsets()[index]
		self.phaseregion_annotation.xy = screen_position_xy
		
		# Get name of four phase region
		for phase_region in self.phase_region_objects:
			if phase_region.centroid_plot.contains(event)[0]:
				phase_region_compound_names = []
				for compound_name in phase_region.name.split(","):
					phase_region_compound_names.append( Compound_Name_Formal(compound_name, "latex") )
				text = ", ".join(phase_region_compound_names)
				break
		
		# Set annotation
		self.phaseregion_annotation.set_text(text)


		# Check if annotation spills out of axis
		axis_bounds = self.composition_phasediagram_plot_drawing.get_window_extent().bounds # Bounds are given in (x0, y0, width, height) format
		axis_x1, axis_y1 = axis_bounds[0]+axis_bounds[2], axis_bounds[1]+axis_bounds[3]
		annotation_bounds = self.phaseregion_annotation.get_window_extent().bounds
		annotation_x1, annotation_y1 = annotation_bounds[0]+annotation_bounds[2], annotation_bounds[1]+annotation_bounds[3]
		if (annotation_x1 > axis_x1):
			self.phaseregion_annotation.set( position=(-axis_bounds[2]/2, 20) )
	
	
	def Hover(self, event):

		is_visible = self.phaseregion_annotation.get_visible()
		
		if event.inaxes == self.composition_phasediagram_plot_drawing:
			
			contains, index = self.centroids_plot.contains(event)
			if contains:
				self.Update_Annotation(index["ind"][0], event)
				self.phaseregion_annotation.set_visible(True)
				self.composition_phasediagram_plot_canvas.draw()
			else:
				if is_visible:
					self.phaseregion_annotation.set_visible(False)
					self.composition_phasediagram_plot_canvas.draw()
	
	
	
	
	
	
	
	
	
	
	
	###############################################################################################
	################## Calculate Chemical Potentials of Elements in Phase Region ##################
	###############################################################################################
	
	def Calculate_ChemicalPotentials_Region(self, event):
		
		print("Calculating chemical potentials...")
		
		# Check to see that a phase region has been selected
		if self.phaseregion_selected is None:
			return
		
		# Get the matrix encoding the stoichiometry of each compound constituting the phase region
		composition_matrix = []
		enthalpies_array = []
		for compound in self.phaseregion_selected.name.replace(" ", "").split(","):
			
			# Get compound stoichiometry and enthalpy
			compound_stoichiometry = []
			compound_enthalpy = self.compounds_info[compound]["dft_total_energy"]
			
			for element in self.elements_list:
				try:
					element_count = self.compounds_info[compound]["dft_"+element]
				except:
					element_count = 0.0
				compound_stoichiometry.append( element_count )
				compound_enthalpy -= element_count * self.compounds_info[element]["mu0"]
			
			composition_matrix.append(compound_stoichiometry)
			enthalpies_array.append(compound_enthalpy)
		
		composition_matrix = np.asarray(composition_matrix)
		enthalpies_array = np.asarray(enthalpies_array)
		
		# Solve for the delta mu values
		deltamus = np.dot( np.linalg.inv(composition_matrix), enthalpies_array )
		print(self.elements_list)
		print(deltamus)
		print(self.phaseregion_selected.name)
		
		# Record delta mu values
		for element, deltamu in zip(self.elements_list, deltamus):
			self.deltamu_values[element] = deltamu





class PhaseRegion(object):
	def __init__(self):
		self.name = None
		self.vertices = None
		self.centroid = None
		self.centroid_plot = None






