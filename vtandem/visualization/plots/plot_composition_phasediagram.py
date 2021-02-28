
__author__ = 'Michael_Lidia_Jiaxing_Benita_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'


import numpy as np
import itertools
import periodictable
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter, PDEntry
from pymatgen.core.composition import Composition

import PyQt5
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *


from vtandem.visualization.plots.save_plot import SaveFigure



class Plot_Composition_PhaseDiagram(SaveFigure):
	
	def __init__(self, type: str, main_compound_info: dict, compounds_info: dict):
		
		# Record type of phase diagram (either ternary or quaternary)
		if type not in ["ternary", "quaternary"]:
			raise ValueError("Argument 'type' can be either 'ternary' or 'quaternary'.")
		self.type = type
		
		# All elements in the periodic table
		self.all_elements = []
		for element in periodictable.elements:
			self.all_elements.append(str(element))
		
		# Font description for phase stability diagram plot
		self.font = {'family': 'sans-serif', 'color':  'black', 'weight': 'bold', 'size': 10 }
		
		# Store all extracted DFT data
		self.main_compound_info = main_compound_info
		self.compounds_info = compounds_info
		
		# Initialize pymatgen phase diagram objects
		self.pmg_phasediagram = None
		self.pmg_phasediagram_plot_object = None
		
		# Phase diagram (in composition space) object
		self.composition_phasediagram_plot_figure = plt.figure()
		self.composition_phasediagram_plot_canvas = FigureCanvas(self.composition_phasediagram_plot_figure)
		if self.type == "ternary":
			self.composition_phasediagram_plot_drawing = self.composition_phasediagram_plot_figure.add_subplot(111)
		if self.type == "quaternary":
			self.composition_phasediagram_plot_drawing = self.composition_phasediagram_plot_figure.add_subplot(111, projection='3d')
		
		# Save figure feature
		SaveFigure.__init__(self, self.composition_phasediagram_plot_figure)
		
		# Store all lines and labels of the compositional phase diagram
		self.lines = []
		self.labels = {}
		
		# Store all plots of phase diagram vertices
		self.vertices_plots = {}
		
		
		# Generate compositional phase diagram
		self.Create_Compositional_PhaseDiagram()
		self.Plot_Compositional_PhaseDiagram()
		
		
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
		
		self.phaseregion_annotation = self.composition_phasediagram_plot_drawing.annotate("", xy=(0,0), xytext=(20,20), textcoords="offset points",
																								bbox = dict(boxstyle="round", fc="w"),
																								arrowprops = dict(arrowstyle = "->") )
	
	
	
	
	
	
	def Create_Compositional_PhaseDiagram(self):
		
		# Record all entries for the phase diagram
		phasediagram_entries = []
		
		for compound in self.compounds_info.keys():
			
			# Disregard elements not included in main compound
			if (compound in self.all_elements) and (compound not in self.elements_list):
				continue
			
			# Get the compound's composition
			compound_composition = {}
			if compound in self.elements_list:	# Elemental material
				compound_composition[compound] = self.compounds_info[compound]["dft_"+compound]
			else:	# Compound material
				for element in self.compounds_info[compound]["elements_list"]:
					compound_composition[element] = self.compounds_info[compound]["dft_"+element]
			
			# Get the compound's total energy
			compound_total_energy = self.compounds_info[compound]["dft_total_energy"]
			
			# Record to list of entries
			#phasediagram_entries.append(PDEntry(compound_composition, compound_total_energy, compound))
			#phasediagram_entries.append(PDEntry(compound_composition, compound_total_energy))
			phasediagram_entries.append(PDEntry(composition=Composition(compound_composition), energy=compound_total_energy, name=compound))
		
		# Calculate compositional phase diagram (using pymatgen)
		#	The output data structure is as follows:
		#		lines --> List of arrays, each array is 2x2 for ternary (3x3 for quaternary, etc.), column vector represents point on phase diagram.
		#					ex: array([ [0.3, 0.5], [1.0, 0.0] ]) is a line that goes from point [x=0.3, y=1.0] to point [x=0.5, y=0.0]
		#		labels --> Dictionary with point-PDEntry pairs.
		self.pmg_phasediagram = PhaseDiagram(phasediagram_entries)
		self.pmg_phasediagram_plot_object = PDPlotter(self.pmg_phasediagram)
		(lines, labels, unstable) = self.pmg_phasediagram_plot_object.pd_plot_data
		
		# Record all lines and points of the compositional phase diagram
		self.lines = lines
		self.labels = labels
	
	
	
	
	
	def Plot_Compositional_PhaseDiagram(self, label_stable=True):
		
		# Plot settings
		linewidth = 1
		markersize = 5
		
		# Plot compositional phase diagram
		count = 1
		newlabels = []
		if self.type == "ternary":
			for x, y in self.lines:
				self.composition_phasediagram_plot_drawing.plot(x, y, "bo-", linewidth=linewidth, markeredgecolor="b", markerfacecolor="r", markersize=markersize)
			for coords in sorted(self.labels.keys()):
				entry = self.labels[coords]
				label = entry.name
				if label_stable:
					if len(entry.composition.elements) == 1:
						self.composition_phasediagram_plot_drawing.text(coords[0], coords[1], label, fontdict=self.font)
					else:
						self.composition_phasediagram_plot_drawing.text(coords[0], coords[1], str(count), fontdict=self.font)
						newlabels.append("{} : {}".format(count, label))
						count += 1
		elif self.type == "quaternary":
			for x, y, z in self.lines:
				self.composition_phasediagram_plot_drawing.plot(x, y, z, "bo-", linewidth=linewidth, markeredgecolor="b", markerfacecolor="r", markersize=markersize)
			for coords in sorted(self.labels.keys()):
				entry = self.labels[coords]
				label = entry.name
				if label_stable:
					if len(entry.composition.elements) == 1:
						self.composition_phasediagram_plot_drawing.text(coords[0], coords[1], coords[2], label, fontdict=self.font)
					else:
						self.composition_phasediagram_plot_drawing.text(coords[0], coords[1], coords[2], str(count), fontdict=self.font)
						newlabels.append("{} : {}".format(count, label))
						count += 1
		
		# Draw compositional phase diagram
		self.composition_phasediagram_legend = self.composition_phasediagram_plot_figure.text(0.01, 0.01, "\n".join(newlabels))
		self.composition_phasediagram_plot_drawing.axis("off")
		self.composition_phasediagram_plot_canvas.draw()
	
	
	
	
	
	
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
				phase_region_vertices.append( pd_coordinates[Composition(compound)] )
			
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
		
		# Get xy position of cursor on screen
		screen_position_xy = self.centroids_plot.get_offsets()[index]
		self.phaseregion_annotation.xy = screen_position_xy
		
		# Get name of four phase region
		for phase_region in self.phase_region_objects:
			if phase_region.centroid_plot.contains(event)[0]:
				text = phase_region.name
				break
		
		# Set annotation
		self.phaseregion_annotation.set_text(text)
		self.phaseregion_annotation.get_bbox_patch().set_facecolor("w")
		self.phaseregion_annotation.get_bbox_patch().set_alpha(1.0)
	
	
	
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
		
		print(self.phaseregion_selected.name)
		
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
		print(deltamus)
		
		# Record delta mu values
		for element, deltamu in zip(self.elements_list, deltamus):
			self.deltamu_values[element] = deltamu





class PhaseRegion(object):
	def __init__(self):
		self.name = None
		self.vertices = None
		self.centroid = None
		self.centroid_plot = None






