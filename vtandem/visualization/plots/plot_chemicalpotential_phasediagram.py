
__author__ = 'Michael_Lidia_Jiaxing_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'


import numpy as np
import itertools
import copy
import periodictable
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

import PyQt5
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

from vtandem.visualization.utils.chemicalpotential_phasediagram import Calculate_PhaseDiagram_Projected2D

from vtandem.visualization.utils.compound_name import Compound_Name_Formal

from vtandem.visualization.plots.save_plot import SaveFigure


class ChemicalPotential_PhaseDiagramProjected2D(QWidget, SaveFigure):
	
	def __init__(self, main_compound = None, type = None):
		
		# Check legitimacy of "type" argument
		if (type != "ternary") and (type != "quaternary"):
			raise Exception("'type' argument must be either 'ternary' or 'quaternary'. Skipping...")
		self.type = type
		
		# All elements in the periodic table
		self.all_elements = []
		for element in periodictable.elements:
			self.all_elements.append(str(element))
		
		# Font description for phase stability diagram plot
		self.font = {'family': 'sans-serif', 'color':  'black',	'weight': 'normal',	'size': 12 }
		
		# Establish the first, second, and third species of the quaternary compound.
		# Note that this list is subject to change, depending on what the user chooses.
		self.main_compound  = main_compound
		
		# Store all extracted DFT data
		self.main_compound_info = {}
		self.compounds_info = {}
		
		# Initialize necessary objects in app
		self.main_compound_plot = None
		self.competing_compound_plots = {}
		self.phase_stability_region = None
		self.PSR_vertices = []
		self.PSR_vertices_plot = None
		
		
		# Phase diagram plot
		self.phase_diagram_plot_figure = plt.figure()
		self.phase_diagram_plot_figure.subplots_adjust(left=0.0, bottom=0.0, right=0.8, top=0.8)
		self.phase_diagram_plot_drawing = self.phase_diagram_plot_figure.add_subplot(111)
		self.phase_diagram_plot_canvas = FigureCanvas(self.phase_diagram_plot_figure)
		
		
		# Save figure feature
		SaveFigure.__init__(self, self.phase_diagram_plot_figure)
	
	
	
	###############################################################################################
	##################################### Phase Diagram ###########################################
	###############################################################################################
	def Update_PhaseDiagram_Plot_Axes(self):
		
		self.phase_diagram_plot_drawing.set_xlim(self.phasediagram_endpoints, 0.0)
		self.phase_diagram_plot_drawing.set_ylim(self.phasediagram_endpoints, 0.0)
		self.phase_diagram_plot_drawing.set_xlabel("$\Delta\mu_{"+self.first_element+"}$ (eV)",fontdict=self.font)
		self.phase_diagram_plot_drawing.set_ylabel("$\Delta\mu_{"+self.second_element+"}$ (eV)",fontdict=self.font,rotation=270,labelpad=20)
		self.phase_diagram_plot_drawing.xaxis.tick_top()
		self.phase_diagram_plot_drawing.yaxis.tick_right()
		self.phase_diagram_plot_drawing.tick_params(axis='both', labelsize=self.font['size']-2)
		self.phase_diagram_plot_drawing.xaxis.set_label_position("top")
		self.phase_diagram_plot_drawing.yaxis.set_label_position("right")
		self.phase_diagram_plot_drawing.spines['left'].set_visible(False)
		self.phase_diagram_plot_drawing.spines['bottom'].set_visible(False)
		self.phase_diagram_plot_drawing.set_aspect("equal")
	
	
	
	def Plot_PhaseDiagram(self):
		
		if self.type == "ternary":
			elements_dict = {1: self.first_element, 2: self.second_element, 3: self.third_element}
		elif self.type == "quaternary":
			elements_dict = {1: self.first_element, 2: self.second_element, 3: self.third_element, 4: self.fourth_element}
		
		main_compound_deltamu_first_element, main_compound_stability_limit, \
			competing_compounds_deltamu_first_element_limit, competing_compounds_deltamu_second_element_limit, \
			main_compound_deltamu_first_element_cutoff, stability_minimum_cutoff, stability_maximum_cutoff \
			= Calculate_PhaseDiagram_Projected2D(self.main_compound, elements_dict, self.compounds_info, self.deltamu, self.main_compound_info)
		
		try:
			self.main_compound_plot.set_data(main_compound_deltamu_first_element, main_compound_stability_limit)
		except:
			self.main_compound_plot, = self.phase_diagram_plot_drawing.plot(main_compound_deltamu_first_element, main_compound_stability_limit, color='k')
		
		for competing_compound in self.compounds_info.keys():
			
			# Skip if compound is either the main compound or one of the elements
			if (competing_compound in self.all_elements) or (competing_compound == self.main_compound):
				continue
			
			try:
				# See if plot exists so it just needs to be updated
				self.competing_compound_plots[competing_compound].set_data(competing_compounds_deltamu_first_element_limit[competing_compound], competing_compounds_deltamu_second_element_limit[competing_compound])
			except:
				# If the plot doesn't exist initially, then create it
				self.competing_compound_plots[competing_compound], = self.phase_diagram_plot_drawing.plot(competing_compounds_deltamu_first_element_limit[competing_compound], competing_compounds_deltamu_second_element_limit[competing_compound], label=Compound_Name_Formal(competing_compound, "latex"))
		
		
		# Phase stability region
		try:
			self.phase_stability_region.remove()
			self.phase_stability_region = self.phase_diagram_plot_drawing.fill_between(main_compound_deltamu_first_element_cutoff, stability_maximum_cutoff, stability_minimum_cutoff, facecolor='0.75')
		except:
			self.phase_stability_region = self.phase_diagram_plot_drawing.fill_between(main_compound_deltamu_first_element_cutoff, stability_maximum_cutoff, stability_minimum_cutoff, facecolor='0.75')
		
		
		# Legend
		self.phase_diagram_plot_drawing.legend(loc=3, fontsize=self.font['size']-2)
		
		
		
		# Phase stability region vertices
		self.PSR_vertices = Find_PhaseStabilityRegion_Vertices(self.phase_stability_region)
		if self.PSR_vertices == []:
			try:
				self.PSR_vertices_plot.remove()
			except:
				pass
		elif self.PSR_vertices != []:
			try:
				self.PSR_vertices_plot.remove()
				self.PSR_vertices_plot = self.phase_diagram_plot_drawing.scatter(*zip(*self.PSR_vertices), s=20, c='black')
			except:
				self.PSR_vertices_plot = self.phase_diagram_plot_drawing.scatter(*zip(*self.PSR_vertices), s=20, c='black')
				pass
		
		self.phase_diagram_plot_canvas.draw()










def Find_PhaseStabilityRegion_Vertices(phase_stability_region):
	
	# This function finds the vertices bounding the phase stability region. It takes
	#	the points of the phase stability region as input.
	
	PSR_Vertices_Unrepeated = []
	
	if phase_stability_region.get_paths() != []:
		
		PSR_Vertices = []
		PSR_Bound_Slope_Previous = None
		PSR_Bounding_Point_Previous = None
		tolerance = 1E-6
		for PSR_Bounds in phase_stability_region.get_paths()[0].iter_segments():
			PSR_Bounding_Point = PSR_Bounds[0]
			try:
				PSR_Bound_Slope = (PSR_Bounding_Point[1]-PSR_Bounding_Point_Previous[1]) / (PSR_Bounding_Point[0]-PSR_Bounding_Point_Previous[0])
			except:
				PSR_Bounding_Point_Previous = PSR_Bounding_Point
				PSR_Bound_Slope_Previous = 0.0
				continue
			if (PSR_Bound_Slope < PSR_Bound_Slope_Previous - tolerance) or (PSR_Bound_Slope > PSR_Bound_Slope_Previous + tolerance):
				PSR_Vertices.append(PSR_Bounding_Point_Previous)
				PSR_Vertices.append(PSR_Bounding_Point)
			PSR_Bound_Slope_Previous = PSR_Bound_Slope
			PSR_Bounding_Point_Previous = PSR_Bounding_Point
		PSR_Vertices_Omit = []
		for PSR_Vertices_Index in range(len(PSR_Vertices)-1):
			if (np.linalg.norm(PSR_Vertices[PSR_Vertices_Index]-PSR_Vertices[PSR_Vertices_Index+1]) < 0.01):
				PSR_Vertices_Omit.append(PSR_Vertices[PSR_Vertices_Index])
		PSR_Vertices_Unrepeated = [x for x in PSR_Vertices if (not any((x is y for y in PSR_Vertices_Omit)))]
	
	return PSR_Vertices_Unrepeated













