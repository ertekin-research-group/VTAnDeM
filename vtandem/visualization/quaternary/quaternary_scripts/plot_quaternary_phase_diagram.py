
__author__ = 'Michael_Lidia_Jiaxing_Benita_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'


###############################################################################################################################
###############################################################################################################################
##################################################### Introduction ############################################################
###############################################################################################################################
###############################################################################################################################

# 	This script stores the phase stability diagram of a single quaternary compound as an object within VTAnDeM.
#
# 	The stability of any compound is dictated by thermodynamic principles. All materials are made of atoms, and multicomponent
#		systems in particular (e.g. systems of two or more types of atoms, like GaAs) can exhibit distinct phases/materials
#		depending on how the atoms are arranged. Therefore, for a given set of atoms, there is a myriad of possible materials 
#		that can be formed by e.g. different synthesis methods. Fundamentally, the material with the lowest formation energy in some
#		environmental condition (e.g. GaAs in a 'Ga-rich' environment) is most and will likely form in a lab. Accordingly, if
#		we are studying some quaternary system, there may be other compounds that "compete" with the main compound to become the
#		most stable. For example, if we are interested in studying the quaternary compound Cu2HgGeTe4, in regions that are Cu-rich
#		or Hg-poor (or any combination of the sorts), another compound such as CuTe may be more stable than Cu2HgGeTe4. It is 
#		therefore crucial for material/device processing purposes that we know the exact environmental conditions necessary to create a 
#		stable version of the quaternary compound of interest. The stability of compounds can fortunately be studied using Density 
#		Functional Theory (DFT) calculations.



###############################################################################################################################
###############################################################################################################################
################################################### Import Libraries ##########################################################
###############################################################################################################################
###############################################################################################################################

import numpy as np
import itertools
import copy
import periodictable
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import PyQt5
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *


class ChemicalPotential_Quaternary_PhaseDiagramProjected2D(QWidget):
	
	def __init__(self, parent = None, main_compound = None, first_element = None, second_element = None, third_element = None, fourth_element = None):
		
		# All elements in the periodic table
		self.all_elements = []
		for element in periodictable.elements:
			self.all_elements.append(str(element))
		
		#QWidget.__init__(self, parent)
		
		# Font description for phase stability diagram plot
		self.font = {'family': 'sans-serif',
				'color':  'black',
				'weight': 'normal',
				'size': 16 }
		
		# Establish the first, second, third, and fourth species of the quaternary compound.
		# Note that this list is subject to change, depending on what the user chooses.
		self.main_compound  = main_compound
		self.first_element	= first_element
		self.second_element	= second_element
		self.third_element	= third_element
		self.fourth_element	= fourth_element
		self.elements_list  = [self.first_element, self.second_element, self.third_element, self.fourth_element]	# List of elements that gets updated as user selects element order
		self.elements_list_original = [self.first_element, self.second_element, self.third_element, self.fourth_element]	# Unchanged list of elements for compound naming purposes (e.g. see Compound_Name_Formal)
		
		# Store all extracted DFT data
		self.compounds_info = {}
		
		# Necessary numerical variables
		self.main_compound_number_first_specie  = 0	# Number of each specie in the main quaternary compound
		self.main_compound_number_second_specie = 0
		self.main_compound_number_third_specie  = 0
		self.main_compound_number_fourth_specie = 0
		self.main_compound_enthalpy = 0.0			# Enthalpy of the main quaternary compound
		self.phasediagram_endpoints = 0.0			# Endpoints for quaternary phase diagram
		self.mu4 = 0.0								# Track the fourth species mu value as the user changes it
		
		# Quaternary phase diagram object
		self.quaternary_phase_diagram_plot_figure = plt.figure()
		#self.quaternary_phase_diagram_plot_figure.subplots_adjust(left=0.001, right=0.85, bottom=0.001, top=0.85)
		self.quaternary_phase_diagram_plot_figure.subplots_adjust(left=0.0, bottom=0.0, right=0.8, top=0.8)
		self.quaternary_phase_diagram_plot_drawing = self.quaternary_phase_diagram_plot_figure.add_subplot(111)
		self.quaternary_phase_diagram_plot_canvas = FigureCanvas(self.quaternary_phase_diagram_plot_figure)
		
		#self.toolbar = NavigationToolbar(self.quaternary_phase_diagram_plot_canvas, self)
		
		# Necessary plot-related variables
		self.main_compound_plot = None				# Main compound plot (holds plot object)
		self.competing_compound_plots = {}			# Competing compound plots (holds plot objects)
		self.competing_compound_colorwheel = {}		# Legend
		self.PSR_vertices_plot = None				# Dots for intersections
		self.PSR_vertices = []						# Coordinates of each intersection
	
	
	
	###############################################################################################
	########################## Rewrite Compound Name Latex-Style ##################################
	###############################################################################################
	
	def Compound_Name_Formal(self, compound_name):
		
		# Go into the compounds_info dictionary and obtains the chemistry and stoichiometry of the compound of choice
		compound_species_info = self.compounds_info[compound_name]
		
		compound_name_formal = ""
		for species in compound_species_info["elements_list"]:				# Loop through the list of possible species that can be contained in the compound
			if compound_species_info[species] == 0:
				continue								# Don't add the species to the name if the compound doesn't contain the species
			elif compound_species_info[species] == 1:
				compound_name_formal += species			# Add the species to the name of the compound
			elif compound_species_info[species] > 1:
				compound_name_formal += species+"$_{"+str(int(compound_species_info[species]))+"}$"	# Add the species to the name of the compound
				#compound_name_formal += species+"<sub>"+str(int(compound_species_info[species]))+"</sub>"	# Add the species to the name of the compound
																											#	with a subscript for the stoichiometry
		
		return compound_name_formal
	
	
	
	def Update_PhaseDiagram_Object(self):
		
		# Number of elements in main compound
		self.main_compound_number_first_specie = self.compounds_info[self.main_compound][self.first_element]	# Number of first species in quaternary compound
		self.main_compound_number_second_specie = self.compounds_info[self.main_compound][self.second_element]	# Number of second species in quaternary compound
		self.main_compound_number_third_specie = self.compounds_info[self.main_compound][self.third_element]	# Number of third species in quaternary compound
		self.main_compound_number_fourth_specie = self.compounds_info[self.main_compound][self.fourth_element]	# Number of fourth species in quaternary compound
		
		# Enthalpy of quaternary compound
		self.main_compound_enthalpy = self.compounds_info[self.main_compound]["enthalpy"]	# Enthalpy of quaternary compound
		
		# Endpoints of phase diagram
		self.phasediagram_endpoints = min(self.main_compound_enthalpy/self.main_compound_number_first_specie, self.main_compound_enthalpy/self.main_compound_number_second_specie, self.main_compound_enthalpy/self.main_compound_number_third_specie, self.main_compound_enthalpy/self.main_compound_number_fourth_specie)
	
	
	def Update_PhaseDiagram_Plot_Axes(self):
		
		self.quaternary_phase_diagram_plot_drawing.set_xlim(self.phasediagram_endpoints, 0.0)
		self.quaternary_phase_diagram_plot_drawing.set_ylim(self.phasediagram_endpoints, 0.0)
		self.quaternary_phase_diagram_plot_drawing.set_xlabel("$\Delta\mu_{"+self.first_element+"}$ (eV)",fontdict=self.font)
		self.quaternary_phase_diagram_plot_drawing.set_ylabel("$\Delta\mu_{"+self.second_element+"}$ (eV)",fontdict=self.font,rotation=270,labelpad=20)
		self.quaternary_phase_diagram_plot_drawing.xaxis.tick_top()
		self.quaternary_phase_diagram_plot_drawing.yaxis.tick_right()
		#self.quaternary_phase_diagram_plot_drawing.tick_params(axis='both', labelsize=9)
		self.quaternary_phase_diagram_plot_drawing.xaxis.set_label_position("top")
		self.quaternary_phase_diagram_plot_drawing.yaxis.set_label_position("right")
		self.quaternary_phase_diagram_plot_drawing.spines['left'].set_visible(False)
		self.quaternary_phase_diagram_plot_drawing.spines['bottom'].set_visible(False)
		self.quaternary_phase_diagram_plot_drawing.set_aspect("equal")
	
	
	
	"""
	def Update_PhaseDiagram_Plot_Axes(self):
		
		# This function changes the axes labels of the quaternary phase diagram when the user selects
		#	the species to be plotted. This will be used in conjunction with the 
		#	Activate_PhaseDiagram_Plot_Axes function, which creates the initial phase diagram plot.
		
		self.quaternary_phase_diagram_plot_drawing.set_xlabel("$\Delta\mu_{"+self.first_element+"}$ (eV)",fontdict=self.font,labelpad=12)
		self.quaternary_phase_diagram_plot_drawing.set_ylabel("$\Delta\mu_{"+self.second_element+"}$ (eV)",fontdict=self.font,rotation=270,labelpad=20)
	"""
	
	def Plot_PhaseDiagram(self):
		
		# Plotting the phase diagram from DFT data is one of the main features of this app. This function
		#	single-handedly plots the phase diagram, so it's arguably one of the most important block of
		#	code in this script.
		# This function will be used in two forms:
		#	1) when the user generates the phase diagram by clicking the "Generate Plot" button
		#	2) when the user changes the mu4 value
		# Because there are different forms in which this function will be used, the programmer should take
		#	care when editing this block of code.
		
		main_compound_enthalpy_adjusted = self.main_compound_enthalpy - self.main_compound_number_fourth_specie*self.mu4	# Enthalpy adjusted for mu4 value (main compound)
		
		main_compound_deltamu_first_element = np.linspace(self.main_compound_enthalpy/self.main_compound_number_first_specie, 0, 1000)	# Array of possible mu1 values in the plot
		
		main_compound_stability_limit = (main_compound_enthalpy_adjusted - self.main_compound_number_first_specie*main_compound_deltamu_first_element) / self.main_compound_number_second_specie	# Array that represents the (diagonal) stability limit of the main compound
		
		# Check if the phase diagram has already been plotted
		try:
			self.main_compound_plot.set_ydata(main_compound_stability_limit)	# If the plot's already there, simply update the stability limit of the main compound
		except:
			self.main_compound_plot, = self.quaternary_phase_diagram_plot_drawing.plot(main_compound_deltamu_first_element, main_compound_stability_limit)	# If the main compound stabilitiy limit hasn't been plotted yet, plot it
		
		# Find the bounds of the quaternary phase stability region (for shading the stability region of the quaternary compound)
		stability_minimum_bound = []
		stability_maximum_bound = []
		
		vertical_left_values = []
		vertical_right_values = []
		
		# Loop through all compounds in the database
		for competing_compound in self.compounds_info.keys():
			
			# Skip if compound is either the main compound or one of the elements
			if (competing_compound in self.all_elements) or (competing_compound == self.main_compound):
				continue
			
			try:
				competing_compound_number_first_specie = self.compounds_info[competing_compound][self.first_element]
			except:
				competing_compound_number_first_specie = 0.0
				pass
			
			try:
				competing_compound_number_second_specie = self.compounds_info[competing_compound][self.second_element]
			except:
				competing_compound_number_second_specie = 0.0
				pass
			
			try:
				competing_compound_number_third_specie = self.compounds_info[competing_compound][self.third_element]
			except:
				competing_compound_number_third_specie = 0.0
				pass
			
			try:
				competing_compound_number_fourth_specie = self.compounds_info[competing_compound][self.fourth_element]
			except:
				competing_compound_number_fourth_specie = 0.0
				pass
			
			competing_compound_enthalpy = self.compounds_info[competing_compound]["enthalpy"]
			competing_compound_enthalpy_adjusted = competing_compound_enthalpy - competing_compound_number_fourth_specie*self.mu4		# Enthalpy adjusted for mu4 value (competing compound)
			
			difference_enthalpy_adjusted = competing_compound_enthalpy_adjusted - (competing_compound_number_third_specie / self.main_compound_number_third_specie) * main_compound_enthalpy_adjusted
			
			coefficient_first_specie = competing_compound_number_first_specie - (self.main_compound_number_first_specie * competing_compound_number_third_specie) / self.main_compound_number_third_specie
			coefficient_second_specie = competing_compound_number_second_specie - (self.main_compound_number_second_specie * competing_compound_number_third_specie) / self.main_compound_number_third_specie
			
			if (coefficient_first_specie == 0.0) and (coefficient_second_specie == 0.0):
				continue
				print("OMG SOMETHING'S WRONG!!!")
				print("Compound may not be stoichiometrically balanced???")
				
			elif (coefficient_first_specie != 0.0) and (coefficient_second_specie == 0.0):		# Vertical line
				constant_deltamu_first_element = difference_enthalpy_adjusted / coefficient_first_specie
				if coefficient_first_specie > 0.0:
					vertical_left_values.append(constant_deltamu_first_element)
				elif coefficient_first_specie < 0.0:
					vertical_right_values.append(constant_deltamu_first_element)
				
				competing_compound_deltamu_first_element = np.ones(len(main_compound_deltamu_first_element)) * constant_deltamu_first_element
				competing_compound_deltamu_second_element = np.linspace(self.main_compound_enthalpy/self.main_compound_number_second_specie, 0, len(competing_compound_deltamu_first_element))
				
				main_compound_stability_limit_vertical = ( main_compound_enthalpy_adjusted - self.main_compound_number_first_specie*constant_deltamu_first_element ) / self.main_compound_number_second_specie
				competing_compound_deltamu_first_element_limit = [competing_compound_deltamu_first_element[i] for i in range(len(competing_compound_deltamu_first_element)) if (main_compound_stability_limit_vertical < competing_compound_deltamu_second_element[i])]
				competing_compound_deltamu_second_element_limit = [competing_compound_deltamu_second_element[i] for i in range(len(competing_compound_deltamu_first_element)) if (main_compound_stability_limit_vertical < competing_compound_deltamu_second_element[i])]
				
			elif (coefficient_first_specie == 0.0) and (coefficient_second_specie != 0.0):
				competing_compound_deltamu_first_element = copy.deepcopy(main_compound_deltamu_first_element)
				constant_deltamu_second_element = difference_enthalpy_adjusted / coefficient_second_specie
				competing_compound_deltamu_second_element = np.ones(len(competing_compound_deltamu_first_element)) * constant_deltamu_second_element
				
				if coefficient_second_specie > 0.0:
					stability_maximum_bound.append(competing_compound_deltamu_second_element)
				elif coefficient_second_specie < 0.0:
					stability_minimum_bound.append(competing_compound_deltamu_second_element)
				
				competing_compound_deltamu_first_element_limit = [competing_compound_deltamu_first_element[i] for i in range(len(competing_compound_deltamu_first_element)) if (main_compound_stability_limit[i] < competing_compound_deltamu_second_element[i])]
				competing_compound_deltamu_second_element_limit = [competing_compound_deltamu_second_element[i] for i in range(len(competing_compound_deltamu_first_element)) if (main_compound_stability_limit[i] < competing_compound_deltamu_second_element[i])]
			
			elif (coefficient_first_specie != 0.0) and (coefficient_second_specie != 0.0):
				competing_compound_deltamu_first_element = copy.deepcopy(main_compound_deltamu_first_element)
				competing_compound_deltamu_second_element = ( difference_enthalpy_adjusted - coefficient_first_specie*competing_compound_deltamu_first_element ) / coefficient_second_specie
				
				if coefficient_second_specie > 0.0:
					stability_maximum_bound.append(competing_compound_deltamu_second_element)
				elif coefficient_second_specie < 0.0:
					stability_minimum_bound.append(competing_compound_deltamu_second_element)
				
				competing_compound_deltamu_first_element_limit = [competing_compound_deltamu_first_element[i] for i in range(len(competing_compound_deltamu_first_element)) if (main_compound_stability_limit[i] < competing_compound_deltamu_second_element[i])]
				competing_compound_deltamu_second_element_limit = [competing_compound_deltamu_second_element[i] for i in range(len(competing_compound_deltamu_first_element)) if (main_compound_stability_limit[i] < competing_compound_deltamu_second_element[i])]
			
			try:
				# See if plot exists so it just needs to be updated
				self.competing_compound_plots[competing_compound].set_data(competing_compound_deltamu_first_element_limit, competing_compound_deltamu_second_element_limit)
			except:
				# If the plot doesn't exist initially, then create it
				self.competing_compound_plots[competing_compound], = self.quaternary_phase_diagram_plot_drawing.plot(competing_compound_deltamu_first_element_limit, competing_compound_deltamu_second_element_limit, label=self.Compound_Name_Formal(competing_compound))
				self.competing_compound_colorwheel[competing_compound] = self.competing_compound_plots[competing_compound].get_color()
		
		
		# Legend
		self.quaternary_phase_diagram_plot_drawing.legend(loc=3)
		
		stability_minimum_bound.append(main_compound_stability_limit)
		stability_maximum_bound.append(np.zeros(len(main_compound_deltamu_first_element)))
		stability_absolute_minimum = np.fromiter(map(max, zip(*itertools.chain(stability_minimum_bound))), dtype=np.float)
		stability_absolute_maximum = np.fromiter(map(min, zip(*itertools.chain(stability_maximum_bound))), dtype=np.float)
		
		self.main_compound_deltamu_first_element_cutoff = []
		self.stability_minimum_cutoff = []
		self.stability_maximum_cutoff = []
		
		for i in range(len(main_compound_deltamu_first_element)):
			if (vertical_left_values != []) and (vertical_right_values != []):
				if (stability_absolute_minimum[i] < stability_absolute_maximum[i]) and (main_compound_deltamu_first_element[i] < min(vertical_left_values)) and (main_compound_deltamu_first_element[i] > max(vertical_right_values)):
					self.main_compound_deltamu_first_element_cutoff.append(main_compound_deltamu_first_element[i])
					self.stability_minimum_cutoff.append(stability_absolute_minimum[i])
					self.stability_maximum_cutoff.append(stability_absolute_maximum[i])
			if (vertical_left_values != []) and (vertical_right_values == []):
				if (stability_absolute_minimum[i] < stability_absolute_maximum[i]) and (main_compound_deltamu_first_element[i] < min(vertical_left_values)):
					self.main_compound_deltamu_first_element_cutoff.append(main_compound_deltamu_first_element[i])
					self.stability_minimum_cutoff.append(stability_absolute_minimum[i])
					self.stability_maximum_cutoff.append(stability_absolute_maximum[i])
			if (vertical_left_values == []) and (vertical_right_values != []):
				if (stability_absolute_minimum[i] < stability_absolute_maximum[i]) and (main_compound_deltamu_first_element[i] > max(vertical_right_values)):
					self.main_compound_deltamu_first_element_cutoff.append(main_compound_deltamu_first_element[i])
					self.stability_minimum_cutoff.append(stability_absolute_minimum[i])
					self.stability_maximum_cutoff.append(stability_absolute_maximum[i])
			if (vertical_left_values == []) and (vertical_right_values == []):
				if (stability_absolute_minimum[i] < stability_absolute_maximum[i]):
					self.main_compound_deltamu_first_element_cutoff.append(main_compound_deltamu_first_element[i])
					self.stability_minimum_cutoff.append(stability_absolute_minimum[i])
					self.stability_maximum_cutoff.append(stability_absolute_maximum[i])
		
		try:
			self.phase_stability_region.remove()
			self.phase_stability_region = self.quaternary_phase_diagram_plot_drawing.fill_between(self.main_compound_deltamu_first_element_cutoff, self.stability_maximum_cutoff, self.stability_minimum_cutoff, facecolor='0.75')
		except:
			self.phase_stability_region = self.quaternary_phase_diagram_plot_drawing.fill_between(self.main_compound_deltamu_first_element_cutoff, self.stability_maximum_cutoff, self.stability_minimum_cutoff, facecolor='0.75')
		
		# Find the vertices of the phase stability region
		if self.phase_stability_region.get_paths() != []:	# Check if phase stability region exists
			self.PSR_vertices = self.Find_PhaseStabilityRegion_Vertices(self.phase_stability_region.get_paths())
			if self.PSR_vertices == []:
				try:
					self.PSR_vertices_plot.remove()
				except:
					pass
			elif self.PSR_vertices != []:
				try:
					self.PSR_vertices_plot.remove()
					self.PSR_vertices_plot = self.quaternary_phase_diagram_plot_drawing.scatter(*zip(*self.PSR_vertices), s=20, c='black')
				except:
					self.PSR_vertices_plot = self.quaternary_phase_diagram_plot_drawing.scatter(*zip(*self.PSR_vertices), s=20, c='black')
					pass
		elif self.phase_stability_region.get_paths() == []:	# If phase stability region does not exist, then done
			self.PSR_vertices = []
			try:
				self.PSR_vertices_plot.remove()
			except:
				pass
		
		# Draw the phase diagram
		self.quaternary_phase_diagram_plot_canvas.draw()
	
	
	def Find_PhaseStabilityRegion_Vertices(self, plots):
		
		# This function finds the vertices bounding the phase stability region. It takes
		#	the points of the phase stability region as input.
		
		PSR_vertices = []					# Keep track of the vertices of the phase stability region
		PSR_bound_slope_previous = None		# Keep track of the slope of a side of the phase stability region
		PSR_bounding_point_previous = None	# Keep track of the phase stability region points
		tolerance = 1E-6
		# Record a point as a vertex if the slope drastically changes
		for PSR_bounds in plots[0].iter_segments():
			PSR_bounding_point = PSR_bounds[0]
			try:
				PSR_bound_slope = (PSR_bounding_point[1]-PSR_bounding_point_previous[1]) / (PSR_bounding_point[0]-PSR_bounding_point_previous[0])
			except:
				PSR_bounding_point_previous = PSR_bounding_point
				PSR_bound_slope_previous = 0.0
				continue
			if (PSR_bound_slope < PSR_bound_slope_previous - tolerance) or (PSR_bound_slope > PSR_bound_slope_previous + tolerance):
				PSR_vertices.append(PSR_bounding_point_previous)
				PSR_vertices.append(PSR_bounding_point)
			PSR_bound_slope_previous = PSR_bound_slope
			PSR_bounding_point_previous = PSR_bounding_point
		# Some vertices may be repeated, so get rid of those
		PSR_vertices_Omit = []
		for PSR_vertices_Index in range(len(PSR_vertices)-1):
			if (np.linalg.norm(PSR_vertices[PSR_vertices_Index]-PSR_vertices[PSR_vertices_Index+1]) < 0.01):
				PSR_vertices_Omit.append(PSR_vertices[PSR_vertices_Index])
		PSR_vertices_unrepeated = [x for x in PSR_vertices if (not any((x is y for y in PSR_vertices_Omit)))]
		return PSR_vertices_unrepeated
	
	
	"""
	def Update_Annotation(self, ind):
		
		position = self.quaternary_phase_diagram_plot_figure.get_offsets()[ind["ind"][0]]
		self.PSR_vertices_annotation.xy = position
		text = "{}, {}".format("hi", "hello")
		self.PSR_vertices_annotation.set_text(text)
	
	
	def PSR_Hover(self, event):
		
		visibility = self.PSR_vertices_annotation.get_visible()
		
		if event.inaxes == self.quaternary_phase_diagram_plot_drawing:
			cont, ind = self.quaternary_phase_diagram_plot_figure.contains(event)
			if cont:
				self.Update_Annotation(ind)
				self.PSR_vertices_annotation.set_visible(True)
				self.quaternary_phase_diagram_plot_canvas.draw()
			else:
				if visibility:
					self.PSR_vertices_annotation.set_visible(False)
					self.quaternary_phase_diagram_plot_canvas.draw()
	"""












































































class ChemicalPotential_Quaternary_PhaseDiagramProjected2D_TripleView(QWidget):
	
	def __init__(self, parent = None, main_compound = None, first_element = None, second_element = None, third_element = None, fourth_element = None):
		
		# All elements in the periodic table
		self.all_elements = []
		for element in periodictable.elements:
			self.all_elements.append(str(element))
		
		# Font description for phase stability diagram plot
		self.font = {'family': 'sans-serif',
				'color':  'black',
				'weight': 'normal',
				'size': 12 }
		
		# Establish the first, second, third, and fourth species of the quaternary compound.
		# Note that this list is subject to change, depending on what the user chooses.
		self.main_compound  = main_compound
		self.first_element	= first_element
		self.second_element	= second_element
		self.third_element	= third_element
		self.fourth_element	= fourth_element
		self.elements_list  = [self.first_element, self.second_element, self.third_element, self.fourth_element]	# List of elements that gets updated as user selects element order
		self.elements_list_original = [self.first_element, self.second_element, self.third_element, self.fourth_element]	# Unchanged list of elements for compound naming purposes (e.g. see Compound_Name_Formal)
		
		# Store all extracted DFT data
		self.compounds_info = {}
		
		# Necessary numerical variables
		self.main_compound_number_first_specie  = 0	# Number of each specie in the main quaternary compound
		self.main_compound_number_second_specie = 0
		self.main_compound_number_third_specie  = 0
		self.main_compound_number_fourth_specie = 0
		self.main_compound_enthalpy = 0.0			# Enthalpy of the main quaternary compound
		self.phasediagram_endpoints = 0.0			# Endpoints for quaternary phase diagram
		self.mu4 = 0.0								# Track the fourth species mu value as the user changes it
		
		# Quaternary phase diagram object
		self.quaternary_phase_diagram_plot_figure = plt.figure(figsize=(7,7))
		#self.quaternary_phase_diagram_plot_figure.subplots_adjust(left=0.025, right=0.7)
		
		"""
		gs = self.quaternary_phase_diagram_plot_figure.add_gridspec(4, 4, wspace=0.0, hspace=0.0)
		self.quaternary_phase_diagram_plot_drawing12 = self.quaternary_phase_diagram_plot_figure.add_subplot(gs[0:1, 0:1])
		self.quaternary_phase_diagram_plot_drawing13 = self.quaternary_phase_diagram_plot_figure.add_subplot(gs[0:1, 2:3])
		self.quaternary_phase_diagram_plot_drawing23 = self.quaternary_phase_diagram_plot_figure.add_subplot(gs[2:3, 1:2])
		"""
		self.quaternary_phase_diagram_plot_drawing12 = self.quaternary_phase_diagram_plot_figure.add_subplot(221)
		self.quaternary_phase_diagram_plot_drawing13 = self.quaternary_phase_diagram_plot_figure.add_subplot(222)
		self.quaternary_phase_diagram_plot_drawing23 = self.quaternary_phase_diagram_plot_figure.add_subplot(223)
		
		self.quaternary_phase_diagram_plot_canvas = FigureCanvas(self.quaternary_phase_diagram_plot_figure)
		
		
		
		# Necessary plot-related variables
		self.main_compound_plot12 = None				# Main compound plot (holds plot object)
		self.main_compound_plot13 = None				# Main compound plot (holds plot object)
		self.main_compound_plot23 = None				# Main compound plot (holds plot object)
		self.competing_compound_plots12 = {}			# Competing compound plots (holds plot objects)
		self.competing_compound_plots13 = {}			# Competing compound plots (holds plot objects)
		self.competing_compound_plots23 = {}			# Competing compound plots (holds plot objects)
		self.competing_compounds_colorwheel = {}		# Legend
		"""
		self.colorscheme = plt.get_cmap("Set3")		# Competing compounds colorwheel
		"""
		
		"""
		self.PSR_vertices_plot = None				# Dots for intersections
		self.PSR_vertices = []						# Coordinates of each intersection
		"""
	
	
	
	###############################################################################################
	########################## Rewrite Compound Name Latex-Style ##################################
	###############################################################################################
	
	def Compound_Name_Formal(self, compound_name):
		
		# Go into the compounds_info dictionary and obtains the chemistry and stoichiometry of the compound of choice
		compound_species_info = self.compounds_info[compound_name]
		
		compound_name_formal = ""
		for species in compound_species_info["elements_list"]:				# Loop through the list of possible species that can be contained in the compound
			if compound_species_info[species] == 0:
				continue								# Don't add the species to the name if the compound doesn't contain the species
			elif compound_species_info[species] == 1:
				compound_name_formal += species			# Add the species to the name of the compound
			elif compound_species_info[species] > 1:
				compound_name_formal += species+"$_{"+str(int(compound_species_info[species]))+"}$"	# Add the species to the name of the compound
				#compound_name_formal += species+"<sub>"+str(int(compound_species_info[species]))+"</sub>"	# Add the species to the name of the compound
																											#	with a subscript for the stoichiometry
		
		return compound_name_formal
	
	
	
	
	###############################################################################################
	################################## Set the dependent element ##################################
	###############################################################################################
	
	def Set_Elements(self, first_element, second_element, third_element, fourth_element):
		
		self.first_element = first_element
		self.second_element = second_element
		self.third_element = third_element
		self.fourth_element = fourth_element
		self.elements_list = [self.first_element, self.second_element, self.third_element, self.fourth_element]
	
	
	
	
	
	def Update_PhaseDiagram_Object(self):
		
		# Number of elements in main compound
		self.main_compound_number_first_specie = self.compounds_info[self.main_compound][self.first_element]	# Number of first species in quaternary compound
		self.main_compound_number_second_specie = self.compounds_info[self.main_compound][self.second_element]	# Number of second species in quaternary compound
		self.main_compound_number_third_specie = self.compounds_info[self.main_compound][self.third_element]	# Number of third species in quaternary compound
		self.main_compound_number_fourth_specie = self.compounds_info[self.main_compound][self.fourth_element]	# Number of fourth species in quaternary compound
		
		# Enthalpy of quaternary compound
		self.main_compound_enthalpy = self.compounds_info[self.main_compound]["enthalpy"]	# Enthalpy of quaternary compound
		
		# Endpoints of phase diagram
		self.phasediagram_endpoints = min(self.main_compound_enthalpy/self.main_compound_number_first_specie, self.main_compound_enthalpy/self.main_compound_number_second_specie, self.main_compound_enthalpy/self.main_compound_number_third_specie, self.main_compound_enthalpy/self.main_compound_number_fourth_specie)
	
	
	def Update_PhaseDiagram_Plot_Axes(self):
		
		self.quaternary_phase_diagram_plot_drawing12.set_xlim(self.phasediagram_endpoints, 0.0)
		self.quaternary_phase_diagram_plot_drawing12.set_ylim(self.phasediagram_endpoints, 0.0)
		self.quaternary_phase_diagram_plot_drawing12.set_xlabel("$\Delta\mu_{"+self.first_element+"}$ (eV)",fontdict=self.font)
		self.quaternary_phase_diagram_plot_drawing12.set_ylabel("$\Delta\mu_{"+self.second_element+"}$ (eV)",fontdict=self.font,rotation=270,labelpad=20)
		self.quaternary_phase_diagram_plot_drawing12.xaxis.tick_top()
		self.quaternary_phase_diagram_plot_drawing12.yaxis.tick_right()
		self.quaternary_phase_diagram_plot_drawing12.xaxis.set_label_position("top")
		self.quaternary_phase_diagram_plot_drawing12.yaxis.set_label_position("right")
		self.quaternary_phase_diagram_plot_drawing12.spines['left'].set_visible(False)
		self.quaternary_phase_diagram_plot_drawing12.spines['bottom'].set_visible(False)
		self.quaternary_phase_diagram_plot_drawing12.set_aspect("equal")
		
		self.quaternary_phase_diagram_plot_drawing13.set_xlim(self.phasediagram_endpoints, 0.0)
		self.quaternary_phase_diagram_plot_drawing13.set_ylim(self.phasediagram_endpoints, 0.0)
		self.quaternary_phase_diagram_plot_drawing13.set_xlabel("$\Delta\mu_{"+self.first_element+"}$ (eV)",fontdict=self.font)
		self.quaternary_phase_diagram_plot_drawing13.set_ylabel("$\Delta\mu_{"+self.third_element+"}$ (eV)",fontdict=self.font,rotation=270,labelpad=20)
		self.quaternary_phase_diagram_plot_drawing13.xaxis.tick_top()
		self.quaternary_phase_diagram_plot_drawing13.yaxis.tick_right()
		self.quaternary_phase_diagram_plot_drawing13.xaxis.set_label_position("top")
		self.quaternary_phase_diagram_plot_drawing13.yaxis.set_label_position("right")
		self.quaternary_phase_diagram_plot_drawing13.spines['left'].set_visible(False)
		self.quaternary_phase_diagram_plot_drawing13.spines['bottom'].set_visible(False)
		self.quaternary_phase_diagram_plot_drawing13.set_aspect("equal")
		
		self.quaternary_phase_diagram_plot_drawing23.set_xlim(self.phasediagram_endpoints, 0.0)
		self.quaternary_phase_diagram_plot_drawing23.set_ylim(self.phasediagram_endpoints, 0.0)
		self.quaternary_phase_diagram_plot_drawing23.set_xlabel("$\Delta\mu_{"+self.second_element+"}$ (eV)",fontdict=self.font)
		self.quaternary_phase_diagram_plot_drawing23.set_ylabel("$\Delta\mu_{"+self.third_element+"}$ (eV)",fontdict=self.font,rotation=270,labelpad=20)
		self.quaternary_phase_diagram_plot_drawing23.xaxis.tick_top()
		self.quaternary_phase_diagram_plot_drawing23.yaxis.tick_right()
		self.quaternary_phase_diagram_plot_drawing23.xaxis.set_label_position("top")
		self.quaternary_phase_diagram_plot_drawing23.yaxis.set_label_position("right")
		self.quaternary_phase_diagram_plot_drawing23.spines['left'].set_visible(False)
		self.quaternary_phase_diagram_plot_drawing23.spines['bottom'].set_visible(False)
		self.quaternary_phase_diagram_plot_drawing23.set_aspect("equal")
		
		self.quaternary_phase_diagram_plot_figure.tight_layout()
	
	
	
	def Establish_CompetingCompounds_Colorwheel(self):
		
		color_counter = 0
		
		# Loop through all compounds in the database
		for competing_compound in self.compounds_info.keys():	
			
			# Skip if compound is either the main compound or one of the elements
			if (competing_compound in self.all_elements) or (competing_compound == self.main_compound):
				continue
			
			self.competing_compounds_colorwheel[competing_compound] = self.colorscheme(color_counter)
			color_counter += 1
	
	
	
	def Plot_PhaseDiagrams(self):
		
		# First-second elements
		main_compound_deltamu_first_element12, main_compound_stability_limit12, \
			competing_compounds_deltamu_first_element_limit12, competing_compounds_deltamu_second_element_limit12, \
			main_compound_deltamu_first_element_cutoff12, stability_minimum_cutoff12, stability_maximum_cutoff12 \
			= self.Calculate_PhaseDiagram(self.first_element, self.second_element, self.third_element)
		
		try:
			self.main_compound_plot12.set_data(main_compound_deltamu_first_element12, main_compound_stability_limit12)
		except:
			self.main_compound_plot12, = self.quaternary_phase_diagram_plot_drawing12.plot(main_compound_deltamu_first_element12, main_compound_stability_limit12, color='k')
		
		for competing_compound in self.compounds_info.keys():
			
			# Skip if compound is either the main compound or one of the elements
			if (competing_compound in self.all_elements) or (competing_compound == self.main_compound):
				continue
			
			try:
				# See if plot exists so it just needs to be updated
				self.competing_compound_plots12[competing_compound].set_data(competing_compounds_deltamu_first_element_limit12[competing_compound], competing_compounds_deltamu_second_element_limit12[competing_compound])
			except:
				# If the plot doesn't exist initially, then create it
				self.competing_compound_plots12[competing_compound], = self.quaternary_phase_diagram_plot_drawing12.plot(competing_compounds_deltamu_first_element_limit12[competing_compound], competing_compounds_deltamu_second_element_limit12[competing_compound], label=self.Compound_Name_Formal(competing_compound), color=self.competing_compounds_colorwheel[self.Compound_Name_Formal(competing_compound)])
		
		try:
			self.phase_stability_region12.remove()
			self.phase_stability_region12 = self.quaternary_phase_diagram_plot_drawing12.fill_between(main_compound_deltamu_first_element_cutoff12, stability_maximum_cutoff12, stability_minimum_cutoff12, facecolor='0.75')
		except:
			self.phase_stability_region12 = self.quaternary_phase_diagram_plot_drawing12.fill_between(main_compound_deltamu_first_element_cutoff12, stability_maximum_cutoff12, stability_minimum_cutoff12, facecolor='0.75')
		
		
		
		
		
		
		# First-third elements
		main_compound_deltamu_first_element13, main_compound_stability_limit13, \
			competing_compounds_deltamu_first_element_limit13, competing_compounds_deltamu_second_element_limit13, \
			main_compound_deltamu_first_element_cutoff13, stability_minimum_cutoff13, stability_maximum_cutoff13 \
			= self.Calculate_PhaseDiagram(self.first_element, self.third_element, self.second_element)
		
		try:
			self.main_compound_plot13.set_data(main_compound_deltamu_first_element13, main_compound_stability_limit13)
		except:
			self.main_compound_plot13, = self.quaternary_phase_diagram_plot_drawing13.plot(main_compound_deltamu_first_element13, main_compound_stability_limit13, color='k')
		
		for competing_compound in self.compounds_info.keys():
			
			# Skip if compound is either the main compound or one of the elements
			if (competing_compound in self.all_elements) or (competing_compound == self.main_compound):
				continue
			
			try:
				# See if plot exists so it just needs to be updated
				self.competing_compound_plots13[competing_compound].set_data(competing_compounds_deltamu_first_element_limit13[competing_compound], competing_compounds_deltamu_second_element_limit13[competing_compound])
			except:
				# If the plot doesn't exist initially, then create it
				self.competing_compound_plots13[competing_compound], = self.quaternary_phase_diagram_plot_drawing13.plot(competing_compounds_deltamu_first_element_limit13[competing_compound], competing_compounds_deltamu_second_element_limit13[competing_compound], label=self.Compound_Name_Formal(competing_compound), color=self.competing_compounds_colorwheel[self.Compound_Name_Formal(competing_compound)])
		
		try:
			self.phase_stability_region13.remove()
			self.phase_stability_region13 = self.quaternary_phase_diagram_plot_drawing13.fill_between(main_compound_deltamu_first_element_cutoff13, stability_maximum_cutoff13, stability_minimum_cutoff13, facecolor='0.75')
		except:
			self.phase_stability_region13 = self.quaternary_phase_diagram_plot_drawing13.fill_between(main_compound_deltamu_first_element_cutoff13, stability_maximum_cutoff13, stability_minimum_cutoff13, facecolor='0.75')
		
		
		
		
		
		
		
		
		# Second-third elements
		main_compound_deltamu_first_element23, main_compound_stability_limit23, \
			competing_compounds_deltamu_first_element_limit23, competing_compounds_deltamu_second_element_limit23, \
			main_compound_deltamu_first_element_cutoff23, stability_minimum_cutoff23, stability_maximum_cutoff23 \
			= self.Calculate_PhaseDiagram(self.second_element, self.third_element, self.first_element)
		
		try:
			self.main_compound_plot23.set_data(main_compound_deltamu_first_element23, main_compound_stability_limit23)
		except:
			self.main_compound_plot23, = self.quaternary_phase_diagram_plot_drawing23.plot(main_compound_deltamu_first_element23, main_compound_stability_limit23, color='k')
		
		for competing_compound in self.compounds_info.keys():
			
			# Skip if compound is either the main compound or one of the elements
			if (competing_compound in self.all_elements) or (competing_compound == self.main_compound):
				continue
			
			try:
				# See if plot exists so it just needs to be updated
				self.competing_compound_plots23[competing_compound].set_data(competing_compounds_deltamu_first_element_limit23[competing_compound], competing_compounds_deltamu_second_element_limit23[competing_compound])
			except:
				# If the plot doesn't exist initially, then create it
				self.competing_compound_plots23[competing_compound], = self.quaternary_phase_diagram_plot_drawing23.plot(competing_compounds_deltamu_first_element_limit23[competing_compound], competing_compounds_deltamu_second_element_limit23[competing_compound], label=self.Compound_Name_Formal(competing_compound), color=self.competing_compounds_colorwheel[self.Compound_Name_Formal(competing_compound)])
		
		try:
			self.phase_stability_region23.remove()
			self.phase_stability_region23 = self.quaternary_phase_diagram_plot_drawing23.fill_between(main_compound_deltamu_first_element_cutoff23, stability_maximum_cutoff23, stability_minimum_cutoff23, facecolor='0.75')
		except:
			self.phase_stability_region23 = self.quaternary_phase_diagram_plot_drawing23.fill_between(main_compound_deltamu_first_element_cutoff23, stability_maximum_cutoff23, stability_minimum_cutoff23, facecolor='0.75')
		
		
		"""
		# Plotting the phase diagram from DFT data is one of the main features of this app. This function
		#	single-handedly plots the phase diagram, so it's arguably one of the most important block of
		#	code in this script.
		# This function will be used in two forms:
		#	1) when the user generates the phase diagram by clicking the "Generate Plot" button
		#	2) when the user changes the mu4 value
		# Because there are different forms in which this function will be used, the programmer should take
		#	care when editing this block of code.
		
		main_compound_enthalpy_adjusted = self.main_compound_enthalpy - self.main_compound_number_fourth_specie*self.mu4	# Enthalpy adjusted for mu4 value (main compound)
		
		# First-second elements
		main_compound_deltamu_first_element12 = np.linspace(self.main_compound_enthalpy/self.main_compound_number_first_specie, 0, 1000)
		main_compound_stability_limit12 = (main_compound_enthalpy_adjusted - self.main_compound_number_first_specie*main_compound_deltamu_first_element12) / self.main_compound_number_second_specie
		try:
			self.main_compound_plot12.set_ydata(main_compound_stability_limit12)
		except:
			self.main_compound_plot12, = self.quaternary_phase_diagram_plot_drawing12.plot(main_compound_deltamu_first_element12, main_compound_stability_limit12)
		
		
		# First-third elements
		main_compound_deltamu_first_element13 = np.linspace(self.main_compound_enthalpy/self.main_compound_number_first_specie, 0, 1000)
		main_compound_stability_limit13 = (main_compound_enthalpy_adjusted - self.main_compound_number_first_specie*main_compound_deltamu_first_element13) / self.main_compound_number_third_specie
		try:
			self.main_compound_plot13.set_ydata(main_compound_stability_limit13)
		except:
			self.main_compound_plot13, = self.quaternary_phase_diagram_plot_drawing13.plot(main_compound_deltamu_first_element13, main_compound_stability_limit13)
		
		
		# Second-third elements
		main_compound_deltamu_second_element23 = np.linspace(self.main_compound_enthalpy/self.main_compound_number_second_specie, 0, 1000)
		main_compound_stability_limit23 = (main_compound_enthalpy_adjusted - self.main_compound_number_second_specie*main_compound_deltamu_second_element23) / self.main_compound_number_third_specie
		try:
			self.main_compound_plot23.set_ydata(main_compound_stability_limit23)
		except:
			self.main_compound_plot23, = self.quaternary_phase_diagram_plot_drawing23.plot(main_compound_deltamu_second_element23, main_compound_stability_limit23)
		
		
		# Find the bounds of the quaternary phase stability region (for shading the stability region of the quaternary compound)
		stability_minimum_bound12 = []
		stability_maximum_bound12 = []
		stability_minimum_bound13 = []
		stability_maximum_bound13 = []
		stability_minimum_bound23 = []
		stability_maximum_bound23 = []
		
		vertical_left_values12 = []
		vertical_right_values12 = []
		vertical_left_values13 = []
		vertical_right_values13 = []
		vertical_left_values23 = []
		vertical_right_values23 = []
		
		color_counter = 0
		
		# Loop through all compounds in the database
		for competing_compound in self.compounds_info.keys():
			
			# Skip if compound is either the main compound or one of the elements
			if (competing_compound in self.all_elements) or (competing_compound == self.main_compound):
				continue
			
			self.competing_compound_colorwheel[competing_compound] = self.colorscheme(color_counter)
			color_counter += 1
			
			
			# First-second elements
			try:
				competing_compound_number_first_specie12 = self.compounds_info[competing_compound][self.first_element]
			except:
				competing_compound_number_first_specie12 = 0.0
				pass
			try:
				competing_compound_number_second_specie12 = self.compounds_info[competing_compound][self.second_element]
			except:
				competing_compound_number_second_specie12 = 0.0
				pass
			try:
				competing_compound_number_third_specie12 = self.compounds_info[competing_compound][self.third_element]
			except:
				competing_compound_number_third_specie12 = 0.0
				pass
			try:
				competing_compound_number_fourth_specie12 = self.compounds_info[competing_compound][self.fourth_element]
			except:
				competing_compound_number_fourth_specie12 = 0.0
				pass
			
			competing_compound_enthalpy12 = self.compounds_info[competing_compound]["enthalpy"]
			competing_compound_enthalpy_adjusted12 = competing_compound_enthalpy12 - competing_compound_number_fourth_specie12*self.mu4		# Enthalpy adjusted for mu4 value (competing compound)
			
			difference_enthalpy_adjusted12 = competing_compound_enthalpy_adjusted12 - (competing_compound_number_third_specie12 / self.main_compound_number_third_specie) * main_compound_enthalpy_adjusted
			
			coefficient_first_specie12 = competing_compound_number_first_specie12 - (self.main_compound_number_first_specie * competing_compound_number_third_specie12) / self.main_compound_number_third_specie
			coefficient_second_specie12 = competing_compound_number_second_specie12 - (self.main_compound_number_second_specie * competing_compound_number_third_specie12) / self.main_compound_number_third_specie
			
			if (coefficient_first_specie12 == 0.0) and (coefficient_second_specie12 == 0.0):
				continue
				print("OMG SOMETHING'S WRONG!!!")
				print("Compound may not be stoichiometrically balanced???")
				
			elif (coefficient_first_specie12 != 0.0) and (coefficient_second_specie12 == 0.0):		# Vertical line
				constant_deltamu_first_element12 = difference_enthalpy_adjusted12 / coefficient_first_specie12
				if coefficient_first_specie12 > 0.0:
					vertical_left_values12.append(constant_deltamu_first_element12)
				elif coefficient_first_specie12 < 0.0:
					vertical_right_values12.append(constant_deltamu_first_element12)
				
				competing_compound_deltamu_first_element = np.ones(len(main_compound_deltamu_first_element)) * constant_deltamu_first_element
				competing_compound_deltamu_second_element = np.linspace(self.main_compound_enthalpy/self.main_compound_number_second_specie, 0, len(competing_compound_deltamu_first_element))
				
				main_compound_stability_limit_vertical = ( main_compound_enthalpy_adjusted - self.main_compound_number_first_specie*constant_deltamu_first_element ) / self.main_compound_number_second_specie
				competing_compound_deltamu_first_element_limit = [competing_compound_deltamu_first_element[i] for i in range(len(competing_compound_deltamu_first_element)) if (main_compound_stability_limit_vertical < competing_compound_deltamu_second_element[i])]
				competing_compound_deltamu_second_element_limit = [competing_compound_deltamu_second_element[i] for i in range(len(competing_compound_deltamu_first_element)) if (main_compound_stability_limit_vertical < competing_compound_deltamu_second_element[i])]
				
			elif (coefficient_first_specie == 0.0) and (coefficient_second_specie != 0.0):
				competing_compound_deltamu_first_element = copy.deepcopy(main_compound_deltamu_first_element)
				constant_deltamu_second_element = difference_enthalpy_adjusted / coefficient_second_specie
				competing_compound_deltamu_second_element = np.ones(len(competing_compound_deltamu_first_element)) * constant_deltamu_second_element
				
				if coefficient_second_specie > 0.0:
					stability_maximum_bound.append(competing_compound_deltamu_second_element)
				elif coefficient_second_specie < 0.0:
					stability_minimum_bound.append(competing_compound_deltamu_second_element)
				
				competing_compound_deltamu_first_element_limit = [competing_compound_deltamu_first_element[i] for i in range(len(competing_compound_deltamu_first_element)) if (main_compound_stability_limit[i] < competing_compound_deltamu_second_element[i])]
				competing_compound_deltamu_second_element_limit = [competing_compound_deltamu_second_element[i] for i in range(len(competing_compound_deltamu_first_element)) if (main_compound_stability_limit[i] < competing_compound_deltamu_second_element[i])]
			
			elif (coefficient_first_specie != 0.0) and (coefficient_second_specie != 0.0):
				competing_compound_deltamu_first_element = copy.deepcopy(main_compound_deltamu_first_element)
				competing_compound_deltamu_second_element = ( difference_enthalpy_adjusted - coefficient_first_specie*competing_compound_deltamu_first_element ) / coefficient_second_specie
				
				if coefficient_second_specie > 0.0:
					stability_maximum_bound.append(competing_compound_deltamu_second_element)
				elif coefficient_second_specie < 0.0:
					stability_minimum_bound.append(competing_compound_deltamu_second_element)
				
				competing_compound_deltamu_first_element_limit = [competing_compound_deltamu_first_element[i] for i in range(len(competing_compound_deltamu_first_element)) if (main_compound_stability_limit[i] < competing_compound_deltamu_second_element[i])]
				competing_compound_deltamu_second_element_limit = [competing_compound_deltamu_second_element[i] for i in range(len(competing_compound_deltamu_first_element)) if (main_compound_stability_limit[i] < competing_compound_deltamu_second_element[i])]
			
			try:
				# See if plot exists so it just needs to be updated
				self.competing_compound_plots12[competing_compound].set_data(competing_compound_deltamu_first_element_limit, competing_compound_deltamu_second_element_limit)
			except:
				# If the plot doesn't exist initially, then create it
				self.competing_compound_plots12[competing_compound], = self.quaternary_phase_diagram_plot_drawing12.plot(competing_compound_deltamu_first_element_limit, competing_compound_deltamu_second_element_limit, label=self.Compound_Name_Formal(competing_compound))
		
		
		# Legend
		#self.quaternary_phase_diagram_plot_drawing.legend(loc=3)
		
		stability_minimum_bound.append(main_compound_stability_limit)
		stability_maximum_bound.append(np.zeros(len(main_compound_deltamu_first_element)))
		stability_absolute_minimum = np.fromiter(map(max, zip(*itertools.chain(stability_minimum_bound))), dtype=np.float)
		stability_absolute_maximum = np.fromiter(map(min, zip(*itertools.chain(stability_maximum_bound))), dtype=np.float)
		
		self.main_compound_deltamu_first_element_cutoff = []
		self.stability_minimum_cutoff = []
		self.stability_maximum_cutoff = []
		
		for i in range(len(main_compound_deltamu_first_element)):
			if (vertical_left_values != []) and (vertical_right_values != []):
				if (stability_absolute_minimum[i] < stability_absolute_maximum[i]) and (main_compound_deltamu_first_element[i] < min(vertical_left_values)) and (main_compound_deltamu_first_element[i] > max(vertical_right_values)):
					self.main_compound_deltamu_first_element_cutoff.append(main_compound_deltamu_first_element[i])
					self.stability_minimum_cutoff.append(stability_absolute_minimum[i])
					self.stability_maximum_cutoff.append(stability_absolute_maximum[i])
			if (vertical_left_values != []) and (vertical_right_values == []):
				if (stability_absolute_minimum[i] < stability_absolute_maximum[i]) and (main_compound_deltamu_first_element[i] < min(vertical_left_values)):
					self.main_compound_deltamu_first_element_cutoff.append(main_compound_deltamu_first_element[i])
					self.stability_minimum_cutoff.append(stability_absolute_minimum[i])
					self.stability_maximum_cutoff.append(stability_absolute_maximum[i])
			if (vertical_left_values == []) and (vertical_right_values != []):
				if (stability_absolute_minimum[i] < stability_absolute_maximum[i]) and (main_compound_deltamu_first_element[i] > max(vertical_right_values)):
					self.main_compound_deltamu_first_element_cutoff.append(main_compound_deltamu_first_element[i])
					self.stability_minimum_cutoff.append(stability_absolute_minimum[i])
					self.stability_maximum_cutoff.append(stability_absolute_maximum[i])
			if (vertical_left_values == []) and (vertical_right_values == []):
				if (stability_absolute_minimum[i] < stability_absolute_maximum[i]):
					self.main_compound_deltamu_first_element_cutoff.append(main_compound_deltamu_first_element[i])
					self.stability_minimum_cutoff.append(stability_absolute_minimum[i])
					self.stability_maximum_cutoff.append(stability_absolute_maximum[i])
		
		try:
			self.phase_stability_region.remove()
			self.phase_stability_region = self.quaternary_phase_diagram_plot_drawing.fill_between(self.main_compound_deltamu_first_element_cutoff, self.stability_maximum_cutoff, self.stability_minimum_cutoff, facecolor='0.75')
		except:
			self.phase_stability_region = self.quaternary_phase_diagram_plot_drawing.fill_between(self.main_compound_deltamu_first_element_cutoff, self.stability_maximum_cutoff, self.stability_minimum_cutoff, facecolor='0.75')
		"""
		
		
		
		
		
		"""
		# Find the vertices of the phase stability region
		if self.phase_stability_region.get_paths() != []:	# Check if phase stability region exists
			self.PSR_vertices = self.Find_PhaseStabilityRegion_Vertices(self.phase_stability_region.get_paths())
			if self.PSR_vertices == []:
				try:
					self.PSR_vertices_plot.remove()
				except:
					pass
			elif self.PSR_vertices != []:
				try:
					self.PSR_vertices_plot.remove()
					self.PSR_vertices_plot = self.quaternary_phase_diagram_plot_drawing.scatter(*zip(*self.PSR_vertices), s=20, c='black')
				except:
					self.PSR_vertices_plot = self.quaternary_phase_diagram_plot_drawing.scatter(*zip(*self.PSR_vertices), s=20, c='black')
					pass
		elif self.phase_stability_region.get_paths() == []:	# If phase stability region does not exist, then done
			self.PSR_vertices = []
			try:
				self.PSR_vertices_plot.remove()
			except:
				pass
		"""
		
		
		
		
		
		
		
		# Draw the phase diagram
		self.quaternary_phase_diagram_plot_canvas.draw()
	
	
	
	
	"""
	def Find_PhaseStabilityRegion_Vertices(self, plots):
		
		# This function finds the vertices bounding the phase stability region. It takes
		#	the points of the phase stability region as input.
		
		PSR_vertices = []					# Keep track of the vertices of the phase stability region
		PSR_bound_slope_previous = None		# Keep track of the slope of a side of the phase stability region
		PSR_bounding_point_previous = None	# Keep track of the phase stability region points
		tolerance = 1E-6
		# Record a point as a vertex if the slope drastically changes
		for PSR_bounds in plots[0].iter_segments():
			PSR_bounding_point = PSR_bounds[0]
			try:
				PSR_bound_slope = (PSR_bounding_point[1]-PSR_bounding_point_previous[1]) / (PSR_bounding_point[0]-PSR_bounding_point_previous[0])
			except:
				PSR_bounding_point_previous = PSR_bounding_point
				PSR_bound_slope_previous = 0.0
				continue
			if (PSR_bound_slope < PSR_bound_slope_previous - tolerance) or (PSR_bound_slope > PSR_bound_slope_previous + tolerance):
				PSR_vertices.append(PSR_bounding_point_previous)
				PSR_vertices.append(PSR_bounding_point)
			PSR_bound_slope_previous = PSR_bound_slope
			PSR_bounding_point_previous = PSR_bounding_point
		# Some vertices may be repeated, so get rid of those
		PSR_vertices_Omit = []
		for PSR_vertices_Index in range(len(PSR_vertices)-1):
			if (np.linalg.norm(PSR_vertices[PSR_vertices_Index]-PSR_vertices[PSR_vertices_Index+1]) < 0.01):
				PSR_vertices_Omit.append(PSR_vertices[PSR_vertices_Index])
		PSR_vertices_unrepeated = [x for x in PSR_vertices if (not any((x is y for y in PSR_vertices_Omit)))]
		return PSR_vertices_unrepeated
	"""
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	def Calculate_PhaseDiagram(self, first_element, second_element, third_element):
		
		# Plotting the phase diagram from DFT data is one of the main features of this app. This function
		#	single-handedly plots the phase diagram, so it's arguably one of the most important block of
		#	code in this script.
		# This function will be used in two forms:
		#	1) when the user generates the phase diagram by clicking the "Generate Plot" button
		#	2) when the user changes the mu4 value
		# Because there are different forms in which this function will be used, the programmer should take
		#	care when editing this block of code.
		
		main_compound_enthalpy_adjusted = self.main_compound_enthalpy - self.main_compound_number_fourth_specie*self.mu4	# Enthalpy adjusted for mu4 value (main compound)
		
		main_compound_deltamu_first_element = np.linspace(self.main_compound_enthalpy/self.compounds_info[self.main_compound][first_element], 0, 1000)	# Array of possible mu1 values in the plot
		
		main_compound_stability_limit = (main_compound_enthalpy_adjusted - self.compounds_info[self.main_compound][first_element]*main_compound_deltamu_first_element) / self.compounds_info[self.main_compound][second_element]	# Array that represents the (diagonal) stability limit of the main compound
		
		
		"""
		# Check if the phase diagram has already been plotted
		try:
			self.main_compound_plot.set_ydata(main_compound_stability_limit)	# If the plot's already there, simply update the stability limit of the main compound
		except:
			self.main_compound_plot, = self.quaternary_phase_diagram_plot_drawing.plot(main_compound_deltamu_first_element, main_compound_stability_limit)	# If the main compound stabilitiy limit hasn't been plotted yet, plot it
		"""
		
		
		# Find the bounds of the quaternary phase stability region (for shading the stability region of the quaternary compound)
		stability_minimum_bound = []
		stability_maximum_bound = []
		
		vertical_left_values = []
		vertical_right_values = []
		
		competing_compounds_deltamu_first_element_limit = {}
		competing_compounds_deltamu_second_element_limit = {}
		
		# Loop through all compounds in the database
		for competing_compound in self.compounds_info.keys():
			
			# Skip if compound is either the main compound or one of the elements
			if (competing_compound in self.all_elements) or (competing_compound == self.main_compound):
				continue
			
			try:
				competing_compound_number_first_specie = self.compounds_info[competing_compound][first_element]
			except:
				competing_compound_number_first_specie = 0.0
				pass
			
			try:
				competing_compound_number_second_specie = self.compounds_info[competing_compound][second_element]
			except:
				competing_compound_number_second_specie = 0.0
				pass
			
			try:
				competing_compound_number_third_specie = self.compounds_info[competing_compound][third_element]
			except:
				competing_compound_number_third_specie = 0.0
				pass
			
			try:
				competing_compound_number_fourth_specie = self.compounds_info[competing_compound][self.fourth_element]
			except:
				competing_compound_number_fourth_specie = 0.0
				pass
			
			competing_compound_enthalpy = self.compounds_info[competing_compound]["enthalpy"]
			competing_compound_enthalpy_adjusted = competing_compound_enthalpy - competing_compound_number_fourth_specie*self.mu4		# Enthalpy adjusted for mu4 value (competing compound)
			
			difference_enthalpy_adjusted = competing_compound_enthalpy_adjusted - (competing_compound_number_third_specie / self.compounds_info[self.main_compound][third_element]) * main_compound_enthalpy_adjusted
			
			coefficient_first_specie = competing_compound_number_first_specie - (self.compounds_info[self.main_compound][first_element] * competing_compound_number_third_specie) / self.compounds_info[self.main_compound][third_element]
			coefficient_second_specie = competing_compound_number_second_specie - (self.compounds_info[self.main_compound][second_element] * competing_compound_number_third_specie) / self.compounds_info[self.main_compound][third_element]
			
			if (coefficient_first_specie == 0.0) and (coefficient_second_specie == 0.0):
				continue
				print("OMG SOMETHING'S WRONG!!!")
				print("Compound may not be stoichiometrically balanced???")
				
			elif (coefficient_first_specie != 0.0) and (coefficient_second_specie == 0.0):		# Vertical line
				constant_deltamu_first_element = difference_enthalpy_adjusted / coefficient_first_specie
				if coefficient_first_specie > 0.0:
					vertical_left_values.append(constant_deltamu_first_element)
				elif coefficient_first_specie < 0.0:
					vertical_right_values.append(constant_deltamu_first_element)
				
				competing_compound_deltamu_first_element = np.ones(len(main_compound_deltamu_first_element)) * constant_deltamu_first_element
				competing_compound_deltamu_second_element = np.linspace(self.main_compound_enthalpy/self.compounds_info[self.main_compound][second_element], 0, len(competing_compound_deltamu_first_element))
				
				main_compound_stability_limit_vertical = ( main_compound_enthalpy_adjusted - self.compounds_info[self.main_compound][first_element]*constant_deltamu_first_element ) / self.compounds_info[self.main_compound][second_element]
				competing_compound_deltamu_first_element_limit = [competing_compound_deltamu_first_element[i] for i in range(len(competing_compound_deltamu_first_element)) if (main_compound_stability_limit_vertical < competing_compound_deltamu_second_element[i])]
				competing_compound_deltamu_second_element_limit = [competing_compound_deltamu_second_element[i] for i in range(len(competing_compound_deltamu_first_element)) if (main_compound_stability_limit_vertical < competing_compound_deltamu_second_element[i])]
				
			elif (coefficient_first_specie == 0.0) and (coefficient_second_specie != 0.0):
				competing_compound_deltamu_first_element = copy.deepcopy(main_compound_deltamu_first_element)
				constant_deltamu_second_element = difference_enthalpy_adjusted / coefficient_second_specie
				competing_compound_deltamu_second_element = np.ones(len(competing_compound_deltamu_first_element)) * constant_deltamu_second_element
				
				if coefficient_second_specie > 0.0:
					stability_maximum_bound.append(competing_compound_deltamu_second_element)
				elif coefficient_second_specie < 0.0:
					stability_minimum_bound.append(competing_compound_deltamu_second_element)
				
				competing_compound_deltamu_first_element_limit = [competing_compound_deltamu_first_element[i] for i in range(len(competing_compound_deltamu_first_element)) if (main_compound_stability_limit[i] < competing_compound_deltamu_second_element[i])]
				competing_compound_deltamu_second_element_limit = [competing_compound_deltamu_second_element[i] for i in range(len(competing_compound_deltamu_first_element)) if (main_compound_stability_limit[i] < competing_compound_deltamu_second_element[i])]
			
			elif (coefficient_first_specie != 0.0) and (coefficient_second_specie != 0.0):
				competing_compound_deltamu_first_element = copy.deepcopy(main_compound_deltamu_first_element)
				competing_compound_deltamu_second_element = ( difference_enthalpy_adjusted - coefficient_first_specie*competing_compound_deltamu_first_element ) / coefficient_second_specie
				
				if coefficient_second_specie > 0.0:
					stability_maximum_bound.append(competing_compound_deltamu_second_element)
				elif coefficient_second_specie < 0.0:
					stability_minimum_bound.append(competing_compound_deltamu_second_element)
				
				competing_compound_deltamu_first_element_limit = [competing_compound_deltamu_first_element[i] for i in range(len(competing_compound_deltamu_first_element)) if (main_compound_stability_limit[i] < competing_compound_deltamu_second_element[i])]
				competing_compound_deltamu_second_element_limit = [competing_compound_deltamu_second_element[i] for i in range(len(competing_compound_deltamu_first_element)) if (main_compound_stability_limit[i] < competing_compound_deltamu_second_element[i])]
			
			competing_compounds_deltamu_first_element_limit[competing_compound] = competing_compound_deltamu_first_element_limit
			competing_compounds_deltamu_second_element_limit[competing_compound] = competing_compound_deltamu_second_element_limit
			
			"""
			try:
				# See if plot exists so it just needs to be updated
				self.competing_compound_plots[competing_compound].set_data(competing_compound_deltamu_first_element_limit, competing_compound_deltamu_second_element_limit)
			except:
				# If the plot doesn't exist initially, then create it
				self.competing_compound_plots[competing_compound], = self.quaternary_phase_diagram_plot_drawing.plot(competing_compound_deltamu_first_element_limit, competing_compound_deltamu_second_element_limit, label=self.Compound_Name_Formal(competing_compound))
			"""
		
		
		"""
		# Legend
		self.quaternary_phase_diagram_plot_drawing.legend(loc=3)
		"""
		
		stability_minimum_bound.append(main_compound_stability_limit)
		stability_maximum_bound.append(np.zeros(len(main_compound_deltamu_first_element)))
		stability_absolute_minimum = np.fromiter(map(max, zip(*itertools.chain(stability_minimum_bound))), dtype=np.float)
		stability_absolute_maximum = np.fromiter(map(min, zip(*itertools.chain(stability_maximum_bound))), dtype=np.float)
		
		main_compound_deltamu_first_element_cutoff = []
		stability_minimum_cutoff = []
		stability_maximum_cutoff = []
		
		for i in range(len(main_compound_deltamu_first_element)):
			if (vertical_left_values != []) and (vertical_right_values != []):
				if (stability_absolute_minimum[i] < stability_absolute_maximum[i]) and (main_compound_deltamu_first_element[i] < min(vertical_left_values)) and (main_compound_deltamu_first_element[i] > max(vertical_right_values)):
					main_compound_deltamu_first_element_cutoff.append(main_compound_deltamu_first_element[i])
					stability_minimum_cutoff.append(stability_absolute_minimum[i])
					stability_maximum_cutoff.append(stability_absolute_maximum[i])
			if (vertical_left_values != []) and (vertical_right_values == []):
				if (stability_absolute_minimum[i] < stability_absolute_maximum[i]) and (main_compound_deltamu_first_element[i] < min(vertical_left_values)):
					main_compound_deltamu_first_element_cutoff.append(main_compound_deltamu_first_element[i])
					stability_minimum_cutoff.append(stability_absolute_minimum[i])
					stability_maximum_cutoff.append(stability_absolute_maximum[i])
			if (vertical_left_values == []) and (vertical_right_values != []):
				if (stability_absolute_minimum[i] < stability_absolute_maximum[i]) and (main_compound_deltamu_first_element[i] > max(vertical_right_values)):
					main_compound_deltamu_first_element_cutoff.append(main_compound_deltamu_first_element[i])
					stability_minimum_cutoff.append(stability_absolute_minimum[i])
					stability_maximum_cutoff.append(stability_absolute_maximum[i])
			if (vertical_left_values == []) and (vertical_right_values == []):
				if (stability_absolute_minimum[i] < stability_absolute_maximum[i]):
					main_compound_deltamu_first_element_cutoff.append(main_compound_deltamu_first_element[i])
					stability_minimum_cutoff.append(stability_absolute_minimum[i])
					stability_maximum_cutoff.append(stability_absolute_maximum[i])
		
		
		
		return main_compound_deltamu_first_element, main_compound_stability_limit, competing_compounds_deltamu_first_element_limit, competing_compounds_deltamu_second_element_limit, main_compound_deltamu_first_element_cutoff, stability_minimum_cutoff, stability_maximum_cutoff
		
		"""
		try:
			self.phase_stability_region.remove()
			self.phase_stability_region = self.quaternary_phase_diagram_plot_drawing.fill_between(self.main_compound_deltamu_first_element_cutoff, self.stability_maximum_cutoff, self.stability_minimum_cutoff, facecolor='0.75')
		except:
			self.phase_stability_region = self.quaternary_phase_diagram_plot_drawing.fill_between(self.main_compound_deltamu_first_element_cutoff, self.stability_maximum_cutoff, self.stability_minimum_cutoff, facecolor='0.75')
		"""
		
		
		
		"""
		# Find the vertices of the phase stability region
		if self.phase_stability_region.get_paths() != []:	# Check if phase stability region exists
			self.PSR_vertices = self.Find_PhaseStabilityRegion_Vertices(self.phase_stability_region.get_paths())
			if self.PSR_vertices == []:
				try:
					self.PSR_vertices_plot.remove()
				except:
					pass
			elif self.PSR_vertices != []:
				try:
					self.PSR_vertices_plot.remove()
					self.PSR_vertices_plot = self.quaternary_phase_diagram_plot_drawing.scatter(*zip(*self.PSR_vertices), s=20, c='black')
				except:
					self.PSR_vertices_plot = self.quaternary_phase_diagram_plot_drawing.scatter(*zip(*self.PSR_vertices), s=20, c='black')
					pass
		elif self.phase_stability_region.get_paths() == []:	# If phase stability region does not exist, then done
			self.PSR_vertices = []
			try:
				self.PSR_vertices_plot.remove()
			except:
				pass
		
		# Draw the phase diagram
		self.quaternary_phase_diagram_plot_canvas.draw()
		"""
	
	
	
	
	
	
	
	"""
	def Find_PhaseStabilityRegion_Vertices(self, plots):
		
		# This function finds the vertices bounding the phase stability region. It takes
		#	the points of the phase stability region as input.
		
		PSR_vertices = []					# Keep track of the vertices of the phase stability region
		PSR_bound_slope_previous = None		# Keep track of the slope of a side of the phase stability region
		PSR_bounding_point_previous = None	# Keep track of the phase stability region points
		tolerance = 1E-6
		# Record a point as a vertex if the slope drastically changes
		for PSR_bounds in plots[0].iter_segments():
			PSR_bounding_point = PSR_bounds[0]
			try:
				PSR_bound_slope = (PSR_bounding_point[1]-PSR_bounding_point_previous[1]) / (PSR_bounding_point[0]-PSR_bounding_point_previous[0])
			except:
				PSR_bounding_point_previous = PSR_bounding_point
				PSR_bound_slope_previous = 0.0
				continue
			if (PSR_bound_slope < PSR_bound_slope_previous - tolerance) or (PSR_bound_slope > PSR_bound_slope_previous + tolerance):
				PSR_vertices.append(PSR_bounding_point_previous)
				PSR_vertices.append(PSR_bounding_point)
			PSR_bound_slope_previous = PSR_bound_slope
			PSR_bounding_point_previous = PSR_bounding_point
		# Some vertices may be repeated, so get rid of those
		PSR_vertices_Omit = []
		for PSR_vertices_Index in range(len(PSR_vertices)-1):
			if (np.linalg.norm(PSR_vertices[PSR_vertices_Index]-PSR_vertices[PSR_vertices_Index+1]) < 0.01):
				PSR_vertices_Omit.append(PSR_vertices[PSR_vertices_Index])
		PSR_vertices_unrepeated = [x for x in PSR_vertices if (not any((x is y for y in PSR_vertices_Omit)))]
		return PSR_vertices_unrepeated
	"""





