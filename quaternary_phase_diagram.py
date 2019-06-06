
__author__ = 'Michael_Lidia_Jiaxing_Benita_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'


import numpy as np
import itertools
import copy
import periodictable
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure


class Quaternary_Phase_Diagram(object):
	
	def __init__(self, parent = None, main_compound = None, first_species = None, second_species = None, third_species = None, fourth_species = None):
		
		# Font description for phase stability diagram plot
		self.font = {'family': 'Arial',
				'color':  'black',
				'weight': 'normal',
				'size': 14 }
		
		# Establish the first, second, third, and fourth species of the quaternary compound.
		# Note that this list is subject to change, depending on what the user chooses.
		self.main_compound  = main_compound
		self.first_species	= first_species
		self.second_species	= second_species
		self.third_species	= third_species
		self.fourth_species	= fourth_species
		self.species_list   = [self.first_species, self.second_species, self.third_species, self.fourth_species]
		
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
		self.quaternary_phase_diagram_plot_figure.subplots_adjust(left=0.05, right=0.8)
		self.quaternary_phase_diagram_plot_drawing = self.quaternary_phase_diagram_plot_figure.add_subplot(111)
		self.quaternary_phase_diagram_plot_canvas = FigureCanvas(self.quaternary_phase_diagram_plot_figure)
		
		"""
		# Phase stability crossover points plot (separate figure)
		self.quaternary_phase_stability_crossovers_plot_figure = plt.figure()
		self.quaternary_phase_stability_crossovers_plot_drawing = self.quaternary_phase_stability_crossovers_plot_figure.add_subplot(111)
		self.quaternary_phase_stability_crossovers_plot_canvas = FigureCanvas(self.quaternary_phase_stability_crossovers_plot_figure)
		"""
		
		# Necessary plot-related variables
		self.main_compound_plot = None				# Main compound plot (holds plot object)
		self.competing_compound_plots = {}			# Competing compound plots (holds plot objects)
		self.competing_compound_colorwheel = {}		# Legend
		self.PSR_vertices_plot = None				# Dots for intersections
		self.PSR_vertices = []						# Coordinates of each intersection
		
		"""
		self.PSR_vertices_annotation = self.quaternary_phase_diagram_plot_drawing.annotate("", xy=(0,0), xytext=(0.05, 0.05), textcoords="offset points", \
																		bbox=dict(boxstyle="round", fc="w"), arrowprops=dict(arrowstyle="->"))
		self.PSR_vertices_annotation.set_visible(False)
		self.quaternary_phase_diagram_plot_figure.canvas.mpl_connect("motion_notify_event", self.PSR_Hover)
		"""
	
	
	
	def Activate_PhaseDiagram_Object(self):
		
		# Initialize the necessary numerical variables. This function is ideally called once
		#	the user has created self.compounds_info.
		
		self.main_compound_number_first_specie = self.compounds_info[self.main_compound][self.first_species]	# Number of first species in quaternary compound
		self.main_compound_number_second_specie = self.compounds_info[self.main_compound][self.second_species]	# Number of second species in quaternary compound
		self.main_compound_number_third_specie = self.compounds_info[self.main_compound][self.third_species]	# Number of third species in quaternary compound
		self.main_compound_number_fourth_specie = self.compounds_info[self.main_compound][self.fourth_species]	# Number of fourth species in quaternary compound
		self.main_compound_enthalpy = self.compounds_info[self.main_compound]["enthalpy"]	# Enthalpy of quaternary compound
		self.phasediagram_endpoints = min(self.main_compound_enthalpy/self.main_compound_number_first_specie, self.main_compound_enthalpy/self.main_compound_number_second_specie, self.main_compound_enthalpy/self.main_compound_number_third_specie, self.main_compound_enthalpy/self.main_compound_number_fourth_specie) # Endpoints for quaternary phase diagram
	
	
	def Activate_PhaseDiagram_Plot_Axes(self):
		
		# This function simply creates the axes of the quaternary phase diagram.
		
		self.quaternary_phase_diagram_plot_drawing.set_xlim(self.phasediagram_endpoints, 0.0)
		self.quaternary_phase_diagram_plot_drawing.set_ylim(self.phasediagram_endpoints, 0.0)
		self.quaternary_phase_diagram_plot_drawing.set_xlabel("$\Delta\mu_{a}$ (eV)",fontdict=self.font,labelpad=12)
		self.quaternary_phase_diagram_plot_drawing.set_ylabel("$\Delta\mu_{b}$ (eV)",fontdict=self.font,rotation=270,labelpad=20)
		self.quaternary_phase_diagram_plot_drawing.xaxis.tick_top()
		self.quaternary_phase_diagram_plot_drawing.yaxis.tick_right()
		self.quaternary_phase_diagram_plot_drawing.tick_params(axis='both', labelsize=6)
		self.quaternary_phase_diagram_plot_drawing.xaxis.set_label_position("top")
		self.quaternary_phase_diagram_plot_drawing.yaxis.set_label_position("right")
		self.quaternary_phase_diagram_plot_drawing.spines['left'].set_visible(False)
		self.quaternary_phase_diagram_plot_drawing.spines['bottom'].set_visible(False)
		self.quaternary_phase_diagram_plot_drawing.set_aspect("equal")
	
	
	def Update_PhaseDiagram_Plot_Axes(self):
		
		# This function changes the axes labels of the quaternary phase diagram when the user selects
		#	the species to be plotted. This will be used in conjunction with the 
		#	Activate_PhaseDiagram_Plot_Axes function, which creates the initial phase diagram plot.
		
		self.quaternary_phase_diagram_plot_drawing.set_xlabel("$\Delta\mu_{"+self.first_species+"}$ (eV)",fontdict=self.font,labelpad=12)
		self.quaternary_phase_diagram_plot_drawing.set_ylabel("$\Delta\mu_{"+self.second_species+"}$ (eV)",fontdict=self.font,rotation=270,labelpad=20)
	
	
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
		
		main_compound_deltamu_first_species = np.linspace(self.main_compound_enthalpy/self.main_compound_number_first_specie, 0, 1000)	# Array of possible mu1 values in the plot
		
		main_compound_stability_limit = (main_compound_enthalpy_adjusted - self.main_compound_number_first_specie*main_compound_deltamu_first_species) / self.main_compound_number_second_specie	# Array that represents the (diagonal) stability limit of the main compound
		
		# Check if the phase diagram has already been plotted
		try:
			self.main_compound_plot.set_ydata(main_compound_stability_limit)	# If the plot's already there, simply update the stability limit of the main compound
		except:
			self.main_compound_plot, = self.quaternary_phase_diagram_plot_drawing.plot(main_compound_deltamu_first_species, main_compound_stability_limit)	# If the main compound stabilitiy limit hasn't been plotted yet, plot it
		
		# Find the bounds of the quaternary phase stability region (for shading the stability region of the quaternary compound)
		stability_minimum_bound = []
		stability_maximum_bound = []
		
		vertical_left_values = []
		vertical_right_values = []
		
		# Loop through all compounds in the database
		for competing_compound_index, competing_compound in enumerate(self.compounds_info.keys()):

			# Check if the compound should be considered a "competing compound"
			foreign_compound = False
			for element in periodictable.elements:
				if str(element) in self.species_list:
					continue
				elif str(element) not in self.compounds_info[competing_compound].keys():
					continue
				else:
					if self.compounds_info[competing_compound][str(element)] != 0.0:
						foreign_compound = True
			
			if foreign_compound:	# Foreign compound (e.g. compound containing at least one element not in the quaternary compound)
				continue
			elif (competing_compound in self.species_list) or (competing_compound == self.main_compound):	# Single-element compounds and the main compound are not considered "competing compounds" of the main quaternary
				continue
			else:
				
				competing_compound_number_first_specie = self.compounds_info[competing_compound][self.first_species]
				competing_compound_number_second_specie = self.compounds_info[competing_compound][self.second_species]
				competing_compound_number_third_specie = self.compounds_info[competing_compound][self.third_species]
				competing_compound_number_fourth_specie = self.compounds_info[competing_compound][self.fourth_species]
				
				competing_compound_enthalpy = self.compounds_info[competing_compound]["enthalpy"]
				competing_compound_enthalpy_adjusted = competing_compound_enthalpy - competing_compound_number_fourth_specie*self.mu4		# Enthalpy adjusted for mu4 value (competing compound)
				
				difference_enthalpy_adjusted = competing_compound_enthalpy_adjusted - (competing_compound_number_third_specie / self.main_compound_number_third_specie) * main_compound_enthalpy_adjusted
				
				coefficient_first_specie = competing_compound_number_first_specie - (self.main_compound_number_first_specie * competing_compound_number_third_specie) / self.main_compound_number_third_specie
				coefficient_second_specie = competing_compound_number_second_specie - (self.main_compound_number_second_specie * competing_compound_number_third_specie) / self.main_compound_number_third_specie
				
				if (coefficient_first_specie == 0.0) and (coefficient_second_specie == 0.0):
					continue
					print("OMG SOMETHING'S WRONG!!!")
					print("Compound may not be stoichiomatrically balanced???")
					
				elif (coefficient_first_specie != 0.0) and (coefficient_second_specie == 0.0):		# Vertical line
					constant_deltamu_first_species = difference_enthalpy_adjusted / coefficient_first_specie
					if coefficient_first_specie > 0.0:
						vertical_left_values.append(constant_deltamu_first_species)
					elif coefficient_first_specie < 0.0:
						vertical_right_values.append(constant_deltamu_first_species)
					
					competing_compound_deltamu_first_species = np.ones(len(main_compound_deltamu_first_species)) * constant_deltamu_first_species
					competing_compound_deltamu_second_species = np.linspace(self.main_compound_enthalpy/self.main_compound_number_second_specie, 0, len(competing_compound_deltamu_first_species))
					
					main_compound_stability_limit_vertical = ( main_compound_enthalpy_adjusted - self.main_compound_number_first_specie*constant_deltamu_first_species ) / self.main_compound_number_second_specie
					competing_compound_deltamu_first_species_limit = [competing_compound_deltamu_first_species[i] for i in range(len(competing_compound_deltamu_first_species)) if (main_compound_stability_limit_vertical < competing_compound_deltamu_second_species[i])]
					competing_compound_deltamu_second_species_limit = [competing_compound_deltamu_second_species[i] for i in range(len(competing_compound_deltamu_first_species)) if (main_compound_stability_limit_vertical < competing_compound_deltamu_second_species[i])]
					
				elif (coefficient_first_specie == 0.0) and (coefficient_second_specie != 0.0):
					competing_compound_deltamu_first_species = copy.deepcopy(main_compound_deltamu_first_species)
					constant_deltamu_second_species = difference_enthalpy_adjusted / coefficient_second_specie
					competing_compound_deltamu_second_species = np.ones(len(competing_compound_deltamu_first_species)) * constant_deltamu_second_species
					
					if coefficient_second_specie > 0.0:
						stability_maximum_bound.append(competing_compound_deltamu_second_species)
					elif coefficient_second_specie < 0.0:
						stability_minimum_bound.append(competing_compound_deltamu_second_species)
					
					competing_compound_deltamu_first_species_limit = [competing_compound_deltamu_first_species[i] for i in range(len(competing_compound_deltamu_first_species)) if (main_compound_stability_limit[i] < competing_compound_deltamu_second_species[i])]
					competing_compound_deltamu_second_species_limit = [competing_compound_deltamu_second_species[i] for i in range(len(competing_compound_deltamu_first_species)) if (main_compound_stability_limit[i] < competing_compound_deltamu_second_species[i])]
				
				elif (coefficient_first_specie != 0.0) and (coefficient_second_specie != 0.0):
					competing_compound_deltamu_first_species = copy.deepcopy(main_compound_deltamu_first_species)
					competing_compound_deltamu_second_species = ( difference_enthalpy_adjusted - coefficient_first_specie*competing_compound_deltamu_first_species ) / coefficient_second_specie
					
					if coefficient_second_specie > 0.0:
						stability_maximum_bound.append(competing_compound_deltamu_second_species)
					elif coefficient_second_specie < 0.0:
						stability_minimum_bound.append(competing_compound_deltamu_second_species)
					
					competing_compound_deltamu_first_species_limit = [competing_compound_deltamu_first_species[i] for i in range(len(competing_compound_deltamu_first_species)) if (main_compound_stability_limit[i] < competing_compound_deltamu_second_species[i])]
					competing_compound_deltamu_second_species_limit = [competing_compound_deltamu_second_species[i] for i in range(len(competing_compound_deltamu_first_species)) if (main_compound_stability_limit[i] < competing_compound_deltamu_second_species[i])]
				
				try:
					# See if plot exists so it just needs to be updated
					self.competing_compound_plots[competing_compound].set_data(competing_compound_deltamu_first_species_limit, competing_compound_deltamu_second_species_limit)
				except:
					# If the plot doesn't exist initially, then create it
					self.competing_compound_plots[competing_compound], = self.quaternary_phase_diagram_plot_drawing.plot(competing_compound_deltamu_first_species_limit, competing_compound_deltamu_second_species_limit, label=competing_compound)
					self.competing_compound_colorwheel[competing_compound] = self.competing_compound_plots[competing_compound].get_color()
		
		# Legend
		self.quaternary_phase_diagram_plot_drawing.legend(loc=3)
		
		stability_minimum_bound.append(main_compound_stability_limit)
		stability_maximum_bound.append(np.zeros(len(main_compound_deltamu_first_species)))
		stability_absolute_minimum = np.fromiter(map(max, zip(*itertools.chain(stability_minimum_bound))), dtype=np.float)
		stability_absolute_maximum = np.fromiter(map(min, zip(*itertools.chain(stability_maximum_bound))), dtype=np.float)
		
		self.main_compound_deltamu_first_species_cutoff = []
		self.stability_minimum_cutoff = []
		self.stability_maximum_cutoff = []
		
		"""
		for i in range(len(main_compound_deltamu_first_species)):
			if vertical_left_values != []:
				if (stability_absolute_minimum[i] < stability_absolute_maximum[i]) and (main_compound_deltamu_first_species[i] < min(vertical_left_values)):
					self.main_compound_deltamu_first_species_cutoff.append(main_compound_deltamu_first_species[i])
					self.stability_minimum_cutoff.append(stability_absolute_minimum[i])
					self.stability_maximum_cutoff.append(stability_absolute_maximum[i])
			if vertical_right_values != []:
				if (stability_absolute_minimum[i] < stability_absolute_maximum[i]) and (main_compound_deltamu_first_species[i] > max(vertical_right_values)):
					self.main_compound_deltamu_first_species_cutoff.append(main_compound_deltamu_first_species[i])
					self.stability_minimum_cutoff.append(stability_absolute_minimum[i])
					self.stability_maximum_cutoff.append(stability_absolute_maximum[i])
		"""
		for i in range(len(main_compound_deltamu_first_species)):
			if (vertical_left_values != []) and (vertical_right_values != []):
				if (stability_absolute_minimum[i] < stability_absolute_maximum[i]) and (main_compound_deltamu_first_species[i] < min(vertical_left_values)) and (main_compound_deltamu_first_species[i] > max(vertical_right_values)):
					self.main_compound_deltamu_first_species_cutoff.append(main_compound_deltamu_first_species[i])
					self.stability_minimum_cutoff.append(stability_absolute_minimum[i])
					self.stability_maximum_cutoff.append(stability_absolute_maximum[i])
			if (vertical_left_values != []) and (vertical_right_values == []):
				if (stability_absolute_minimum[i] < stability_absolute_maximum[i]) and (main_compound_deltamu_first_species[i] < min(vertical_left_values)):
					self.main_compound_deltamu_first_species_cutoff.append(main_compound_deltamu_first_species[i])
					self.stability_minimum_cutoff.append(stability_absolute_minimum[i])
					self.stability_maximum_cutoff.append(stability_absolute_maximum[i])
			if (vertical_left_values == []) and (vertical_right_values != []):
				if (stability_absolute_minimum[i] < stability_absolute_maximum[i]) and (main_compound_deltamu_first_species[i] > max(vertical_right_values)):
					self.main_compound_deltamu_first_species_cutoff.append(main_compound_deltamu_first_species[i])
					self.stability_minimum_cutoff.append(stability_absolute_minimum[i])
					self.stability_maximum_cutoff.append(stability_absolute_maximum[i])
			if (vertical_left_values == []) and (vertical_right_values == []):
				if (stability_absolute_minimum[i] < stability_absolute_maximum[i]):
					self.main_compound_deltamu_first_species_cutoff.append(main_compound_deltamu_first_species[i])
					self.stability_minimum_cutoff.append(stability_absolute_minimum[i])
					self.stability_maximum_cutoff.append(stability_absolute_maximum[i])
		
		try:
			self.phase_stability_region.remove()
			self.phase_stability_region = self.quaternary_phase_diagram_plot_drawing.fill_between(self.main_compound_deltamu_first_species_cutoff, self.stability_maximum_cutoff, self.stability_minimum_cutoff, facecolor='0.75')
		except:
			self.phase_stability_region = self.quaternary_phase_diagram_plot_drawing.fill_between(self.main_compound_deltamu_first_species_cutoff, self.stability_maximum_cutoff, self.stability_minimum_cutoff, facecolor='0.75')
		
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
		
		# Draw the new phase diagram
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






