
__author__ = 'Michael_Lidia_Jiaxing_Benita_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'


import numpy as np
import itertools
import copy
import periodictable
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure


class Ternary_Phase_Diagram(object):
	
	def __init__(self, parent = None, main_compound = None, first_species = None, second_species = None, third_species = None):
		
		# Font description for phase stability diagram plot
		self.font = {'family': 'Arial',
				'color':  'black',
				'weight': 'normal',
				'size': 16 }
		
		# Label the first, second, and third species of the atoms in the ternary compound
		self.main_compound  = main_compound
		self.first_species  = first_species
		self.second_species = second_species
		self.third_species  = third_species
		self.species_list   = [self.first_species, self.second_species, self.third_species]					# Species list (order MAY change)
		
		#Store all extracted DFT data
		self.compounds_info = {}
		
		# Information about main ternary compound
		self.main_compound_number_first_specie  = 0	# Number of each specie in ternary compound
		self.main_compound_number_second_specie = 0
		self.main_compound_number_third_specie  = 0
		self.main_compound_enthalpy = 0.0			# Enthalpy of ternary compound
		self.phasediagram_endpoints = 0.0			# Endpoints for ternary phase diagram
		
		
		
		
		# Initialize necessary objects in app
		self.main_compound_plot = None
		self.competing_compound_plots = {}
		self.competing_compound_colorwheel = {}
		self.PSR_vertices_plot = None
		self.PSR_vertices = []
		
		
		
		# Ternary phase diagram plot
		self.ternary_phase_diagram_plot_figure = plt.figure()
		self.ternary_phase_diagram_plot_figure.subplots_adjust(left=0.05, right=0.8)
		self.ternary_phase_diagram_plot_drawing = self.ternary_phase_diagram_plot_figure.add_subplot(111)
		self.ternary_phase_diagram_plot_canvas = FigureCanvas(self.ternary_phase_diagram_plot_figure)
	
	
	
	def Activate_PhaseDiagram_Object(self):
		
		# Initialize the necessary numerical variables. This function is ideally called once
		#	the user has created self.compounds_info.
		
		self.main_compound_number_first_specie = self.compounds_info[self.main_compound][self.first_species]	# Number of first species in ternary compound
		self.main_compound_number_second_specie = self.compounds_info[self.main_compound][self.second_species]	# Number of second species in ternary compound
		self.main_compound_number_third_specie = self.compounds_info[self.main_compound][self.third_species]	# Number of third species in ternary compound
		self.main_compound_enthalpy = self.compounds_info[self.main_compound]["enthalpy"]	# Enthalpy of ternary compound
		self.phasediagram_endpoints = min(self.main_compound_enthalpy/self.main_compound_number_first_specie, self.main_compound_enthalpy/self.main_compound_number_second_specie, self.main_compound_enthalpy/self.main_compound_number_third_specie) # Endpoints for ternary phase diagram
	
	
	###############################################################################################
	##################################### Phase Diagram ###########################################
	###############################################################################################
	
	def Activate_PhaseDiagram_Plot_Axes(self):
		
		# This function simply creates the axes of the ternary phase diagram.
		
		self.ternary_phase_diagram_plot_drawing.set_xlim(self.phasediagram_endpoints, 0.0)
		self.ternary_phase_diagram_plot_drawing.set_ylim(self.phasediagram_endpoints, 0.0)
		self.ternary_phase_diagram_plot_drawing.set_xlabel("$\Delta\mu_{a}$ (eV)",fontdict=self.font,labelpad=12)
		self.ternary_phase_diagram_plot_drawing.set_ylabel("$\Delta\mu_{b}$ (eV)",fontdict=self.font,rotation=270,labelpad=20)
		self.ternary_phase_diagram_plot_drawing.xaxis.tick_top()
		self.ternary_phase_diagram_plot_drawing.yaxis.tick_right()
		self.ternary_phase_diagram_plot_drawing.tick_params(axis='both', labelsize=6)
		self.ternary_phase_diagram_plot_drawing.xaxis.set_label_position("top")
		self.ternary_phase_diagram_plot_drawing.yaxis.set_label_position("right")
		self.ternary_phase_diagram_plot_drawing.spines['left'].set_visible(False)
		self.ternary_phase_diagram_plot_drawing.spines['bottom'].set_visible(False)
		self.ternary_phase_diagram_plot_drawing.set_aspect("equal")
	
	
	def Update_PhaseDiagram_Plot_Axes(self):
		
		# This function resets the axes labels of the ternary phase diagram. This will be used when the 
		#	user changes the species to be plotted, since the phase stability diagram (the triangular
		#	shape) is not always the same between different species.
		
		self.ternary_phase_diagram_plot_drawing.set_xlabel("$\Delta\mu_{"+self.first_species+"}$ (eV)",fontdict=self.font,labelpad=12)
		self.ternary_phase_diagram_plot_drawing.set_ylabel("$\Delta\mu_{"+self.second_species+"}$ (eV)",fontdict=self.font,rotation=270,labelpad=20)
	
	
	
	
	def Plot_Ternary_PhaseDiagram(self):
		
		
		main_compound_deltamu_first_species = np.linspace(self.main_compound_enthalpy/self.main_compound_number_first_specie, 0, 1000)
		
		main_compound_stability_limit = (self.main_compound_enthalpy - self.main_compound_number_first_specie * main_compound_deltamu_first_species) / self.main_compound_number_second_specie
		
		
		try:
			self.main_compound_plot.set_ydata(main_compound_stability_limit)
		except:
			self.main_compound_plot, = self.ternary_phase_diagram_plot_drawing.plot(main_compound_deltamu_first_species, main_compound_stability_limit)
		
		
		
		stability_minimum_bound = []
		stability_maximum_bound = []
		
		vertical_left_values = []
		vertical_right_values = []
		
		
		for competing_compound_index, competing_compound in enumerate(self.compounds_info.keys()):
			
			foreign_compound = False
			for element in periodictable.elements:
				if str(element) in self.species_list:
					continue
				elif str(element) not in self.compounds_info[competing_compound].keys():
					continue
				else:
					if self.compounds_info[competing_compound][str(element)] != 0.0:
						foreign_compound = True
			
			
			if foreign_compound:
				continue
			elif (competing_compound in self.species_list) or (competing_compound == self.main_compound):
				continue
			else:
				#print "Analyzing the competing compound '"+competing_compound+"'..."
				
				competing_compound_number_first_specie = self.compounds_info[competing_compound][self.first_species]
				competing_compound_number_second_specie = self.compounds_info[competing_compound][self.second_species]
				competing_compound_number_third_specie = self.compounds_info[competing_compound][self.third_species]
				
				competing_compound_enthalpy = self.compounds_info[competing_compound]["enthalpy"]
				
				coefficient_first_specie = competing_compound_number_first_specie - float(competing_compound_number_third_specie)*float(self.main_compound_number_first_specie)/float(self.main_compound_number_third_specie)
				coefficient_second_specie = competing_compound_number_second_specie - float(competing_compound_number_third_specie)*float(self.main_compound_number_second_specie)/float(self.main_compound_number_third_specie)
				
				competing_compound_deltamu_first_species = copy.deepcopy(main_compound_deltamu_first_species)
				competing_compound_deltamu_second_species = ( competing_compound_enthalpy - float(competing_compound_number_third_specie)/float(self.main_compound_number_third_specie) * self.main_compound_enthalpy \
															- (coefficient_first_specie * main_compound_deltamu_first_species) ) / coefficient_second_specie
				
				if coefficient_second_specie > 0.0:
					stability_maximum_bound.append(competing_compound_deltamu_second_species)
				elif coefficient_second_specie < 0.0:
					stability_minimum_bound.append(competing_compound_deltamu_second_species)
				
				competing_compound_deltamu_first_species_limit = [competing_compound_deltamu_first_species[i] for i in range(len(main_compound_deltamu_first_species)) if (main_compound_stability_limit[i] < competing_compound_deltamu_second_species[i])]
				competing_compound_deltamu_second_species_limit = [competing_compound_deltamu_second_species[i] for i in range(len(main_compound_deltamu_first_species)) if (main_compound_stability_limit[i] < competing_compound_deltamu_second_species[i])]
				
				try:
					self.competing_compound_plots[competing_compound].set_data(competing_compound_deltamu_first_species_limit, competing_compound_deltamu_second_species_limit)
				except:
					# If the plot doesn't exist initially, then create it
					self.competing_compound_plots[competing_compound], = self.ternary_phase_diagram_plot_drawing.plot(competing_compound_deltamu_first_species_limit, competing_compound_deltamu_second_species_limit, label=competing_compound)
					self.competing_compound_colorwheel[competing_compound] = self.competing_compound_plots[competing_compound].get_color()
		
		self.ternary_phase_diagram_plot_drawing.legend(loc=3)
		
		stability_minimum_bound.append(main_compound_stability_limit)
		stability_maximum_bound.append(np.zeros(len(main_compound_deltamu_first_species)))
		stability_absolute_minimum = np.fromiter(map(max, zip(*itertools.chain(stability_minimum_bound))), dtype=np.float)
		stability_absolute_maximum = np.fromiter(map(min, zip(*itertools.chain(stability_maximum_bound))), dtype=np.float)
		
		self.main_compound_deltamu_first_species_cutoff = []
		self.stability_minimum_cutoff = []
		self.stability_maximum_cutoff = []
		
		for i in range(len(main_compound_deltamu_first_species)):
			if stability_absolute_minimum[i] < stability_absolute_maximum[i]:
				self.main_compound_deltamu_first_species_cutoff.append(main_compound_deltamu_first_species[i])
				self.stability_minimum_cutoff.append(stability_absolute_minimum[i])
				self.stability_maximum_cutoff.append(stability_absolute_maximum[i])
		
		try:
			self.phase_stability_region.remove()
			self.phase_stability_region = self.ternary_phase_diagram_plot_drawing.fill_between(self.main_compound_deltamu_first_species_cutoff, self.stability_maximum_cutoff, self.stability_minimum_cutoff, facecolor='0.75')
		except:
			self.phase_stability_region = self.ternary_phase_diagram_plot_drawing.fill_between(self.main_compound_deltamu_first_species_cutoff, self.stability_maximum_cutoff, self.stability_minimum_cutoff, facecolor='0.75')
		
		
		if self.phase_stability_region.get_paths() != []:
			PSR_Vertices = []
			PSR_Bound_Slope_Previous = None
			PSR_Bounding_Point_Previous = None
			tolerance = 1E-6
			for PSR_Bounds in self.phase_stability_region.get_paths()[0].iter_segments():
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
			self.PSR_vertices = PSR_Vertices_Unrepeated
			if PSR_Vertices_Unrepeated == []:
				try:
					self.PSR_vertices_plot.remove()
				except:
					pass
			elif PSR_Vertices_Unrepeated != []:
				try:
					self.PSR_vertices_plot.remove()
					self.PSR_vertices_plot = self.ternary_phase_diagram_plot_drawing.scatter(*zip(*PSR_Vertices_Unrepeated), s=20, c='black')
				except:
					self.PSR_vertices_plot = self.ternary_phase_diagram_plot_drawing.scatter(*zip(*PSR_Vertices_Unrepeated), s=20, c='black')
					pass
		elif self.phase_stability_region.get_paths() == []:
			self.PSR_vertices = []
			try:
				self.PSR_vertices_plot.remove()
			except:
				pass
		
		
		
		self.ternary_phase_diagram_plot_figure.canvas.draw_idle()


















