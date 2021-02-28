
__author__ = 'Michael_Lidia_Jiaxing_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'


import numpy as np
import matplotlib.path as mpltPath
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# Import compositional phase diagram object
from vtandem.visualization.plots.plot_composition_phasediagram import Plot_Composition_PhaseDiagram, PhaseRegion

from pymatgen.core.composition import Composition



class Plot_Composition_Quaternary_PhaseDiagram(Plot_Composition_PhaseDiagram):
	
	def __init__(self, main_compound = None, first_element = None, second_element = None, third_element = None, fourth_element = None, compounds_info = None, main_compound_info = None):
		
		# Establish the first, second, third, and fourth species of the quaternary compound.
		self.main_compound  = main_compound
		self.first_element	= first_element
		self.second_element	= second_element
		self.third_element	= third_element
		self.fourth_element	= fourth_element
		self.elements_list  = [self.first_element, self.second_element, self.third_element, self.fourth_element]
		
		
		
		# Calculate deltamu for third element
		try:
			main_compound_enthalpy = main_compound_info["dft_BulkEnergy"]
			for element in self.elements_list:
				main_compound_enthalpy -= main_compound_info["dft_"+element] * compounds_info[element]["mu0"]
			deltamu_fourth_element = main_compound_enthalpy / main_compound_info["dft_"+self.fourth_element]
		except:
			main_compound_enthalpy = compounds_info[self.main_compound]["dft_total_energy"]
			for element in self.elements_list:
				main_compound_enthalpy -= compounds_info[self.main_compound]["dft_"+element] * compounds_info[element]["mu0"]
			deltamu_fourth_element = main_compound_enthalpy / compounds_info[self.main_compound]["dft_"+self.fourth_element]
			pass
		
		# Record delta mu values
		self.deltamu_values = {}
		self.deltamu_values[first_element] = 0.0
		self.deltamu_values[second_element] = 0.0
		self.deltamu_values[third_element] = 0.0
		self.deltamu_values[fourth_element] = deltamu_fourth_element
		
		
		
		
		
		
		# Inherit all variables (plot object, etc.) from parent object (Composition_PhaseDiagram)
		super().__init__(type = "quaternary", main_compound_info = main_compound_info, compounds_info = compounds_info)
		
		
		
		
		
		
		
		"""
		# Generate compositional phase diagram
		self.Create_Compositional_PhaseDiagram()
		self.Plot_Compositional_PhaseDiagram()
		
		
		
		### Find all four-phase regions after compositional phase diagram is drawn.
		self.four_phase_region_objects = []
		"""
		# Find all four-phase regions in the quaternary composition space.
		self.Find_All_PhaseRegions()
		
		# Plot all the centroids as a scatter plot (MUST come after generating four-phase regions).
		self.Plot_Centroids()
		
		
		
		"""
		# Plot selected four-phase region.
		#	Store list of the four points constituting the four-phase region in self.fourphaseregion_selected.
		#	Store the polygon plot object in self.fourphaseregion_shade.
		self.fourphaseregion_selected = None
		self.fourphaseregion_shade = []
		"""
		self.composition_phasediagram_plot_figure.canvas.mpl_connect('button_press_event', self.Shade_FourPhaseRegion)
		
		"""
		self.composition_phasediagram_plot_figure.canvas.mpl_connect('button_press_event', self.Calculate_ChemicalPotentials_FourPhaseRegion)
		self.composition_phasediagram_plot_figure.canvas.mpl_connect('motion_notify_event', self.Hover)
		
		self.fourphaseregion_annotation = self.composition_phasediagram_plot_drawing.annotate("", xy=(0,0), xytext=(20,20), textcoords="offset points",
																								bbox = dict(boxstyle="round", fc="w"),
																								arrowprops = dict(arrowstyle = "->") )
		"""
	
	
	
	
	
	
	
	
	###############################################################################################
	########################## Plot Centroids of All Four-Phase Regions ###########################
	###############################################################################################
	
	def Plot_Centroids_Quaternary(self):
		
		# Check that all four-phase regions have been found
		if self.phase_region_objects == []:
			return
		
		# Set color of scatter plot points
		scatterplot_color = 'k'
		scatterplot_marker = '*'
		
		# Plot all centroids individually in scatter plot
		centroids = []
		for four_phase_region in self.phase_region_objects:
			centroid = four_phase_region.centroid
			centroids.append(centroid)
			four_phase_region.centroid_plot = self.composition_phasediagram_plot_drawing.scatter(	centroid[0],
																									centroid[1],
																									centroid[2],
																									color = scatterplot_color,
																									marker = scatterplot_marker	)
		
		# Plot all centroids together in scatter plot
		centroids = np.asarray(centroids)
		self.centroids_plot = self.composition_phasediagram_plot_drawing.scatter(	centroids[:,0],
																					centroids[:,1],
																					centroids[:,2],
																					color = scatterplot_color,
																					marker = scatterplot_marker )
	
	
	"""
	def Update_Annotation(self, index, event):
		
		# Get xy position of cursor on screen
		screen_position_xy = self.centroids_plot.get_offsets()[index]
		self.fourphaseregion_annotation.xy = screen_position_xy
		
		# Get name of four phase region
		for four_phase_region in self.four_phase_region_objects:
			if four_phase_region.centroid_plot.contains(event)[0]:
				text = four_phase_region.name
				break
		
		# Set annotation
		self.fourphaseregion_annotation.set_text(text)
		self.fourphaseregion_annotation.get_bbox_patch().set_facecolor("w")
		self.fourphaseregion_annotation.get_bbox_patch().set_alpha(1.0)
	
	
	
	def Hover(self, event):
		is_visible = self.fourphaseregion_annotation.get_visible()
		
		if event.inaxes == self.composition_phasediagram_plot_drawing:
			contains, index = self.centroids_plot.contains(event)
			if contains:
				self.Update_Annotation(index["ind"][0], event)
				self.fourphaseregion_annotation.set_visible(True)
				self.composition_phasediagram_plot_canvas.draw()
			else:
				if is_visible:
					self.fourphaseregion_annotation.set_visible(False)
					self.composition_phasediagram_plot_canvas.draw()
	"""
	
	
	
	###############################################################################################
	############# Shade Four-Phase Region in Compositional Phase Diagram when Clicked #############
	###############################################################################################
	
	def Shade_FourPhaseRegion(self, event):
		
		print("Four phase region selected...")
		
		# Check to see that all four-phase regions have been calculated
		if self.phase_region_objects == []:
			return
		
		# Read coordinates of clicked point
		contains, index = self.centroids_plot.contains(event)
		if not contains:
			# Remove selection
			if self.phaseregion_shade != []:
				for triangle_plot in self.phaseregion_shade:
					triangle_plot.remove()
			self.phaseregion_shade = []
			self.composition_phasediagram_plot_canvas.draw()
			return
		
		# Clear previous four-phase region
		if self.phaseregion_shade != []:
			for triangle_plot in self.phaseregion_shade:
				triangle_plot.remove()
		
		# Get selected four-phase region
		for four_phase_region in self.phase_region_objects:
			if four_phase_region.centroid_plot.contains(event)[0]:
				self.phaseregion_selected = four_phase_region
				break
		
		triangle1 = [self.phaseregion_selected.vertices[0], self.phaseregion_selected.vertices[1], self.phaseregion_selected.vertices[2]]
		triangle2 = [self.phaseregion_selected.vertices[0], self.phaseregion_selected.vertices[1], self.phaseregion_selected.vertices[3]]
		triangle3 = [self.phaseregion_selected.vertices[0], self.phaseregion_selected.vertices[2], self.phaseregion_selected.vertices[3]]
		triangle4 = [self.phaseregion_selected.vertices[1], self.phaseregion_selected.vertices[2], self.phaseregion_selected.vertices[3]]
		
		alpha = 0.75
		triangle_plot1 = Poly3DCollection(triangle1, fc='gray', ec='k', alpha=alpha)
		triangle_plot2 = Poly3DCollection(triangle2, fc='gray', ec='k', alpha=alpha)
		triangle_plot3 = Poly3DCollection(triangle3, fc='gray', ec='k', alpha=alpha)
		triangle_plot4 = Poly3DCollection(triangle4, fc='gray', ec='k', alpha=alpha)
		
		self.composition_phasediagram_plot_drawing.add_collection3d(triangle_plot1)
		self.composition_phasediagram_plot_drawing.add_collection3d(triangle_plot2)
		self.composition_phasediagram_plot_drawing.add_collection3d(triangle_plot3)
		self.composition_phasediagram_plot_drawing.add_collection3d(triangle_plot4)
		
		self.composition_phasediagram_plot_canvas.draw()
		
		self.phaseregion_shade = [triangle_plot1, triangle_plot2, triangle_plot3, triangle_plot4]
	
	
	
	
	"""
	###############################################################################################
	################ Calculate Chemical Potentials of Elements in Four-Phase Region ###############
	###############################################################################################
	
	def Calculate_ChemicalPotentials_FourPhaseRegion(self, event):
		
		# Check to see that a four-phase region has been selected
		if self.fourphaseregion_selected is None:
			return
		
		# Get the names of compounds constituting the four-phase region
		fourphaseregion_compounds = self.fourphaseregion_selected.name.replace(" ", "").split(",")
		
		# Get the matrix encoding the stoichiometry of each compound constituting the four-phase region
		composition_matrix = np.zeros((4, 4))
		for compound_index, compound in enumerate(fourphaseregion_compounds):
			for element_index, element in enumerate(self.elements_list):
				try:
					composition_matrix[compound_index][element_index] = self.compounds_info[compound][element]
				except:
					composition_matrix[compound_index][element_index] = 0.0
		
		# Get the enthalpies of formation of each compound in the four-phase region
		enthalpies_array = np.zeros((4, 1))
		for compound_index, compound in enumerate(fourphaseregion_compounds):
			competing_compound_enthalpy = self.compounds_info[compound]["dft_total_energy"] / self.compounds_info[compound]["formula_units"]
			for element in self.elements_list:
				if element not in self.compounds_info[compound].keys():
					continue
				competing_compound_enthalpy -= self.compounds_info[compound][element] * self.compounds_info[element]["mu0"]
			enthalpies_array[compound_index] = competing_compound_enthalpy
		
		# Solve for the delta mu values
		deltamus = np.dot( np.linalg.inv(composition_matrix), enthalpies_array )
		
		# Record delta mu values
		for element, deltamu in zip(self.elements_list, deltamus):
			self.deltamu_values[element] = deltamu[0]
	"""




"""
class FourPhaseRegion(object):
	def __init__(self):
		self.name = None
		self.vertices = None
		self.centroid = None
		self.centroid_plot = None
"""




