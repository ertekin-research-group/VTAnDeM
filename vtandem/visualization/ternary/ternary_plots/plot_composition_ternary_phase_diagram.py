
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'
__author__ = 'Michael_Lidia_Jiaxing_Elif'

import numpy as np
import matplotlib.path as mpltPath
import matplotlib.pyplot as plt

# Import compositional phase diagram object
from vtandem.visualization.plots.plot_composition_phasediagram import Plot_Composition_PhaseDiagram, PhaseRegion

from pymatgen.core.composition import Composition


class Plot_Composition_Ternary_PhaseDiagram(Plot_Composition_PhaseDiagram):

	def __init__(self, main_compound = None, first_element = None, second_element = None, third_element = None, compounds_info = None, main_compound_info = None):
		
		# Establish the first, second, and third species of the ternary compound.
		self.main_compound  = main_compound
		self.first_element	= first_element
		self.second_element	= second_element
		self.third_element	= third_element
		self.elements_list  = [self.first_element, self.second_element, self.third_element]
		
		"""
		# Record DFT data
		self.compounds_info = compounds_info
		"""
		
		
		
		
		# Calculate deltamu for third element
		try:
			main_compound_enthalpy = main_compound_info["dft_BulkEnergy"]
			for element in self.elements_list:
				main_compound_enthalpy -= main_compound_info["dft_"+element] * compounds_info[element]["mu0"]
			deltamu_third_element = main_compound_enthalpy / main_compound_info["dft_"+self.third_element]
		except:
			main_compound_enthalpy = compounds_info[self.main_compound]["dft_total_energy"]
			for element in self.elements_list:
				main_compound_enthalpy -= compounds_info[self.main_compound]["dft_"+element] * compounds_info[element]["mu0"]
			deltamu_third_element = main_compound_enthalpy / compounds_info[self.main_compound]["dft_"+self.third_element]
			pass
		
		
		# Keep track of mu values of the species in the ternary compound
		self.deltamu_values = {}
		self.deltamu_values[first_element] = 0.0
		self.deltamu_values[second_element] = 0.0
		self.deltamu_values[third_element] = deltamu_third_element
		
		
		
		
		
		# Inherit all variables (plot object, etc.) from parent object (Composition_PhaseDiagram)
		super().__init__(type = "ternary", main_compound_info = main_compound_info, compounds_info = compounds_info)
		
		
		
		
		
		
		
		
		"""
		# Generate compositional phase diagram
		self.Create_Compositional_PhaseDiagram()
		self.Plot_Compositional_PhaseDiagram()
		
		### Find all three-phase regions after compositional phase diagram is drawn
		self.threephaseregions = []
		# Store names of the three-phase regions (in the same order as self.threephaseregions).
		self.threephaseregion_names = []
		# Store the centroid of each three-phase region (in the same order as self.threephaseregions).
		self.threephaseregion_centroids = []
		"""
		# Find all three-phase regions in the quaternary composition space.
		self.Find_All_PhaseRegions()
		
		# Plot all the centroids as a scatter plot (MUST come after generating three-phase regions).
		self.Plot_Centroids()
		
		
		"""
		# Plot selected three-phase region.
		#	Store list of the three points constituting the three-phase region in self.threephaseregion_selected.
		#	Store the polygon plot object in self.threephaseregion_shade.
		self.threephaseregion_selected = None
		self.threephaseregion_shade = None
		"""
		self.composition_phasediagram_plot_figure.canvas.mpl_connect('button_press_event', self.Shade_ThreePhaseRegion)
		
		"""
		self.composition_phasediagram_plot_figure.canvas.mpl_connect('button_press_event', self.Calculate_ChemicalPotentials_ThreePhaseRegion)
		self.composition_phasediagram_plot_figure.canvas.mpl_connect('motion_notify_event', self.Hover)
		
		self.threephaseregion_annotation = self.composition_phasediagram_plot_drawing.annotate("", xy=(0,0), xytext=(20,20), textcoords="offset points",
																								bbox = dict(boxstyle="round", fc="w"),
																								arrowprops = dict(arrowstyle = "->") )
		"""
	
	
	
	
	
	
	###############################################################################################
	########################## Plot Centroids of All Three-Phase Regions ###########################
	###############################################################################################
	
	def Plot_Centroids_Ternary(self):
		
		# Check that all three-phase regions have been found
		if self.phase_region_objects == []:
			return
		
		# Set color of scatter plot points
		scatterplot_color = 'k'
		scatterplot_marker = '*'
		
		# Plot all centroids as scatter plot
		self.threephaseregion_centroids = np.asarray(self.threephaseregion_centroids)
		self.centroids_plot = self.composition_phasediagram_plot_drawing.scatter(	self.threephaseregion_centroids[:,0],
																					self.threephaseregion_centroids[:,1],
																					color = scatterplot_color,
																					marker = scatterplot_marker	)
	
	"""
	def Update_Annotation(self, index):
		
		position = self.centroids_plot.get_offsets()[index]
		self.threephaseregion_annotation.xy = position
		
		text = self.threephaseregion_names[index]
		
		self.threephaseregion_annotation.set_text(text)
		self.threephaseregion_annotation.get_bbox_patch().set_facecolor("w")
		self.threephaseregion_annotation.get_bbox_patch().set_alpha(1.0)
	
	
	def Hover(self, event):
		is_visible = self.threephaseregion_annotation.get_visible()
		
		if event.inaxes == self.composition_phasediagram_plot_drawing:
			contains, index = self.centroids_plot.contains(event)
			if contains:
				self.Update_Annotation(index["ind"][0])
				self.threephaseregion_annotation.set_visible(True)
				self.composition_phasediagram_plot_canvas.draw()
			else:
				if is_visible:
					self.threephaseregion_annotation.set_visible(False)
					self.composition_phasediagram_plot_canvas.draw()
	"""
	
	
	
	
	
	
	
	
	
	
	
	###############################################################################################
	############ Shade Three-Phase Region in Compositional Phase Diagram when Clicked #############
	###############################################################################################
	
	def Shade_ThreePhaseRegion(self, event):
		
		# Check to see that all three-phase regions have been calculated
		#if self.threephaseregions is None:
		if self.phase_region_objects == []:
			return
		
		# Read coordinates of clicked point
		point_x = event.xdata
		point_y = event.ydata
		
		# Check that the user presses somewhere on the plot (and not anywhere else)
		if (not isinstance(point_x, float)) or (not isinstance(point_y, float)):
			return
		
		# Shade in the selected three-phase region
		for three_phase_region in self.phase_region_objects:
			is_in_threephaseregion = mpltPath.Path(three_phase_region.vertices).contains_point([point_x, point_y])
			if is_in_threephaseregion:
				self.phaseregion_selected = three_phase_region
				if self.phaseregion_shade != []:
					self.phaseregion_shade.remove()
				self.phaseregion_shade = plt.Polygon(three_phase_region.vertices, color='gray', alpha=0.5)
				self.composition_phasediagram_plot_drawing.add_patch(self.phaseregion_shade)
				self.composition_phasediagram_plot_canvas.draw()
	
	
	"""
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
			entry = self.labels[tuple(point)]	# PDEntry object
			threephaseregion_compounds.append( entry.name )
		
		# Get the matrix encoding the stoichiometry of each compound constituting the three-phase region
		composition_matrix = np.zeros((3, 3))
		for compound_index, compound in enumerate(threephaseregion_compounds):
			for element_index, element in enumerate(self.elements_list):
				try:
					composition_matrix[compound_index][element_index] = self.compounds_info[compound]["dft_"+element]
				except:
					composition_matrix[compound_index][element_index] = 0.0
		
		# Get the enthalpies of formation of each compound in the three-phase region
		enthalpies_array = np.zeros((3, 1))
		for compound_index, compound in enumerate(threephaseregion_compounds):
			compound_enthalpy = self.compounds_info[compound]["dft_total_energy"]
			for element in self.compounds_info[compound]["elements_list"]:
				compound_enthalpy -= self.compounds_info[compound]["dft_"+element] * self.compounds_info[element]["mu0"]
			#enthalpies_array[compound_index] = self.compounds_info[compound]["enthalpy"]
			enthalpies_array[compound_index] = compound_enthalpy
		
		print(composition_matrix)
		print(enthalpies_array)
		
		# Solve for the delta mu values
		deltamus = np.dot( np.linalg.inv(composition_matrix), enthalpies_array )
		print(deltamus)
		
		# Record delta mu values
		for element, deltamu in zip(self.elements_list, deltamus):
			self.deltamu_values[element] = deltamu
	"""


"""
class ThreePhaseRegion(object):
	def __init__(self):
		self.name = None
		self.vertices = None
		self.centroid = None
		self.centroid_plot = None
"""




