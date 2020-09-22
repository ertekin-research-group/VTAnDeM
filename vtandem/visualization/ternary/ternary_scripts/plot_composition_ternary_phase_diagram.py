
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'
__author__ = 'Michael_Lidia_Jiaxing_Elif'

import numpy as np
import matplotlib.path as mpltPath
import matplotlib.pyplot as plt

# Import compositional phase diagram object
from vtandem.visualization.composition_phase_diagram import Composition_PhaseDiagram


class Composition_Ternary_PhaseDiagram(Composition_PhaseDiagram):

	def __init__(self, main_compound = None, first_element = None, second_element = None, third_element = None, compounds_info = None):
		
		# Inherit all variables (plot object, etc.) from parent object (Composition_PhaseDiagram)
		super().__init__(type = "ternary")
		
		# Establish the first, second, and third species of the ternary compound.
		self.main_compound  = main_compound
		self.first_element	= first_element
		self.second_element	= second_element
		self.third_element	= third_element
		self.elements_list  = [self.first_element, self.second_element, self.third_element]
		
		# Record DFT data
		self.compounds_info = compounds_info
		
		# Keep track of mu values of the species in the ternary compound
		self.deltamu_values = {}
		self.deltamu_values[first_element] = 0.0
		self.deltamu_values[second_element] = 0.0
		self.deltamu_values[third_element] = self.compounds_info[main_compound]["enthalpy"] / self.compounds_info[main_compound][self.third_element]
		
		# Generate compositional phase diagram
		self.Create_Compositional_PhaseDiagram()
		self.Plot_Compositional_PhaseDiagram()
		
		### Find all three-phase regions after compositional phase diagram is drawn.
		self.threephaseregions = []
		# Store names of the three-phase regions (in the same order as self.threephaseregions).
		self.threephaseregion_names = []
		# Store the centroid of each three-phase region (in the same order as self.threephaseregions).
		self.threephaseregion_centroids = []
		# Find all three-phase regions in the quaternary composition space.
		self.Find_All_ThreePhaseRegions()
		# Plot all the centroids as a scatter plot in the compositional phase diagram.
		self.Plot_Centroids()
		
		
		
		# Plot selected three-phase region.
		#	Store list of the three points constituting the three-phase region in self.threephaseregion_selected.
		#	Store the polygon plot object in self.threephaseregion_shade.
		self.threephaseregion_selected = None
		self.threephaseregion_shade = None
		self.composition_phasediagram_plot_figure.canvas.mpl_connect('button_press_event', self.Shade_ThreePhaseRegion)
		self.composition_phasediagram_plot_figure.canvas.mpl_connect('button_press_event', self.Calculate_ChemicalPotentials_ThreePhaseRegion)
		
		
		
		self.composition_phasediagram_plot_figure.canvas.mpl_connect('motion_notify_event', self.Hover)
		
		self.threephaseregion_annotation = self.composition_phasediagram_plot_drawing.annotate("", xy=(0,0), xytext=(20,20), textcoords="offset points",
																								bbox = dict(boxstyle="round", fc="w"),
																								arrowprops = dict(arrowstyle = "->") )
	
	
	
	
	
	
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
		for line in self.lines:
			endpoint1, endpoint2 = np.transpose(line)
			for point in self.labels.keys():
				point = np.asarray(point)
				if 		( any(np.array_equal(np.transpose(np.vstack((point, endpoint1))), i) for i in self.lines) or any(np.array_equal(np.transpose(np.vstack((endpoint1, point))), i) for i in self.lines) ) \
					and ( any(np.array_equal(np.transpose(np.vstack((point, endpoint2))), i) for i in self.lines) or any(np.array_equal(np.transpose(np.vstack((endpoint2, point))), i) for i in self.lines) ):
					
					# Check if three phase region has already been recorded
					if any( all( list(x) in threephaseregion for x in [point, endpoint1, endpoint2] ) for threephaseregion in self.threephaseregions):
						continue
					
					# Record three phase region
					three_phase_region = [list(point), list(endpoint1), list(endpoint2)]
					self.threephaseregions.append(three_phase_region)
					
					# Record names of compounds constituting the three-phase region
					three_phase_region_compound_names = []
					for pt in [point, endpoint1, endpoint2]:
						three_phase_region_compound_names.append( self.labels[tuple(pt)].name )
					self.threephaseregion_names.append(', '.join(three_phase_region_compound_names))
					
					# Find and record the centroid of the three-phase region
					centroid = ( point + endpoint1 + endpoint2 ) / 3.
					self.threephaseregion_centroids.append(centroid)
	
	
	
	
	
	###############################################################################################
	########################## Plot Centroids of All Three-Phase Regions ###########################
	###############################################################################################
	
	def Plot_Centroids(self):
		
		# Check that all three-phase regions have been found
		if self.threephaseregion_centroids == []:
			return
		
		# Plot all centroids as scatter plot
		self.threephaseregion_centroids = np.asarray(self.threephaseregion_centroids)
		self.centroids_plot = self.composition_phasediagram_plot_drawing.scatter(	self.threephaseregion_centroids[:,0],
																					self.threephaseregion_centroids[:,1]	)
	
	
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
				self.composition_phasediagram_plot_drawing.add_patch(self.threephaseregion_shade)
				self.composition_phasediagram_plot_canvas.draw()
	
	
	
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
					composition_matrix[compound_index][element_index] = self.compounds_info[compound][element]
				except:
					composition_matrix[compound_index][element_index] = 0.0
		
		# Get the enthalpies of formation of each compound in the three-phase region
		enthalpies_array = np.zeros((3, 1))
		for compound_index, compound in enumerate(threephaseregion_compounds):
			enthalpies_array[compound_index] = self.compounds_info[compound]["enthalpy"]
		
		# Solve for the delta mu values
		deltamus = np.dot( np.linalg.inv(composition_matrix), enthalpies_array )
		
		# Record delta mu values
		for element, deltamu in zip(self.elements_list, deltamus):
			self.deltamu_values[element] = deltamu








