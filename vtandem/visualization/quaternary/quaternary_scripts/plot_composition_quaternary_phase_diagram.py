
__author__ = 'Michael_Lidia_Jiaxing_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'


import numpy as np
import matplotlib.path as mpltPath
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

# Import compositional phase diagram object
from vtandem.visualization.composition_phase_diagram import Composition_PhaseDiagram


class Composition_Quaternary_PhaseDiagram(Composition_PhaseDiagram):

	def __init__(self, main_compound = None, first_element = None, second_element = None, third_element = None, fourth_element = None, compounds_info = None):
		
		# Inherit all variables (plot object, etc.) from parent object (Composition_PhaseDiagram)
		super().__init__(type = "quaternary")
		
		# Establish the first, second, third, and fourth species of the quaternary compound.
		self.main_compound  = main_compound
		self.first_element	= first_element
		self.second_element	= second_element
		self.third_element	= third_element
		self.fourth_element	= fourth_element
		self.elements_list  = [self.first_element, self.second_element, self.third_element, self.fourth_element]
		
		# Record DFT data
		self.compounds_info = compounds_info
		
		# Record delta mu values
		self.deltamu_values = {}
		self.deltamu_values[first_element] = 0.0
		self.deltamu_values[second_element] = 0.0
		self.deltamu_values[third_element] = 0.0
		self.deltamu_values[fourth_element] = self.compounds_info[main_compound]["enthalpy"] / self.compounds_info[main_compound][self.fourth_element]
		
		
		# Generate compositional phase diagram
		self.Create_Compositional_PhaseDiagram()
		self.Plot_Compositional_PhaseDiagram()
		
		
		### Find all four-phase regions after compositional phase diagram is drawn.
		
		# Store lists of four points constituting the four-phase regions.
		self.fourphaseregions = []
		# Store names of the four-phase regions (in the same order as self.fourphaseregions).
		self.fourphaseregion_names = []
		# Store the centroid of each four-phase region (in the same order as self.fourphaseregions).
		self.fourphaseregion_centroids = []
		# Find all four-phase regions in the quaternary composition space.
		self.Find_All_FourPhaseRegions()
		# Plot all the centroids as a scatter plot in the compositional phase diagram.
		self.Plot_Centroids()
		
		# Plot selected four-phase region.
		#	Store list of the four points constituting the four-phase region in self.fourphaseregion_selected.
		#	Store the polygon plot object in self.fourphaseregion_shade.
		self.fourphaseregion_selected = None
		self.fourphaseregion_shade = []
		self.composition_phasediagram_plot_figure.canvas.mpl_connect('button_press_event', self.Shade_FourPhaseRegion)
		self.composition_phasediagram_plot_figure.canvas.mpl_connect('button_press_event', self.Calculate_ChemicalPotentials_FourPhaseRegion)
		
		
		
		self.composition_phasediagram_plot_figure.canvas.mpl_connect('motion_notify_event', self.Hover)
		
		
		
		self.fourphaseregion_annotation = self.composition_phasediagram_plot_drawing.annotate("", xy=(0,0), xytext=(20,20), textcoords="offset points",
																								bbox = dict(boxstyle="round", fc="w"),
																								arrowprops = dict(arrowstyle = "->") )
	
	
	
	
	###############################################################################################
	################## Find All Four-Phase Regions in Compositional Phase Diagram #################
	###############################################################################################
	
	def Find_All_FourPhaseRegions(self):
		
		# Obtain all four-phase regions.
		#	The data structure is as follows:
		#		lines --> List of arrays, each array is 2x2 for ternary (3x2 for quaternary, etc.), column vector represents point on phase diagram.
		#					ex: array([ [0.3, 0.5], [1.0, 0.0] ]) is a line that goes from point [x=0.3, y=1.0] to point [x=0.5, y=0.0]
		#		labels --> Dictionary with point-PDEntry pairs.
		#	Algorithm works as follows:
		#		1) Loop through all lines and get the endpoints.
		#		2) Loop through all points.
		#		3) If the endpoints both draw lines to the third point, then a four-phase region exists between the four points.
		#		4) Loop through all points.
		#		5) If the four points all draw lines to the fourth point, then a four-phase region exists between the four points.
		for line in self.lines:
			endpoint1, endpoint2 = np.transpose(line)
			for point in self.labels.keys():
				point = np.asarray(point)
				
				# Check if the four points form a connected triangle
				if 		( any(np.array_equal(np.transpose(np.vstack((point, endpoint1))), i) for i in self.lines) or any(np.array_equal(np.transpose(np.vstack((endpoint1, point))), i) for i in self.lines) ) \
					and ( any(np.array_equal(np.transpose(np.vstack((point, endpoint2))), i) for i in self.lines) or any(np.array_equal(np.transpose(np.vstack((endpoint2, point))), i) for i in self.lines) ):
					
					
					for point2 in self.labels.keys():
						point2 = np.asarray(point2)
						
						# Check if the four points form a connected tetrahedron
						if 		( any(np.array_equal(np.transpose(np.vstack((point2, endpoint1))), i) for i in self.lines) or any(np.array_equal(np.transpose(np.vstack((endpoint1, point2))), i) for i in self.lines) ) \
							and	( any(np.array_equal(np.transpose(np.vstack((point2, endpoint2))), i) for i in self.lines) or any(np.array_equal(np.transpose(np.vstack((endpoint2, point2))), i) for i in self.lines) ) \
							and ( any(np.array_equal(np.transpose(np.vstack((point2, point))), i) for i in self.lines) or any(np.array_equal(np.transpose(np.vstack((point, point2))), i) for i in self.lines) ):
							
							# Check if four-phase region has already been recorded
							if any( all( list(x) in fourphaseregion for x in [point, point2, endpoint1, endpoint2] ) for fourphaseregion in self.fourphaseregions ):
								continue
							
							# Record four phase region
							four_phase_region = [list(point), list(point2), list(endpoint1), list(endpoint2)]
							self.fourphaseregions.append(four_phase_region)
							
							# Record names of compounds constituting the four-phase region
							four_phase_region_compound_names = []
							for pt in [point, point2, endpoint1, endpoint2]:
								four_phase_region_compound_names.append( self.labels[tuple(pt)].name )
							self.fourphaseregion_names.append(', '.join(four_phase_region_compound_names))
							
							# Find and record the centroid of the four-phase region
							centroid = ( point + point2 + endpoint1 + endpoint2 ) / 4.
							self.fourphaseregion_centroids.append(centroid)
	
	
	
	
	###############################################################################################
	########################## Plot Centroids of All Four-Phase Regions ###########################
	###############################################################################################
	
	def Plot_Centroids(self):
		
		# Check that all four-phase regions have been found
		if self.fourphaseregion_centroids == []:
			return
		
		# Plot all centroids as scatter plot
		self.fourphaseregion_centroids = np.asarray(self.fourphaseregion_centroids)
		self.centroids_plot = self.composition_phasediagram_plot_drawing.scatter(	self.fourphaseregion_centroids[:,0],
																					self.fourphaseregion_centroids[:,1],
																					self.fourphaseregion_centroids[:,2]	)
	
	
	
	
	def Update_Annotation(self, index):
		
		position = self.centroids_plot.get_offsets()[index]
		self.fourphaseregion_annotation.xy = position
		
		text = self.fourphaseregion_names[index]
		
		self.fourphaseregion_annotation.set_text(text)
		self.fourphaseregion_annotation.get_bbox_patch().set_facecolor("w")
		self.fourphaseregion_annotation.get_bbox_patch().set_alpha(1.0)
	
	
	def Hover(self, event):
		is_visible = self.fourphaseregion_annotation.get_visible()
		
		if event.inaxes == self.composition_phasediagram_plot_drawing:
			contains, index = self.centroids_plot.contains(event)
			if contains:
				self.Update_Annotation(index["ind"][0])
				self.fourphaseregion_annotation.set_visible(True)
				self.composition_phasediagram_plot_canvas.draw()
			else:
				if is_visible:
					self.fourphaseregion_annotation.set_visible(False)
					self.composition_phasediagram_plot_canvas.draw()
	
	
	
	
	###############################################################################################
	############# Shade Four-Phase Region in Compositional Phase Diagram when Clicked #############
	###############################################################################################
	
	def Shade_FourPhaseRegion(self, event):
		
		# Check to see that all four-phase regions have been calculated
		if self.fourphaseregions == []:
			return
		
		# Read coordinates of clicked point
		contains, index = self.centroids_plot.contains(event)
		if not contains:
			return
		
		point_index = index["ind"][0]
		
		# Clear previous four-phase region
		if self.fourphaseregion_shade != []:
			for triangle_plot in self.fourphaseregion_shade:
				triangle_plot.remove()
		
		self.fourphaseregion_selected = self.fourphaseregions[point_index]
		
		triangle1 = [self.fourphaseregion_selected[0], self.fourphaseregion_selected[1], self.fourphaseregion_selected[2]]
		triangle2 = [self.fourphaseregion_selected[0], self.fourphaseregion_selected[1], self.fourphaseregion_selected[3]]
		triangle3 = [self.fourphaseregion_selected[0], self.fourphaseregion_selected[2], self.fourphaseregion_selected[3]]
		triangle4 = [self.fourphaseregion_selected[1], self.fourphaseregion_selected[2], self.fourphaseregion_selected[3]]
		
		triangle_plot1 = Poly3DCollection(triangle1, fc='gray', ec='k', alpha=1.0)
		triangle_plot2 = Poly3DCollection(triangle2, fc='gray', ec='k', alpha=1.0)
		triangle_plot3 = Poly3DCollection(triangle3, fc='gray', ec='k', alpha=1.0)
		triangle_plot4 = Poly3DCollection(triangle4, fc='gray', ec='k', alpha=1.0)
		
		self.fourphaseregion_triangle1 = self.composition_phasediagram_plot_drawing.add_collection3d(triangle_plot1)
		self.fourphaseregion_triangle2 = self.composition_phasediagram_plot_drawing.add_collection3d(triangle_plot2)
		self.fourphaseregion_triangle3 = self.composition_phasediagram_plot_drawing.add_collection3d(triangle_plot3)
		self.fourphaseregion_triangle4 = self.composition_phasediagram_plot_drawing.add_collection3d(triangle_plot4)
		
		self.composition_phasediagram_plot_canvas.draw()
		
		self.fourphaseregion_shade = [triangle_plot1, triangle_plot2, triangle_plot3, triangle_plot4]
	
	
	
	###############################################################################################
	################ Calculate Chemical Potentials of Elements in Four-Phase Region ###############
	###############################################################################################
	
	def Calculate_ChemicalPotentials_FourPhaseRegion(self, event):
		
		# Check to see that a four-phase region has been selected
		if self.fourphaseregion_selected is None:
			return
		
		# Get the names of compounds constituting the four-phase region
		fourphaseregion_compounds = []
		for point in self.fourphaseregion_selected:
			entry = self.labels[tuple(point)]	# PDEntry object
			fourphaseregion_compounds.append( entry.name )
		
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
			enthalpies_array[compound_index] = self.compounds_info[compound]["enthalpy"]
		
		# Solve for the delta mu values
		deltamus = np.dot( np.linalg.inv(composition_matrix), enthalpies_array )
		
		# Record delta mu values
		for element, deltamu in zip(self.elements_list, deltamus):
			self.deltamu_values[element] = deltamu








