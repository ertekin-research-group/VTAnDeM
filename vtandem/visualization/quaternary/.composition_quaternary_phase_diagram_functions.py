
__author__ = 'Michael_Lidia_Jiaxing_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'


import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection


class Composition_Quaternary_PhaseDiagram_Functions:
	
	def __init__(self):
		
		pass
	

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
	
	###############################################################################################
	############# Shade Four-Phase Region in Compositional Phase Diagram when Clicked #############
	###############################################################################################
	
	def Shade_FourPhaseRegion(self, event):
		
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
			self.phaseregion_selected = None
			self.composition_phasediagram_plot_canvas.draw()
			return
		
		# Clear previous four-phase region
		if self.phaseregion_shade != []:
			for triangle_plot in self.phaseregion_shade:
				triangle_plot.remove()
		
		# Get selected four-phase region, if any were selected
		self.phaseregion_selected = None
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

