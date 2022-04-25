
__author__ = 'Michael_Lidia_Jiaxing_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'


from mpl_toolkits.mplot3d.art3d import Poly3DCollection


# Import compositional phase diagram object
from vtandem.visualization.plots.plot_composition_phasediagram import Plot_Composition_PhaseDiagram



class Plot_Composition_Quaternary_PhaseDiagram(Plot_Composition_PhaseDiagram):
	
	def __init__(self, main_compound = None, first_element = None, second_element = None, third_element = None, fourth_element = None, compounds_info = None, main_compound_info = None):
		
		# Establish the first, second, third, and fourth species of the quaternary compound.
		self.main_compound  = main_compound
		self.first_element	= first_element
		self.second_element	= second_element
		self.third_element	= third_element
		self.fourth_element	= fourth_element
		self.elements_list  = [self.first_element, self.second_element, self.third_element, self.fourth_element]
		
		# Initialize deltamu values of all species in the compound
		self.deltamu_values = {}
		for element in self.elements_list:
			self.deltamu_values[element] = 0.0
		
		# Inherit all variables (plot object, etc.) and functions from parent objects
		Plot_Composition_PhaseDiagram.__init__(self, type = "quaternary", main_compound_info = main_compound_info, compounds_info = compounds_info)

		# Find all four-phase regions in the quaternary composition space.
		self.Find_All_PhaseRegions()
		
		# Plot all the centroids as a scatter plot (MUST come after generating four-phase regions).
		self.Plot_Centroids()
		
		self.composition_phasediagram_plot_figure.canvas.mpl_connect('button_press_event', self.Shade_FourPhaseRegion)



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


