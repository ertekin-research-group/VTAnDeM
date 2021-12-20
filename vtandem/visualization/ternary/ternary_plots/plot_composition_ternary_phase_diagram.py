
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
		
		
		
		# Find all three-phase regions in the quaternary composition space.
		self.Find_All_PhaseRegions()
		
		# Plot all the centroids as a scatter plot (MUST come after generating three-phase regions).
		self.Plot_Centroids()
		
		
		
		self.composition_phasediagram_plot_figure.canvas.mpl_connect('button_press_event', self.Shade_ThreePhaseRegion)
	
	
	
	
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
	
	
	
	###############################################################################################
	############ Shade Three-Phase Region in Compositional Phase Diagram when Clicked #############
	###############################################################################################
	
	def Shade_ThreePhaseRegion(self, event):
		
		# Check to see that all three-phase regions have been calculated
		if self.phase_region_objects == []:
			return
		
		# Read coordinates of clicked point
		point_x = event.xdata
		point_y = event.ydata
		
		# Check that the user presses somewhere on the plot (and not anywhere else)
		if (not isinstance(point_x, float)) or (not isinstance(point_y, float)):
			return
		
		# Delete any previous shadings
		if self.phaseregion_shade != []:
			self.phaseregion_shade.remove()
			self.phaseregion_shade = []
			self.composition_phasediagram_plot_canvas.draw()
		
		# Get selected three-phase region, if any were selected
		self.phaseregion_selected = None
		for three_phase_region in self.phase_region_objects:
			is_in_threephaseregion = mpltPath.Path(three_phase_region.vertices).contains_point([point_x, point_y])
			if is_in_threephaseregion:
				self.phaseregion_selected = three_phase_region
				break
		
		# Shade in the selected three-phase region
		if self.phaseregion_selected is not None:
			self.phaseregion_shade = plt.Polygon(self.phaseregion_selected.vertices, color='gray', alpha=0.5)
			self.composition_phasediagram_plot_drawing.add_patch(self.phaseregion_shade)
			self.composition_phasediagram_plot_canvas.draw()




