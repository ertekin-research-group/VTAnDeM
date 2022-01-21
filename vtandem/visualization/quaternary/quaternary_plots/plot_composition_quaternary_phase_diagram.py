
__author__ = 'Michael_Lidia_Jiaxing_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'


# Import compositional phase diagram object
from vtandem.visualization.plots.plot_composition_phasediagram import Plot_Composition_PhaseDiagram
from vtandem.visualization.quaternary.composition_quaternary_phase_diagram_functions import Composition_Quaternary_PhaseDiagram_Functions



class Plot_Composition_Quaternary_PhaseDiagram(Plot_Composition_PhaseDiagram, Composition_Quaternary_PhaseDiagram_Functions):
	
	def __init__(self, main_compound = None, first_element = None, second_element = None, third_element = None, fourth_element = None, compounds_info = None, main_compound_info = None):
		
		# Establish the first, second, third, and fourth species of the quaternary compound.
		self.main_compound  = main_compound
		self.first_element	= first_element
		self.second_element	= second_element
		self.third_element	= third_element
		self.fourth_element	= fourth_element
		self.elements_list  = [self.first_element, self.second_element, self.third_element, self.fourth_element]
		
		"""
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
		"""

		
		# Initialize deltamu values of all species in the compound
		self.deltamu_values = {}
		for element in self.elements_list:
			self.deltamu_values[element] = 0.0
		
		
		print(main_compound_info, "quat script")
		
		# Inherit all variables (plot object, etc.) and functions from parent objects
		Plot_Composition_PhaseDiagram.__init__(self, type = "quaternary", main_compound_info = main_compound_info, compounds_info = compounds_info)
		Composition_Quaternary_PhaseDiagram_Functions.__init__(self)

		# Find all four-phase regions in the quaternary composition space.
		self.Find_All_PhaseRegions()
		
		# Plot all the centroids as a scatter plot (MUST come after generating four-phase regions).
		self.Plot_Centroids()
		
		self.composition_phasediagram_plot_figure.canvas.mpl_connect('button_press_event', self.Shade_FourPhaseRegion)


