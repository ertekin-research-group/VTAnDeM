
__author__ = 'Michael_Lidia_Jiaxing_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'


###############################################################################################################################
###############################################################################################################################
##################################################### Introduction ############################################################
###############################################################################################################################
###############################################################################################################################

# 	This script stores the phase stability diagram of a single quaternary compound as an object within VTAnDeM.
#
# 	The stability of any compound is dictated by thermodynamic principles. All materials are made of atoms, and multicomponent
#		systems in particular (e.g. systems of two or more types of atoms, like GaAs) can exhibit distinct phases/materials
#		depending on how the atoms are arranged. Therefore, for a given set of atoms, there is a myriad of possible materials 
#		that can be formed by e.g. different synthesis methods. Fundamentally, the material with the lowest formation energy in some
#		environmental condition (e.g. GaAs in a 'Ga-rich' environment) is most and will likely form in a lab. Accordingly, if
#		we are studying some quaternary system, there may be other compounds that "compete" with the main compound to become the
#		most stable. For example, if we are interested in studying the quaternary compound Cu2HgGeTe4, in regions that are Cu-rich
#		or Hg-poor (or any combination of the sorts), another compound such as CuTe may be more stable than Cu2HgGeTe4. It is 
#		therefore crucial for material/device processing purposes that we know the exact environmental conditions necessary to create a 
#		stable version of the quaternary compound of interest. The stability of compounds can fortunately be studied using Density 
#		Functional Theory (DFT) calculations.



###############################################################################################################################
###############################################################################################################################
################################################### Import Libraries ##########################################################
###############################################################################################################################
###############################################################################################################################

import numpy as np
import itertools
import copy
import periodictable
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

import PyQt5
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

from vtandem.visualization.utils.chemicalpotential_phasediagram import Calculate_PhaseDiagram_Projected2D
from vtandem.visualization.utils.compound_name import Compound_Name_Formal

from vtandem.visualization.plots.plot_chemicalpotential_phasediagram import ChemicalPotential_PhaseDiagramProjected2D



class ChemicalPotential_Quaternary_PhaseDiagramProjected2D(ChemicalPotential_PhaseDiagramProjected2D):
	
	def __init__(self, parent = None, main_compound = None, first_element = None, second_element = None, third_element = None, fourth_element = None):
		
		super().__init__(main_compound, "quaternary")
		
		# Establish the first, second, third, and fourth species of the quaternary compound.
		# Note that this list is subject to change, depending on what the user chooses.
		self.main_compound  = main_compound
		self.first_element	= first_element
		self.second_element	= second_element
		self.third_element	= third_element
		self.fourth_element	= fourth_element
		self.elements_list  = [self.first_element, self.second_element, self.third_element, self.fourth_element]	# List of elements that gets updated as user selects element order
		self.elements_list_original = [self.first_element, self.second_element, self.third_element, self.fourth_element]	# Unchanged list of elements for compound naming purposes (e.g. see Compound_Name_Formal)
		
		# Necessary numerical variables
		self.main_compound_elements_count = {}
		self.main_compound_enthalpy = 0.0			# Enthalpy of the main quaternary compound
		self.phasediagram_endpoints = 0.0			# Endpoints for quaternary phase diagram
		#self.mu4 = 0.0								# Track the fourth species mu value as the user changes it
		self.deltamu = {1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0}




class ChemicalPotential_Quaternary_PhaseDiagramProjected2D_TripleView(QWidget):
	
	def __init__(self, parent = None, main_compound = None, first_element = None, second_element = None, third_element = None, fourth_element = None):
		
		"""
		# All elements in the periodic table
		self.all_elements = []
		for element in periodictable.elements:
			self.all_elements.append(str(element))
		
		# Font description for phase stability diagram plot
		self.font = {'family': 'sans-serif',
				'color':  'black',
				'weight': 'normal',
				'size': 12 }
		"""
		
		# Establish the first, second, third, and fourth species of the quaternary compound.
		# Note that this list is subject to change, depending on what the user chooses.
		self.main_compound  = main_compound
		self.first_element	= first_element
		self.second_element	= second_element
		self.third_element	= third_element
		self.fourth_element	= fourth_element
		self.elements_list  = [self.first_element, self.second_element, self.third_element, self.fourth_element]	# List of elements that gets updated as user selects element order
		self.elements_list_original = [self.first_element, self.second_element, self.third_element, self.fourth_element]	# Unchanged list of elements for compound naming purposes (e.g. see Compound_Name_Formal)
		
		"""
		# Store all extracted DFT data
		self.compounds_info = {}
		"""
		
		# Necessary numerical variables
		"""
		self.main_compound_number_first_specie  = 0	# Number of each specie in the main quaternary compound
		self.main_compound_number_second_specie = 0
		self.main_compound_number_third_specie  = 0
		self.main_compound_number_fourth_specie = 0
		"""
		
		"""
		self.main_compound_elements_count = {}
		self.main_compound_enthalpy = 0.0			# Enthalpy of the main quaternary compound
		self.phasediagram_endpoints = 0.0			# Endpoints for quaternary phase diagram
		#self.mu4 = 0.0								# Track the fourth species mu value as the user changes it
		"""
		self.deltamu = {1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0}								# Track the fourth species mu value as the user changes it
		
		"""
		# Quaternary phase diagram object
		self.tripleview_phase_diagram_plot_figure = plt.figure(figsize=(7,7))
		
		self.tripleview_phase_diagram_plot_drawing12 = self.tripleview_phase_diagram_plot_figure.add_subplot(221)
		self.tripleview_phase_diagram_plot_drawing13 = self.tripleview_phase_diagram_plot_figure.add_subplot(222)
		self.tripleview_phase_diagram_plot_drawing23 = self.tripleview_phase_diagram_plot_figure.add_subplot(223)
		
		self.tripleview_phase_diagram_plot_canvas = FigureCanvas(self.tripleview_phase_diagram_plot_figure)
		
		
		
		# Necessary plot-related variables
		self.main_compound_plot12 = None				# Main compound plot (holds plot object)
		self.main_compound_plot13 = None				# Main compound plot (holds plot object)
		self.main_compound_plot23 = None				# Main compound plot (holds plot object)
		self.competing_compound_plots12 = {}			# Competing compound plots (holds plot objects)
		self.competing_compound_plots13 = {}			# Competing compound plots (holds plot objects)
		self.competing_compound_plots23 = {}			# Competing compound plots (holds plot objects)
		self.competing_compounds_colorwheel = {}		# Legend
		"""
	
	
	
	###############################################################################################
	################################## Set the dependent element ##################################
	###############################################################################################
	
	def Set_Elements(self, first_element, second_element, third_element, fourth_element):
		
		self.first_element = first_element
		self.second_element = second_element
		self.third_element = third_element
		self.fourth_element = fourth_element
		self.elements_list = [self.first_element, self.second_element, self.third_element, self.fourth_element]
	
	
	
	
	
	def Update_PhaseDiagram_Object(self):
		
		# Number of elements in main compound
		"""
		self.main_compound_number_first_specie = self.compounds_info[self.main_compound][self.first_element]	# Number of first species in quaternary compound
		self.main_compound_number_second_specie = self.compounds_info[self.main_compound][self.second_element]	# Number of second species in quaternary compound
		self.main_compound_number_third_specie = self.compounds_info[self.main_compound][self.third_element]	# Number of third species in quaternary compound
		self.main_compound_number_fourth_specie = self.compounds_info[self.main_compound][self.fourth_element]	# Number of fourth species in quaternary compound
		"""
		
		"""
		for element in self.elements_list:
			main_compound_elements_count[element] = 
		"""
		
		"""
		# Enthalpy of quaternary compound
		#self.main_compound_enthalpy = self.compounds_info[self.main_compound]["enthalpy"]	# Enthalpy of quaternary compound
		enthalpy_tracker = defects_data["Bulk"]["dft_BulkEnergy"]
		for element in self.elements_list:
			enthalpy_tracker -= self.main_compound_elements_count[element] * self.compounds_info[element]["mu0"]
		self.main_compound_enthalpy = enthalpy_tracker
		
		# Endpoints of phase diagram
		self.phasediagram_endpoints = min(self.main_compound_enthalpy/self.main_compound_number_first_specie, self.main_compound_enthalpy/self.main_compound_number_second_specie, self.main_compound_enthalpy/self.main_compound_number_third_specie, self.main_compound_enthalpy/self.main_compound_number_fourth_specie)
		"""
	
	"""
	def Update_PhaseDiagram_Plot_Axes(self):
		
		self.tripleview_phase_diagram_plot_drawing12.set_xlim(self.phasediagram_endpoints, 0.0)
		self.tripleview_phase_diagram_plot_drawing12.set_ylim(self.phasediagram_endpoints, 0.0)
		self.tripleview_phase_diagram_plot_drawing12.set_xlabel("$\Delta\mu_{"+self.first_element+"}$ (eV)",fontdict=self.font)
		self.tripleview_phase_diagram_plot_drawing12.set_ylabel("$\Delta\mu_{"+self.second_element+"}$ (eV)",fontdict=self.font,rotation=270,labelpad=20)
		self.tripleview_phase_diagram_plot_drawing12.xaxis.tick_top()
		self.tripleview_phase_diagram_plot_drawing12.yaxis.tick_right()
		self.tripleview_phase_diagram_plot_drawing12.xaxis.set_label_position("top")
		self.tripleview_phase_diagram_plot_drawing12.yaxis.set_label_position("right")
		self.tripleview_phase_diagram_plot_drawing12.spines['left'].set_visible(False)
		self.tripleview_phase_diagram_plot_drawing12.spines['bottom'].set_visible(False)
		self.tripleview_phase_diagram_plot_drawing12.set_aspect("equal")
		
		self.tripleview_phase_diagram_plot_drawing13.set_xlim(self.phasediagram_endpoints, 0.0)
		self.tripleview_phase_diagram_plot_drawing13.set_ylim(self.phasediagram_endpoints, 0.0)
		self.tripleview_phase_diagram_plot_drawing13.set_xlabel("$\Delta\mu_{"+self.first_element+"}$ (eV)",fontdict=self.font)
		self.tripleview_phase_diagram_plot_drawing13.set_ylabel("$\Delta\mu_{"+self.third_element+"}$ (eV)",fontdict=self.font,rotation=270,labelpad=20)
		self.tripleview_phase_diagram_plot_drawing13.xaxis.tick_top()
		self.tripleview_phase_diagram_plot_drawing13.yaxis.tick_right()
		self.tripleview_phase_diagram_plot_drawing13.xaxis.set_label_position("top")
		self.tripleview_phase_diagram_plot_drawing13.yaxis.set_label_position("right")
		self.tripleview_phase_diagram_plot_drawing13.spines['left'].set_visible(False)
		self.tripleview_phase_diagram_plot_drawing13.spines['bottom'].set_visible(False)
		self.tripleview_phase_diagram_plot_drawing13.set_aspect("equal")
		
		self.tripleview_phase_diagram_plot_drawing23.set_xlim(self.phasediagram_endpoints, 0.0)
		self.tripleview_phase_diagram_plot_drawing23.set_ylim(self.phasediagram_endpoints, 0.0)
		self.tripleview_phase_diagram_plot_drawing23.set_xlabel("$\Delta\mu_{"+self.second_element+"}$ (eV)",fontdict=self.font)
		self.tripleview_phase_diagram_plot_drawing23.set_ylabel("$\Delta\mu_{"+self.third_element+"}$ (eV)",fontdict=self.font,rotation=270,labelpad=20)
		self.tripleview_phase_diagram_plot_drawing23.xaxis.tick_top()
		self.tripleview_phase_diagram_plot_drawing23.yaxis.tick_right()
		self.tripleview_phase_diagram_plot_drawing23.xaxis.set_label_position("top")
		self.tripleview_phase_diagram_plot_drawing23.yaxis.set_label_position("right")
		self.tripleview_phase_diagram_plot_drawing23.spines['left'].set_visible(False)
		self.tripleview_phase_diagram_plot_drawing23.spines['bottom'].set_visible(False)
		self.tripleview_phase_diagram_plot_drawing23.set_aspect("equal")
		
		self.tripleview_phase_diagram_plot_figure.tight_layout()
	"""
	
	
	def Establish_CompetingCompounds_Colorwheel(self):
		
		color_counter = 0
		
		# Loop through all compounds in the database
		for competing_compound in self.compounds_info.keys():	
			
			# Skip if compound is either the main compound or one of the elements
			if (competing_compound in self.all_elements) or (competing_compound == self.main_compound):
				continue
			
			self.competing_compounds_colorwheel[competing_compound] = self.colorscheme(color_counter)
			color_counter += 1
	
	
	
	def Plot_PhaseDiagrams(self):
		
		# First-second elements
		main_compound_deltamu_first_element12, main_compound_stability_limit12, \
			competing_compounds_deltamu_first_element_limit12, competing_compounds_deltamu_second_element_limit12, \
			main_compound_deltamu_first_element_cutoff12, stability_minimum_cutoff12, stability_maximum_cutoff12 \
			= Calculate_PhaseDiagram_Projected2D(self.main_compound, {1: self.first_element, 2: self.second_element, 3: self.third_element, 4: self.fourth_element}, self.compounds_info, self.deltamu)
		
		#print(self.competing_compounds_colorwheel)
		
		try:
			self.main_compound_plot12.set_data(main_compound_deltamu_first_element12, main_compound_stability_limit12)
		except:
			self.main_compound_plot12, = self.tripleview_phase_diagram_plot_drawing12.plot(main_compound_deltamu_first_element12, main_compound_stability_limit12, color='k')
		
		for competing_compound in self.compounds_info.keys():
			
			# Skip if compound is either the main compound or one of the elements
			if (competing_compound in self.all_elements) or (competing_compound == self.main_compound):
				continue
			
			try:
				# See if plot exists so it just needs to be updated
				self.competing_compound_plots12[competing_compound].set_data(competing_compounds_deltamu_first_element_limit12[competing_compound], competing_compounds_deltamu_second_element_limit12[competing_compound])
			except:
				# If the plot doesn't exist initially, then create it
				self.competing_compound_plots12[competing_compound], = self.tripleview_phase_diagram_plot_drawing12.plot(competing_compounds_deltamu_first_element_limit12[competing_compound], competing_compounds_deltamu_second_element_limit12[competing_compound], label=Compound_Name_Formal(competing_compound, self.compounds_info, "latex"), color=self.competing_compounds_colorwheel[Compound_Name_Formal(competing_compound, self.compounds_info, "unicode")])
		
		try:
			self.phase_stability_region12.remove()
			self.phase_stability_region12 = self.tripleview_phase_diagram_plot_drawing12.fill_between(main_compound_deltamu_first_element_cutoff12, stability_maximum_cutoff12, stability_minimum_cutoff12, facecolor='0.75')
		except:
			self.phase_stability_region12 = self.tripleview_phase_diagram_plot_drawing12.fill_between(main_compound_deltamu_first_element_cutoff12, stability_maximum_cutoff12, stability_minimum_cutoff12, facecolor='0.75')
		
		
		
		
		
		
		# First-third elements
		main_compound_deltamu_first_element13, main_compound_stability_limit13, \
			competing_compounds_deltamu_first_element_limit13, competing_compounds_deltamu_second_element_limit13, \
			main_compound_deltamu_first_element_cutoff13, stability_minimum_cutoff13, stability_maximum_cutoff13 \
			= Calculate_PhaseDiagram_Projected2D(self.main_compound, {1: self.first_element, 2: self.third_element, 3: self.second_element, 4: self.fourth_element}, self.compounds_info, self.deltamu)
		
		try:
			self.main_compound_plot13.set_data(main_compound_deltamu_first_element13, main_compound_stability_limit13)
		except:
			self.main_compound_plot13, = self.tripleview_phase_diagram_plot_drawing13.plot(main_compound_deltamu_first_element13, main_compound_stability_limit13, color='k')
		
		for competing_compound in self.compounds_info.keys():
			
			# Skip if compound is either the main compound or one of the elements
			if (competing_compound in self.all_elements) or (competing_compound == self.main_compound):
				continue
			
			try:
				# See if plot exists so it just needs to be updated
				self.competing_compound_plots13[competing_compound].set_data(competing_compounds_deltamu_first_element_limit13[competing_compound], competing_compounds_deltamu_second_element_limit13[competing_compound])
			except:
				# If the plot doesn't exist initially, then create it
				self.competing_compound_plots13[competing_compound], = self.tripleview_phase_diagram_plot_drawing13.plot(competing_compounds_deltamu_first_element_limit13[competing_compound], competing_compounds_deltamu_second_element_limit13[competing_compound], label=Compound_Name_Formal(competing_compound, self.compounds_info, "latex"), color=self.competing_compounds_colorwheel[Compound_Name_Formal(competing_compound, self.compounds_info, "unicode")])
		
		try:
			self.phase_stability_region13.remove()
			self.phase_stability_region13 = self.tripleview_phase_diagram_plot_drawing13.fill_between(main_compound_deltamu_first_element_cutoff13, stability_maximum_cutoff13, stability_minimum_cutoff13, facecolor='0.75')
		except:
			self.phase_stability_region13 = self.tripleview_phase_diagram_plot_drawing13.fill_between(main_compound_deltamu_first_element_cutoff13, stability_maximum_cutoff13, stability_minimum_cutoff13, facecolor='0.75')
		
		
		
		
		
		# Second-third elements
		main_compound_deltamu_first_element23, main_compound_stability_limit23, \
			competing_compounds_deltamu_first_element_limit23, competing_compounds_deltamu_second_element_limit23, \
			main_compound_deltamu_first_element_cutoff23, stability_minimum_cutoff23, stability_maximum_cutoff23 \
			= Calculate_PhaseDiagram_Projected2D(self.main_compound, {1: self.second_element, 2: self.third_element, 3: self.first_element, 4: self.fourth_element}, self.compounds_info, self.deltamu)
		
		try:
			self.main_compound_plot23.set_data(main_compound_deltamu_first_element23, main_compound_stability_limit23)
		except:
			self.main_compound_plot23, = self.tripleview_phase_diagram_plot_drawing23.plot(main_compound_deltamu_first_element23, main_compound_stability_limit23, color='k')
		
		for competing_compound in self.compounds_info.keys():
			
			# Skip if compound is either the main compound or one of the elements
			if (competing_compound in self.all_elements) or (competing_compound == self.main_compound):
				continue
			
			try:
				# See if plot exists so it just needs to be updated
				self.competing_compound_plots23[competing_compound].set_data(competing_compounds_deltamu_first_element_limit23[competing_compound], competing_compounds_deltamu_second_element_limit23[competing_compound])
			except:
				# If the plot doesn't exist initially, then create it
				self.competing_compound_plots23[competing_compound], = self.tripleview_phase_diagram_plot_drawing23.plot(competing_compounds_deltamu_first_element_limit23[competing_compound], competing_compounds_deltamu_second_element_limit23[competing_compound], label=Compound_Name_Formal(competing_compound, self.compounds_info, "latex"), color=self.competing_compounds_colorwheel[Compound_Name_Formal(competing_compound, self.compounds_info, "unicode")])
		
		try:
			self.phase_stability_region23.remove()
			self.phase_stability_region23 = self.tripleview_phase_diagram_plot_drawing23.fill_between(main_compound_deltamu_first_element_cutoff23, stability_maximum_cutoff23, stability_minimum_cutoff23, facecolor='0.75')
		except:
			self.phase_stability_region23 = self.tripleview_phase_diagram_plot_drawing23.fill_between(main_compound_deltamu_first_element_cutoff23, stability_maximum_cutoff23, stability_minimum_cutoff23, facecolor='0.75')
		
		
		
		# Draw the phase diagram
		self.tripleview_phase_diagram_plot_canvas.draw()








