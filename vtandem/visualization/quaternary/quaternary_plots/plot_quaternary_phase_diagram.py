
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

from vtandem.visualization.plots.plot_chemicalpotential_phasediagram import ChemicalPotential_PhaseDiagramProjected2D
from vtandem.visualization.plots.plot_chemicalpotential_phasediagram_projectedtripleview import ChemicalPotential_PhaseDiagramProjected2D_TripleView



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
		self.elements_list_original = [self.first_element, self.second_element, self.third_element, self.fourth_element]	# Unchanged list of elements for compound naming purposes

		# Necessary numerical variables
		#self.main_compound_elements_count = {}
		self.main_compound_enthalpy = 0.0			# Enthalpy of the main quaternary compound
		self.phasediagram_endpoints = 0.0			# Endpoints for quaternary phase diagram
		#self.mu4 = 0.0								# Track the fourth species mu value as the user changes it
		self.deltamu = {1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0}




class ChemicalPotential_Quaternary_PhaseDiagramProjected2D_TripleView(ChemicalPotential_PhaseDiagramProjected2D_TripleView):
	
	def __init__(self, parent = None, main_compound = None, first_element = None, second_element = None, third_element = None, fourth_element = None):
		
		# Establish the first, second, third, and fourth species of the quaternary compound.
		# Note that this list is subject to change, depending on what the user chooses.
		self.main_compound  = main_compound
		self.first_element	= first_element
		self.second_element	= second_element
		self.third_element	= third_element
		self.fourth_element	= fourth_element
		self.elements_list  = [self.first_element, self.second_element, self.third_element, self.fourth_element]	# List of elements that gets updated as user selects element order
		self.elements_list_original = [self.first_element, self.second_element, self.third_element, self.fourth_element]	# Unchanged list of elements for compound naming purposes
		
		self.deltamu = {1: 0.0, 2: 0.0, 3: 0.0, 4: 0.0}								# Track the fourth species mu value as the user changes it
		
		
		super().__init__("quaternary")
	
	
	
	###############################################################################################
	################################## Set the dependent element ##################################
	###############################################################################################
	
	def Set_Elements(self, first_element, second_element, third_element, fourth_element):
		
		self.first_element = first_element
		self.second_element = second_element
		self.third_element = third_element
		self.fourth_element = fourth_element
		self.elements_list = [self.first_element, self.second_element, self.third_element, self.fourth_element]

