
__author__ = 'Michael_Lidia_Jiaxing_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'


import numpy as np
import itertools
import copy
import periodictable
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

import PyQt5
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

from vtandem.visualization.plots.plot_chemicalpotential_phasediagram import ChemicalPotential_PhaseDiagramProjected2D
from vtandem.visualization.plots.plot_chemicalpotential_phasediagram_projectedtripleview import ChemicalPotential_PhaseDiagramProjected2D_TripleView



class ChemicalPotential_Ternary_PhaseDiagramProjected2D(ChemicalPotential_PhaseDiagramProjected2D):
	
	def __init__(self, parent = None, main_compound = None, first_element = None, second_element = None, third_element = None):
		
		super().__init__(main_compound, "ternary")
		
		# Establish the first, second, and third species of the quaternary compound.
		# Note that this list is subject to change, depending on what the user chooses.
		self.main_compound  = main_compound
		self.first_element  = first_element
		self.second_element = second_element
		self.third_element  = third_element
		self.elements_list  = [self.first_element, self.second_element, self.third_element]		# Species list (order MAY change)
		
		# Information about main ternary compound
		self.main_compound_number_first_specie  = 0	# Number of each specie in ternary compound
		self.main_compound_number_second_specie = 0
		self.main_compound_number_third_specie  = 0
		self.main_compound_enthalpy = 0.0			# Enthalpy of ternary compound
		self.phasediagram_endpoints = 0.0			# Endpoints for ternary phase diagram
		self.deltamu = {1: 0.0, 2: 0.0, 3: 0.0}





class ChemicalPotential_Ternary_PhaseDiagramProjected2D_TripleView(ChemicalPotential_PhaseDiagramProjected2D_TripleView):
	
	def __init__(self, parent = None, main_compound = None, first_element = None, second_element = None, third_element = None):
		
		# Establish the first, second, and third species of the quaternary compound.
		# Note that this list is subject to change, depending on what the user chooses.
		self.main_compound  = main_compound
		self.first_element	= first_element
		self.second_element	= second_element
		self.third_element	= third_element
		self.elements_list  = [self.first_element, self.second_element, self.third_element]	# List of elements that gets updated as user selects element order
		self.elements_list_original = [self.first_element, self.second_element, self.third_element]	# Unchanged list of elements for compound naming purposes
		
		self.deltamu = {1: 0.0, 2: 0.0, 3: 0.0}
		
		super().__init__("ternary")
	
	
	###############################################################################################
	################################## Set the dependent element ##################################
	###############################################################################################
	
	def Set_Elements(self, first_element, second_element, third_element):
		
		self.first_element = first_element
		self.second_element = second_element
		self.third_element = third_element
		self.elements_list = [self.first_element, self.second_element, self.third_element]




