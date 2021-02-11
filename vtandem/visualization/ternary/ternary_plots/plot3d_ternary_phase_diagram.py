
__author__ = 'Michael_Lidia_Jiaxing_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'


###############################################################################################################################
###############################################################################################################################
##################################################### Introduction ############################################################
###############################################################################################################################
###############################################################################################################################

# 	This script stores the phase stability diagram of a single ternary compound as an object within VTAnDeM.
#
# 	The stability of any compound is dictated by thermodynamic principles. All materials are made of atoms, and multicomponent
#		systems in particular (e.g. systems of two or more types of atoms, like GaAs) can exhibit distinct phases/materials
#		depending on how the atoms are arranged. Therefore, for a given set of atoms, there is a myriad of possible materials 
#		that can be formed by e.g. different synthesis methods. Fundamentally, the material with the lowest formation energy in some
#		environmental condition (e.g. GaAs in a 'Ga-rich' environment) is most and will likely form in a lab. Accordingly, if
#		we are studying some ternary system, there may be other compounds that "compete" with the main compound to become the
#		most stable. For example, if we are interested in studying the ternary compound Hg2GeTe4, in regions that are Cu-rich
#		or Hg-poor (or any combination of the sorts), another compound such as HgTe may be more stable than Hg2GeTe4. It is 
#		therefore crucial for material/device processing purposes that we know the exact environmental conditions necessary to create a 
#		stable version of the ternary compound of interest. The stability of compounds can fortunately be studied using Density 
#		Functional Theory (DFT) calculations.
#	REWRITE ALL OF THIS


###############################################################################################################################
###############################################################################################################################
################################################### Import Libraries ##########################################################
###############################################################################################################################
###############################################################################################################################

import numpy as np
import periodictable
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Path3DCollection, Line3DCollection
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from polyhedron import Hrep

from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

from vtandem.visualization.utils.compound_name import Compound_Name_Formal
from vtandem.visualization.plots.plot_chemicalpotential_phasediagram3d import Plot_ChemicalPotential_PhaseDiagram3D

class ChemicalPotential_Ternary_PhaseDiagram3D(Plot_ChemicalPotential_PhaseDiagram3D):
	
	def __init__(self, parent = None, main_compound = None, first_element = None, second_element = None, third_element = None):
		
		super().__init__(type = "ternary")
		
		# Establish the first, second, and third elements along with the dependent element of the ternary compound.
		self.main_compound  = main_compound
		self.element_x = first_element
		self.element_y = second_element
		self.dependent_element = third_element
		self.elements_list = [self.element_x, self.element_y, self.dependent_element]	# List of elements that gets updated as user selects element order
		self.elements_list_original = [first_element, second_element, third_element]	# Unchanged list of elements for compound naming purposes (e.g. see Compound_Name_Formal)
	
	
	
	###############################################################################################
	################################## Set the dependent element ##################################
	###############################################################################################
	
	def Set_Elements(self, element_x, element_y, dependent_element):
		
		self.element_x = element_x
		self.element_y = element_y
		self.dependent_element = dependent_element
		self.elements_list = [self.element_x, self.element_y, self.dependent_element]







