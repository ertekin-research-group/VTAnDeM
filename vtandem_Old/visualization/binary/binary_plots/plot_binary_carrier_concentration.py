
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'
__author__ = 'Michael_Lidia_Jiaxing_Elif'


"""
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy import integrate
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

import PyQt5
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

# Import functions for calculating carrier concentration
from vtandem.visualization.utils.carrier_concentration import Calculate_CarrierConcentration
"""

# Import carrier concentration plot object
from vtandem.visualization.plots.plot_carrier_concentration import Plot_CarrierConcentration



class Plot_Binary_Carrier_Concentration(Plot_CarrierConcentration):
	
	def __init__(self, parent = None, main_compound = None, first_element = None, second_element = None):
		
		"""
		# Font description for phase stability diagram plot
		self.font = {'family': 'sans-serif',
				'color':  'black',
				'weight': 'normal',
				'size': 14 }
		"""
		
		# Inherit all variables (plot object, etc.) from parent object (Plot_CarrierConcentration)
		super().__init__()
		
		# Establish the first and second species of the binary compound.
		# Note that this list is subject to change, depending on what the user chooses.
		self.main_compound  = main_compound
		self.first_element	= first_element
		self.second_element	= second_element
		self.elements_list  = [self.first_element, self.second_element]
		
		
		# Keep track of chemical potential values
		self.mu_elements = {self.first_element: {"mu0": 0.0, "deltamu": 0.0},
							self.second_element: {"mu0": 0.0, "deltamu": 0.0} }
		
		# Number of each specie in the main binary compound
		self.main_compound_number_first_specie  = 0
		self.main_compound_number_second_specie = 0
		self.number_species = {	self.first_element: self.main_compound_number_first_specie,
								self.second_element: self.main_compound_number_second_specie }
		
		# Store all extracted DFT data
		self.first_element_mu0 = 0.0
		self.second_element_mu0 = 0.0








