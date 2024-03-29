
__author__ = 'Michael_Lidia_Jiaxing_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'


import numpy as np

import PyQt5
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

from vtandem.visualization.ternary.ternary_plots.plot_ternary_phase_diagram import ChemicalPotential_Ternary_PhaseDiagramProjected2D_TripleView
from vtandem.visualization.ternary.ternary_plots.plot3d_ternary_phase_diagram import ChemicalPotential_Ternary_PhaseDiagram3D



class Tab_Ternary_PhaseDiagram3D(object):
	
	def __init__(self, parent = None, main_compound = None, first_element = None, second_element = None, third_element = None, compounds_info = None, main_compound_info = None):	# User specifies the main compound and its constituents
		
		self.main_compound = main_compound
		self.first_element = first_element
		self.second_element = second_element
		self.third_element = third_element
		self.elements_list = [self.first_element, self.second_element, self.third_element]					# Species list (order MAY change)
		
		self.compounds_info = compounds_info
		self.main_compound_info = main_compound_info
		
		
		
		# Get enthalpy of main compound (for fourth element slider bar)
		enthalpy_tracker = self.main_compound_info["dft_BulkEnergy"]
		for element in self.elements_list:
			enthalpy_tracker -= self.main_compound_info["dft_"+element] * self.compounds_info[element]["mu0"]
		self.main_compound_enthalpy = enthalpy_tracker
		
		
		
		
		self.PhaseDiagram3D = ChemicalPotential_Ternary_PhaseDiagram3D(self, main_compound = self.main_compound, first_element = self.first_element, second_element = self.second_element, third_element = self.third_element)
		self.PhaseDiagram3D.main_compound_enthalpy = self.main_compound_enthalpy
		self.PhaseDiagram3D.compounds_info = self.compounds_info
		self.PhaseDiagram3D.main_compound_info = self.main_compound_info
		
		
		
		self.PhaseDiagram2D_TripleView = ChemicalPotential_Ternary_PhaseDiagramProjected2D_TripleView(self, main_compound = main_compound, first_element = self.second_element, second_element = self.third_element, third_element = self.first_element)
		self.PhaseDiagram2D_TripleView.main_compound_enthalpy = self.main_compound_enthalpy
		self.PhaseDiagram2D_TripleView.compounds_info = self.compounds_info
		self.PhaseDiagram2D_TripleView.main_compound_info = self.main_compound_info
		self.PhaseDiagram2D_TripleView.phasediagram_endpoints = min(self.main_compound_enthalpy/self.main_compound_info["dft_"+self.first_element], self.main_compound_enthalpy/self.main_compound_info["dft_"+self.second_element], self.main_compound_enthalpy/self.main_compound_info["dft_"+self.third_element])
		#self.PhaseDiagram2D_TripleView.Update_PhaseDiagram_Object()
		self.PhaseDiagram2D_TripleView.Update_PhaseDiagram_Plot_Axes()
		
		
		
		###############################################################################################
		###############################################################################################
		#################################### Initialize second tab ####################################
		###############################################################################################
		###############################################################################################
		
		self.tab2 = QWidget()
		self.tab2_layout = QHBoxLayout(self.tab2)
		
		
		###############################################################################################
		############################# 3D Chemical Potential Phase Diagram #############################
		###############################################################################################
		
		###### Main chemical potential phase diagram window widget
		self.chemicalpotential_phasediagram3d_window = QWidget()															# One of the main sub-widgets is where the user defines the settings of the plots.
		self.chemicalpotential_phasediagram3d_window_layout = QVBoxLayout(self.chemicalpotential_phasediagram3d_window)		# The settings should be placed on top of one another, i.e. vertically.
		
		"""
		# (WIDGET) Title of compound
		compound_title_formal = Compound_Name_Formal(main_compound, self.compounds_info, "unicode")		# Generate Latex-readable version of compound name
		self.compound_title = QLabel(compound_title_formal)									# QLabel is a widget that displays text
		self.compound_title.setAlignment(Qt.AlignCenter)									# Align the text to center
		self.compound_title_font = QFont("sans-serif", 24, QFont.Bold) 						# Declare font
		self.compound_title.setFont(self.compound_title_font)								# Set the font for the QLabel text
		self.chemicalpotential_phasediagram3d_window_layout.addWidget(self.compound_title)	# Add the widget to the "main" widget grid layout
		"""
		# (WIDGET) Title
		self.chemicalpotentialPD_name = QLabel("3D Phase Diagram")							# QLabel is a widget that displays text
		self.chemicalpotentialPD_name.setAlignment(Qt.AlignCenter)								# Align the text to center
		self.chemicalpotentialPD_name_font = QFont("sans-serif", 24, QFont.Bold) 					# Declare font
		self.chemicalpotentialPD_name.setFont(self.chemicalpotentialPD_name_font)		# Set the font for the QLabel text
		self.chemicalpotential_phasediagram3d_window_layout.addWidget(self.chemicalpotentialPD_name)	# Add the widget to the "main" widget grid layout
		
		# (WIDGET) Chemical potential phase diagram plot object
		self.chemicalpotential_phase_diagram_plot = self.PhaseDiagram3D.chemicalpotential_phasediagram_plot_canvas
		self.chemicalpotential_phasediagram3d_window_layout.addWidget(self.chemicalpotential_phase_diagram_plot)
		
		# (WIDGET) Save 3D phase diagram as figure
		self.phasediagram3d_savefigure_button = QPushButton("Save 3D Phase Diagram Figure")
		self.phasediagram3d_savefigure_button.clicked[bool].connect(lambda: self.PhaseDiagram3D.SaveFigure())
		self.chemicalpotential_phasediagram3d_window_layout.addWidget(self.phasediagram3d_savefigure_button)
		
		# (WIDGET) Generate a 3D rotating phase diagram animation
		self.generate_phasediagram3d_animation_button = QPushButton("360 Rotation Animation")
		self.generate_phasediagram3d_animation_button.clicked[bool].connect(self.PhaseDiagram3D.Make_SpacePotato_Rotation_Animation)
		self.chemicalpotential_phasediagram3d_window_layout.addWidget(self.generate_phasediagram3d_animation_button)
		
		# Add 3D phase diagram window to tab2
		self.tab2_layout.addWidget(self.chemicalpotential_phasediagram3d_window)
		
		
		###############################################################################################
		###################### Three Projected Chemical Potential Phase Diagrams ######################
		###############################################################################################
		
		
		self.chemicalpotentialPD_projected_window = QWidget()
		self.chemicalpotentialPD_projected_window_layout = QVBoxLayout(self.chemicalpotentialPD_projected_window)
		
		# (WIDGET) Title
		self.chemicalpotentialPD_projected_name = QLabel("Projected Phase Diagram")							# QLabel is a widget that displays text
		self.chemicalpotentialPD_projected_name.setAlignment(Qt.AlignCenter)								# Align the text to center
		self.chemicalpotentialPD_projected_name_font = QFont("sans-serif", 24, QFont.Bold) 					# Declare font
		self.chemicalpotentialPD_projected_name.setFont(self.chemicalpotentialPD_projected_name_font)		# Set the font for the QLabel text
		self.chemicalpotentialPD_projected_window_layout.addWidget(self.chemicalpotentialPD_projected_name)	# Add the widget to the "main" widget grid layout
		
		# (WIDGET) Three projected chemical potential phase diagrams
		self.phase_diagram_plot_2d_tripleview = self.PhaseDiagram2D_TripleView.tripleview_phase_diagram_plot_canvas
		self.chemicalpotentialPD_projected_window_layout.addWidget(self.phase_diagram_plot_2d_tripleview)
		
		# Add window to tab2
		self.tab2_layout.addWidget(self.chemicalpotentialPD_projected_window)
		
		
		# Draw phase diagrams
		self.Generate_PhaseDiagram3D_Function()
	
	
	
	def Generate_PhaseDiagram3D_Function(self):
		
		# This function specifies what happens when the user clicks the "Generate Plot" button. The 
		#	only thing it needs to check is whether the elements chosen as the first, second, and 
		#	third species are unique.
		
		if len(self.elements_list) > len(set(self.elements_list)):	# Every time someone clicks the species' buttons, the self.elements_list gets updated.
																	#	The "set" function checks the self.elements_list and omits any that are repeated.
																	#	If any are repeated, then the chosen species are not unique.
			QMessageBox.about(self, "WARNING", "Pick UNIQUE elements!")
			return
		
		
		self.PhaseDiagram3D.chemicalpotential_phasediagram_plot_axes.clear()
		
		#self.PhaseDiagram3D.Set_Elements(element_x = self.first_element, element_y = self.second_element, dependent_element = self.third_element)
		self.PhaseDiagram3D.Set_Elements(element_x = self.first_element, element_y = self.second_element, element_z = self.third_element)
		
		self.PhaseDiagram3D.Draw_PhaseDiagram3D()
		
		self.PhaseDiagram2D_TripleView.Set_Elements(first_element = self.first_element, second_element = self.second_element, third_element = self.third_element)
		#self.PhaseDiagram2D_TripleView.Update_PhaseDiagram_Object()
		self.PhaseDiagram2D_TripleView.Update_PhaseDiagram_Plot_Axes()
		
		self.PhaseDiagram2D_TripleView.competing_compounds_colorwheel = self.PhaseDiagram3D.competing_compounds_colorwheel
		
		self.PhaseDiagram2D_TripleView.Plot_PhaseDiagrams()




