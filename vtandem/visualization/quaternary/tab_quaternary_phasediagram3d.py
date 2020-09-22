
__author__ = 'Michael_Lidia_Jiaxing_Benita_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'


import numpy as np

import PyQt5
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

from vtandem.visualization.quaternary.quaternary_scripts.plot_quaternary_phase_diagram import ChemicalPotential_Quaternary_PhaseDiagramProjected2D_TripleView
from vtandem.visualization.quaternary.quaternary_scripts.plot3d_quaternary_phase_diagram import ChemicalPotential_Quaternary_PhaseDiagram3D

#from vtandem.visualization.widgets.phasediagram_window_chemicalpotential import ChemicalPotential_PhaseDiagram_Window

from vtandem.visualization.compound_name import Compound_Name_Formal



class Tab_PhaseDiagram3D(object):
	
	def __init__(self, parent = None, main_compound = None, first_element = None, second_element = None, third_element = None, fourth_element = None, compounds_info = None):	# User specifies the main compound and its constituents
		
		self.main_compound = main_compound
		self.first_element = first_element
		self.second_element = second_element
		self.third_element = third_element
		self.fourth_element = fourth_element
		self.elements_list = [self.first_element, self.second_element, self.third_element, self.fourth_element]					# Species list (order MAY change)
		
		self.compounds_info = compounds_info
		
		self.PhaseDiagram3D = ChemicalPotential_Quaternary_PhaseDiagram3D(self, main_compound = self.main_compound, first_element = self.first_element, second_element = self.second_element, third_element = self.third_element, fourth_element = self.fourth_element)
		self.PhaseDiagram3D.compounds_info = self.compounds_info
		
		self.PhaseDiagram2D_TripleView = ChemicalPotential_Quaternary_PhaseDiagramProjected2D_TripleView(self, main_compound = main_compound, first_element = self.second_element, second_element = self.third_element, third_element = self.first_element, fourth_element = self.fourth_element)
		self.PhaseDiagram2D_TripleView.compounds_info = self.compounds_info
		self.PhaseDiagram2D_TripleView.Update_PhaseDiagram_Object()
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
		
		# (WIDGET) Title of compound
		compound_title_formal = Compound_Name_Formal(main_compound, self.compounds_info, "unicode")		# Generate Latex-readable version of compound name
		self.compound_title = QLabel(compound_title_formal)									# QLabel is a widget that displays text
		self.compound_title.setAlignment(Qt.AlignCenter)									# Align the text to center
		self.compound_title_font = QFont("sans-serif", 24, QFont.Bold) 						# Declare font
		self.compound_title.setFont(self.compound_title_font)								# Set the font for the QLabel text
		self.chemicalpotential_phasediagram3d_window_layout.addWidget(self.compound_title)	# Add the widget to the "main" widget grid layout
		
		# (WIDGET) Chemical potential phase diagram plot object
		self.chemicalpotential_phase_diagram_plot = self.PhaseDiagram3D.chemicalpotential_phasediagram_plot_canvas
		self.chemicalpotential_phasediagram3d_window_layout.addWidget(self.chemicalpotential_phase_diagram_plot)
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		# (WIDGET) Fourth element slider (tunes delta mu of fourth element)
		self.fourth_element_slider_widget = QWidget()
		self.fourth_element_slider_layout = QHBoxLayout(self.fourth_element_slider_widget)
		
		self.fourth_element_selection_box = QComboBox()
		self.fourth_element_selection_box.addItem(self.first_element)
		self.fourth_element_selection_box.addItem(self.second_element)
		self.fourth_element_selection_box.addItem(self.third_element)
		self.fourth_element_selection_box.addItem(self.fourth_element)
		self.fourth_element_selection_box.setCurrentIndex(3)										# Set initial state of the drop-down menu to the fourth species
		self.fourth_element_selection_box.currentIndexChanged.connect(self.Update_Fourth_Species)	# When the user selects a species to be the fourth species, update the species list internally
		self.fourth_element_slider_layout.addWidget(self.fourth_element_selection_box)				# Add the drop-down menu widget to the slider widget
		
		self.fourth_element_slider_label = QLabel(u"\u0394"+"\u03BC"+"<sub>"+self.fourth_element+"</sub>")	# Create a label widget to display the "\Delta\mu_4" text next to the slider (unicode format)
		self.fourth_element_slider_label.setFont(QFont("sans-serif", 12))				# Set the font for the label
		self.fourth_element_slider_label.setAlignment(Qt.AlignCenter)					# Align the text to center
		self.fourth_element_slider_layout.addWidget(self.fourth_element_slider_label)	# Place the label widget in the main slider widget layout
		
		self.fourth_element_slider = QSlider(Qt.Horizontal)						# QSlider is the actual name of the widget for the slider object
		self.fourth_element_slider.setMinimum(0)								# Set minimum value
		self.fourth_element_slider.setMaximum(1000)								# Set maximum value
		self.fourth_element_slider.setValue(1000)								# Set the initial value
		self.fourth_element_slider.setSingleStep(1)								# Set the step value (how much one "slide" is worth)
		self.fourth_element_slider.setTickInterval(1)							# 
																				# *Note: The slider widget cannot handle floating values, i.e. you can only define
																				#	integers for the maximum, minimum, etc. The way we get around that is by
																				#	saying that the mu4 value is chosen from an array (self.mu4_value_array), and 
																				#	the integer value of the slider is the INDEX of the array.
		self.fourth_element_slider.setEnabled(False)
		self.fourth_element_slider.valueChanged.connect(self.Update_Fourth_Species_Slider)	# How to update the slider in the case that someone moves it
		self.fourth_element_slider_layout.addWidget(self.fourth_element_slider)	# Add the slider widget to the main slider widget
		
		self.fourth_element_slider_value_label = QLabel("{0:.4f}".format(-0.0))				# Create a label widget to display the current value of mu4 up to 4 digits
		self.fourth_element_slider_value_label.setAlignment(Qt.AlignCenter)					# Align the text to center
		self.fourth_element_slider_layout.addWidget(self.fourth_element_slider_value_label)	# Add this label widget to the main slider widget layout
		
		self.chemicalpotential_phasediagram3d_window_layout.addWidget(self.fourth_element_slider_widget)
		
		
		
		
		
		
		
		
		
		
		
		
		
		
		# (WIDGET) Save 3D phase diagram as figure
		self.phasediagram3d_savefigure_button = QPushButton("Save 3D Phase Diagram Figure")
		self.phasediagram3d_savefigure_button.clicked[bool].connect(lambda: self.SaveFigure("Phase Diagram 3D"))
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
		
		# (WIDGET) Three projected chemical potential phase diagrams
		self.phase_diagram_plot_2d_tripleview = self.PhaseDiagram2D_TripleView.tripleview_phase_diagram_plot_canvas
		
		# Add the three projected chemical potential phase diagrams to tab2
		self.tab2_layout.addWidget(self.phase_diagram_plot_2d_tripleview)
		
		
		self.Generate_PhaseDiagram3D_Function()
	
	
	
	
	"""
	def Update_Species(self, species_number):
		if species_number == 1:
			self.first_element = str(self.first_element_selection_box.currentText())
		elif species_number == 2:
			self.second_element = str(self.second_element_selection_box.currentText())
		elif species_number == 3:
			self.third_element = str(self.third_element_selection_box.currentText())
		elif species_number == 4:
			self.fourth_element = str(self.fourth_element_selection_box.currentText())
		self.elements_list = [self.first_element, self.second_element, self.third_element, self.fourth_element]
	"""
	def Update_Fourth_Species(self):
		
		# Change order of elements
		if str(self.fourth_element_selection_box.currentText()) == self.first_element:
			self.elements_list = [self.fourth_element, self.second_element, self.third_element, self.first_element]
			self.first_element = self.elements_list[0]
			self.fourth_element = self.elements_list[3]
		elif str(self.fourth_element_selection_box.currentText()) == self.second_element:
			self.elements_list = [self.first_element, self.fourth_element, self.third_element, self.second_element]
			self.second_element = self.elements_list[1]
			self.fourth_element = self.elements_list[3]
		elif str(self.fourth_element_selection_box.currentText()) == self.third_element:
			self.elements_list = [self.first_element, self.second_element, self.fourth_element, self.third_element]
			self.third_element = self.elements_list[2]
			self.fourth_element = self.elements_list[3]
		self.fourth_element_slider_label.setText(u"\u0394"+"\u03BC"+"<sub>"+self.fourth_element+"</sub>")
		
		# Update plots
		self.Generate_PhaseDiagram3D_Function()
	
	
	
	
	
	
	def Update_Fourth_Species_Slider(self):
		
		fourth_element_mu4_value_index = self.fourth_element_slider.value()
		fourth_element_mu4_value = self.mu4_value_array[fourth_element_mu4_value_index]
		
		self.PhaseDiagram3D.mu4 = fourth_element_mu4_value
		self.PhaseDiagram3D.Draw_Mu4_Outline()
		
		self.PhaseDiagram2D_TripleView.deltamu[4] = fourth_element_mu4_value
		self.PhaseDiagram2D_TripleView.Plot_PhaseDiagrams()
		
		
		mu4_rounded = round(fourth_element_mu4_value, 4)
		self.fourth_element_slider_value_label.setText("{0:.4f}".format(mu4_rounded))
	
	
	
	
	
	
	
	def Generate_PhaseDiagram3D_Function(self):
		
		# This function specifies what happens when the user clicks the "Generate Plot" button. The 
		#	only thing it needs to check is whether the elements chosen as the first, second, third, 
		#	and fourth species are unique.
		
		if len(self.elements_list) > len(set(self.elements_list)):	# Every time someone clicks the species' buttons, the self.elements_list gets updated.
																	#	The "set" function checks the self.elements_list and omits any that are repeated.
																	#	If any are repeated, then the chosen species are not unique.
			QMessageBox.about(self, "WARNING", "Pick UNIQUE elements!")
			return
		
		
		# Reset the slide bar
		self.fourth_element_slider.setEnabled(True)
		endpoint_slidebar = self.compounds_info[self.main_compound]["enthalpy"] / self.compounds_info[self.main_compound][self.fourth_element]
		self.mu4_value_array = np.linspace(endpoint_slidebar, -0.0, 1001)
		self.fourth_element_slider.setValue(1000)
		
		self.PhaseDiagram3D.chemicalpotential_phasediagram_plot_axes.clear()
		
		self.PhaseDiagram3D.Set_Elements(element_x = self.first_element, element_y = self.second_element, element_z = self.third_element, dependent_element = self.fourth_element)
		
		self.PhaseDiagram3D.Draw_PhaseDiagram3D()
		self.PhaseDiagram3D.Draw_Mu4_Outline()
		
		self.PhaseDiagram2D_TripleView.Set_Elements(first_element = self.first_element, second_element = self.second_element, third_element = self.third_element, fourth_element = self.fourth_element)
		self.PhaseDiagram2D_TripleView.Update_PhaseDiagram_Object()
		self.PhaseDiagram2D_TripleView.Update_PhaseDiagram_Plot_Axes()
		
		self.PhaseDiagram2D_TripleView.competing_compounds_colorwheel = self.PhaseDiagram3D.competing_compounds_colorwheel
		
		self.PhaseDiagram2D_TripleView.Plot_PhaseDiagrams()









