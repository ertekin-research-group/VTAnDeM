
__author__ = 'Michael_Lidia_Jiaxing_Benita_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'







import numpy as np

import PyQt5
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

from vtandem.visualization.ternary.ternary_scripts.plot_ternary_phase_diagram import ChemicalPotential_Ternary_PhaseDiagramProjected2D_TripleView
from vtandem.visualization.ternary.ternary_scripts.plot3d_ternary_phase_diagram import ChemicalPotential_Ternary_PhaseDiagram3D


class Tab_Ternary_PhaseDiagram3D(QWidget):
	
	def __init__(self, parent = None, main_compound = None, first_element = None, second_element = None, third_element = None, compounds_info = None):	# User specifies the main compound and its constituents
		
		QWidget.__init__(self)
		
		
		self.main_compound = main_compound
		self.first_element = first_element
		self.second_element = second_element
		self.third_element = third_element
		self.elements_list = [self.first_element, self.second_element, self.third_element]					# Species list (order MAY change)
		
		self.compounds_info = compounds_info
		
		
		self.tab2 = QWidget()
		self.tab2_layout = QHBoxLayout(self.tab2)
		
		
		
		
		###### "Settings" Widget
		self.phasediagram3d_widget = QWidget()									# One of the main sub-widgets is where the user defines the settings of the plots.
		self.phasediagram3d_widget_layout = QVBoxLayout(self.phasediagram3d_widget)		# The settings should be placed on top of one another, i.e. vertically.
		
		
		
		# (WIDGET) Title of compound
		ternary_compound_title_formal = self.Compound_Name_Formal(self.main_compound)	# Generate Latex-readable version of ternary compound name
		self.ternary_compound_title = QLabel(ternary_compound_title_formal)			# QLabel is a widget that displays text
		self.ternary_compound_title.setAlignment(Qt.AlignCenter)							# Align the text to center
		self.ternary_compound_title_font = QFont("sans-serif", 24, QFont.Bold) 				# Declare font
		self.ternary_compound_title.setFont(self.ternary_compound_title_font)			# Set the font for the QLabel text
		self.phasediagram3d_widget_layout.addWidget(self.ternary_compound_title)				# Add the widget to the "main" widget grid layout
		
		
		
		
		
		
		
		
		
		self.PhaseDiagram3D = ChemicalPotential_Ternary_PhaseDiagram3D(self, main_compound = self.main_compound, first_element = self.first_element, second_element = self.second_element, third_element = self.third_element)
		self.PhaseDiagram3D.compounds_info = self.compounds_info
		self.ternary_phase_diagram_plot_3d = self.PhaseDiagram3D.ternary_phasediagram_3d_plot_canvas
		self.phasediagram3d_widget_layout.addWidget(self.ternary_phase_diagram_plot_3d)
		
		
		
		
		
		
		
		
		self.elements_selector = QWidget()
		self.elements_selector_layout = QHBoxLayout(self.elements_selector)
		
		self.first_element_selection_box = QComboBox()
		self.first_element_selection_box.addItem(self.first_element)
		self.first_element_selection_box.addItem(self.second_element)
		self.first_element_selection_box.addItem(self.third_element)
		self.first_element_selection_box.setCurrentIndex(0)
		self.first_element_selection_box.activated.connect(lambda: self.Update_Species(1))
		self.elements_selector_layout.addWidget(self.first_element_selection_box)
		
		self.second_element_selection_box = QComboBox()
		self.second_element_selection_box.addItem(self.first_element)
		self.second_element_selection_box.addItem(self.second_element)
		self.second_element_selection_box.addItem(self.third_element)
		self.second_element_selection_box.setCurrentIndex(1)
		self.second_element_selection_box.activated.connect(lambda: self.Update_Species(2))
		self.elements_selector_layout.addWidget(self.second_element_selection_box)
		
		self.third_element_selection_box = QComboBox()
		self.third_element_selection_box.addItem(self.first_element)
		self.third_element_selection_box.addItem(self.second_element)
		self.third_element_selection_box.addItem(self.third_element)
		self.third_element_selection_box.setCurrentIndex(2)
		self.third_element_selection_box.activated.connect(lambda: self.Update_Species(3))
		self.elements_selector_layout.addWidget(self.third_element_selection_box)
		
		self.phasediagram3d_widget_layout.addWidget(self.elements_selector)
		
		
		
		
		
		# (WIDGET) Button to generate phase diagram
		self.generate_ternary_phasediagram3d_button_widget = QPushButton("Generate 3D Phase Diagram")
		self.generate_ternary_phasediagram3d_button_widget.clicked[bool].connect(self.Generate_PhaseDiagram3D_Function)
		self.phasediagram3d_widget_layout.addWidget(self.generate_ternary_phasediagram3d_button_widget)
		
		
		# (WIDGET) Save 3D phase diagram as figure
		self.phasediagram3d_savefigure_button = QPushButton("Save 3D Phase Diagram Figure")
		self.phasediagram3d_savefigure_button.clicked[bool].connect(lambda: self.SaveFigure("Phase Diagram 3D"))
		self.phasediagram3d_widget_layout.addWidget(self.phasediagram3d_savefigure_button)
		
		
		
		
		
		self.generate_phasediagram3d_animation_button = QPushButton("360 Rotation Animation")
		self.generate_phasediagram3d_animation_button.clicked[bool].connect(self.PhaseDiagram3D.Make_SpacePotato_Rotation_Animation)
		self.phasediagram3d_widget_layout.addWidget(self.generate_phasediagram3d_animation_button)
		
		
		
		
		
		
		self.tab2_layout.addWidget(self.phasediagram3d_widget)
		
		
		
		self.PhaseDiagram2D_TripleView = ChemicalPotential_Ternary_PhaseDiagramProjected2D_TripleView(self, main_compound = main_compound, first_element = self.second_element, second_element = self.third_element, third_element = self.first_element)
		self.PhaseDiagram2D_TripleView.compounds_info = self.compounds_info
		self.PhaseDiagram2D_TripleView.Update_PhaseDiagram_Object()
		self.PhaseDiagram2D_TripleView.Update_PhaseDiagram_Plot_Axes()
		self.ternary_phase_diagram_plot_2d_tripleview = self.PhaseDiagram2D_TripleView.ternary_phase_diagram_plot_canvas
		
		self.tab2_layout.addWidget(self.ternary_phase_diagram_plot_2d_tripleview)
	
	
	
	
	
	
	
	###############################################################################################
	########################## Rewrite Compound Name Latex-Style ##################################
	###############################################################################################
	
	def Compound_Name_Formal(self, compound_name):
		
		# Go into the compounds_info dictionary and obtains the chemistry and stoichiometry of the compound of choice
		compound_species_info = self.compounds_info[compound_name]
		
		compound_name_formal = ""
		for species in self.elements_list:				# Loop through the list of possible species that can be contained in the compound
			if compound_species_info[species] == 0:
				continue								# Don't add the species to the name if the compound doesn't contain the species
			elif compound_species_info[species] == 1:
				compound_name_formal += species			# Add the species to the name of the compound
			elif compound_species_info[species] > 1:
				compound_name_formal += species+"<sub>"+str(int(compound_species_info[species]))+"</sub>"	# Add the species to the name of the compound
																											#	with a subscript for the stoichiometry
		
		return compound_name_formal
	
	
	
	
	
	
	def Update_Species(self, species_number):
		
		if species_number == 1:
			self.first_element = str(self.first_element_selection_box.currentText())
		elif species_number == 2:
			self.second_element = str(self.second_element_selection_box.currentText())
		elif species_number == 3:
			self.third_element = str(self.third_element_selection_box.currentText())
		self.elements_list = [self.first_element, self.second_element, self.third_element]
	
	
	
	
	
	
	def Generate_PhaseDiagram3D_Function(self, event):
		
		# This function specifies what happens when the user clicks the "Generate Plot" button. The 
		#	only thing it needs to check is whether the elements chosen as the first, second, and 
		#	third species are unique.
		
		if len(self.elements_list) > len(set(self.elements_list)):	# Every time someone clicks the species' buttons, the self.elements_list gets updated.
																	#	The "set" function checks the self.elements_list and omits any that are repeated.
																	#	If any are repeated, then the chosen species are not unique.
			QMessageBox.about(self, "WARNING", "Pick UNIQUE elements!")
			return
		
		
		self.PhaseDiagram3D.ternary_phasediagram_3d_plot_axes.clear()
		
		self.PhaseDiagram3D.Set_Elements(element_x = self.first_element, element_y = self.second_element, element_z = self.third_element)
		
		self.PhaseDiagram3D.Draw_PhaseDiagram3D()
		
		self.PhaseDiagram2D_TripleView.Set_Elements(first_element = self.first_element, second_element = self.second_element, third_element = self.third_element)
		self.PhaseDiagram2D_TripleView.Update_PhaseDiagram_Object()
		self.PhaseDiagram2D_TripleView.Update_PhaseDiagram_Plot_Axes()
		
		self.PhaseDiagram2D_TripleView.competing_compounds_colorwheel = self.PhaseDiagram3D.competing_compounds_colorwheel
		
		self.PhaseDiagram2D_TripleView.Plot_PhaseDiagrams()
	
	
	
	###############################################################################################
	###################################### Save Figure ############################################
	###############################################################################################
	
	def SaveFigure(self, figure_type):
		
		options = QFileDialog.Options()
		options |= QFileDialog.DontUseNativeDialog
		filename, extension_type = QFileDialog.getSaveFileName(self, "Save "+figure_type+" Figure", "", "Portable Network Graphics (*.png);;" \
																										+"Portable Document Format (*.pdf);;" \
																										+"Scalable Vector Graphics (*.svg);;" \
																										+"Encapsulated PostScript (*.eps)", options=options)
		if filename:
			extension = extension_type.split(".")[-1].split(")")[0]
			if filename.split(".")[-1] == extension:
				if figure_type == "Phase Diagram 3D":
					self.PhaseDiagram3D.ternary_phasediagram_3d_plot_figure.savefig(filename, bbox_inches='tight')
			else:
				if figure_type == "Phase Diagram 3D":
					self.PhaseDiagram3D.ternary_phasediagram_3d_plot_figure.savefig(filename+"."+extension, bbox_inches='tight')






