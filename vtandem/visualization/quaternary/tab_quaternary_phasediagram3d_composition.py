
__author__ = 'Michael_Lidia_Jiaxing_Benita_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'




import numpy as np

import PyQt5
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

from vtandem.visualization.quaternary.quaternary_scripts.plot3d_composition_quaternary_phase_diagram import Composition_Quaternary_PhaseDiagram3D


class Tab_Compositional_PhaseDiagram3D(QWidget):
	
	def __init__(self, parent = None, main_compound = None, first_element = None, second_element = None, third_element = None, fourth_element = None, compounds_info = None):	# User specifies the main compound and its constituents
		
		QWidget.__init__(self)
		
		
		self.main_compound = main_compound
		self.first_element = first_element
		self.second_element = second_element
		self.third_element = third_element
		self.fourth_element = fourth_element
		self.elements_list = [self.first_element, self.second_element, self.third_element, self.fourth_element]					# Species list (order MAY change)
		
		self.compounds_info = compounds_info
		
		
		self.tab3 = QWidget()
		self.tab3_layout = QHBoxLayout(self.tab3)
		
		
		
		
		###### "Settings" Widget
		self.settings_widget = QWidget()									# One of the main sub-widgets is where the user defines the settings of the plots.
		self.settings_widget_layout = QVBoxLayout(self.settings_widget)		# The settings should be placed on top of one another, i.e. vertically.
		
		
		
		# (WIDGET) Title of compound
		quaternary_compound_title_formal = self.Compound_Name_Formal(self.main_compound)	# Generate Latex-readable version of quaternary compound name
		self.quaternary_compound_title = QLabel(quaternary_compound_title_formal)			# QLabel is a widget that displays text
		self.quaternary_compound_title.setAlignment(Qt.AlignCenter)							# Align the text to center
		self.quaternary_compound_title_font = QFont("sans-serif", 24, QFont.Bold) 				# Declare font
		self.quaternary_compound_title.setFont(self.quaternary_compound_title_font)			# Set the font for the QLabel text
		self.settings_widget_layout.addWidget(self.quaternary_compound_title)				# Add the widget to the "main" widget grid layout
		
		
		
		
		
		# (WIDGET) Button to generate phase diagram
		self.generate_quaternary_phasediagram3d_button_widget = QPushButton("Generate 3D Phase Diagram")
		self.generate_quaternary_phasediagram3d_button_widget.clicked[bool].connect(self.Generate_Compositional_PhaseDiagram3D_Function)
		self.settings_widget_layout.addWidget(self.generate_quaternary_phasediagram3d_button_widget)
		
		
		# (WIDGET) Save 3D phase diagram as figure
		self.phasediagram3d_savefigure_button = QPushButton("Save 3D Phase Diagram Figure")
		self.phasediagram3d_savefigure_button.clicked[bool].connect(lambda: self.SaveFigure("Phase Diagram 3D Composition"))
		self.settings_widget_layout.addWidget(self.phasediagram3d_savefigure_button)
		
		
		self.tab3_layout.addWidget(self.settings_widget)
		
		
		
		
		
		self.Compositional_PhaseDiagram3D = Composition_Quaternary_PhaseDiagram3D(self, main_compound = self.main_compound, first_element = self.first_element, second_element = self.second_element, third_element = self.third_element, fourth_element = self.fourth_element)
		self.Compositional_PhaseDiagram3D.compounds_info = self.compounds_info
		
		
		
		
		self.composition_phase_diagram_plot_3d = self.Compositional_PhaseDiagram3D.composition_phasediagram_3d_plot_canvas
		
		
		
		
		
		
		self.tab3_layout.addWidget(self.composition_phase_diagram_plot_3d)
	
	
	
	
	
	
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
	
	
	
	def Generate_Compositional_PhaseDiagram3D_Function(self, event):
		self.Compositional_PhaseDiagram3D.Create_Compositional_PhaseDiagram3D()
		self.Compositional_PhaseDiagram3D.Plot_Compositional_PhaseDiagram3D()
	
	
	
	
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
				if figure_type == "Phase Diagram 3D Composition":
					self.PhaseDiagram3D.quaternary_phasediagram_3d_plot_figure.savefig(filename, bbox_inches='tight')
			else:
				if figure_type == "Phase Diagram 3D Composition":
					self.PhaseDiagram3D.quaternary_phasediagram_3d_plot_figure.savefig(filename+"."+extension, bbox_inches='tight')


















