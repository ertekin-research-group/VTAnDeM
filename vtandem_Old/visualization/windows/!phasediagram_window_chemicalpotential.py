
__author__ = 'Michael_Lidia_Jiaxing_Benita_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'


import numpy as np

import PyQt5
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

from vtandem.visualization.utils.compound_name import Compound_Name_Formal



class ChemicalPotential_PhaseDiagram_Window(QWidget):
	
	def __init__(self, main_compound: str, PhaseDiagram3D, PhaseDiagram2D_TripleView):
		
		super().__init__()
		
		# Initialize 3D phase diagram object
		self.PhaseDiagram3D = PhaseDiagram3D
		
		###### Main chemical potential phase diagram window widget
		self.chemicalpotential_phasediagram3d_window = QWidget()															# One of the main sub-widgets is where the user defines the settings of the plots.
		self.chemicalpotential_phasediagram3d_window_layout = QVBoxLayout(self.chemicalpotential_phasediagram3d_window)		# The settings should be placed on top of one another, i.e. vertically.
		
		# (WIDGET) Title of compound
		compound_title_formal = Compound_Name_Formal(main_compound, self.PhaseDiagram3D.compounds_info, "unicode")	# Generate Latex-readable version of compound name
		self.compound_title = QLabel(compound_title_formal)					# QLabel is a widget that displays text
		self.compound_title.setAlignment(Qt.AlignCenter)					# Align the text to center
		self.compound_title_font = QFont("sans-serif", 24, QFont.Bold) 		# Declare font
		self.compound_title.setFont(self.compound_title_font)				# Set the font for the QLabel text
		#self.phasediagram3d_widget_layout.addWidget(self.compound_title)	# Add the widget to the "main" widget grid layout
		
		# (WIDGET) Chemical potential phase diagram plot object
		self.chemicalpotential_phase_diagram_plot = self.PhaseDiagram3D.chemicalpotential_phasediagram_plot_canvas
		#self.phasediagram3d_widget_layout.addWidget(self.chemicalpotential_phase_diagram_plot)
		
		
		# (WIDGET) Button to generate phase diagram
		self.generate_phasediagram3d_button_widget = QPushButton("Generate 3D Phase Diagram")
		self.generate_phasediagram3d_button_widget.clicked[bool].connect(self.Generate_PhaseDiagram3D_Function)
		#self.phasediagram3d_widget_layout.addWidget(self.generate_phasediagram3d_button_widget)
		
		# (WIDGET) Save 3D phase diagram as figure
		self.phasediagram3d_savefigure_button = QPushButton("Save 3D Phase Diagram Figure")
		self.phasediagram3d_savefigure_button.clicked[bool].connect(lambda: self.SaveFigure("Phase Diagram 3D"))
		#self.phasediagram3d_widget_layout.addWidget(self.phasediagram3d_savefigure_button)
		
		# (WIDGET) Generate a 3D rotating phase diagram animation
		self.generate_phasediagram3d_animation_button = QPushButton("360 Rotation Animation")
		self.generate_phasediagram3d_animation_button.clicked[bool].connect(self.PhaseDiagram3D.Make_SpacePotato_Rotation_Animation)
		#self.phasediagram3d_widget_layout.addWidget(self.generate_phasediagram3d_animation_button)
		
		# (WIDGET) Three projected chemical potential phase diagrams
		self.PhaseDiagram2D_TripleView = PhaseDiagram2D_TripleView
		self.phase_diagram_plot_2d_tripleview = self.PhaseDiagram2D_TripleView.tripleview_phase_diagram_plot_canvas
	
	
	
	
	
	
	
	"""
	###############################################################################################
	########################## Rewrite Compound Name Latex-Style ##################################
	###############################################################################################
	
	def Compound_Name_Formal(self, compound_name):
		
		# Go into the compounds_info dictionary and obtains the chemistry and stoichiometry of the compound of choice
		compound_species_info = self.PhaseDiagram3D.compounds_info[compound_name]
		
		compound_name_formal = ""
		for species in self.PhaseDiagram3D.elements_list:				# Loop through the list of possible species that can be contained in the compound
			if compound_species_info[species] == 0:
				continue								# Don't add the species to the name if the compound doesn't contain the species
			elif compound_species_info[species] == 1:
				compound_name_formal += species			# Add the species to the name of the compound
			elif compound_species_info[species] > 1:
				compound_name_formal += species+"<sub>"+str(int(compound_species_info[species]))+"</sub>"	# Add the species to the name of the compound
																											#	with a subscript for the stoichiometry
		
		return compound_name_formal
	"""
	
	
	
	
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
	
	
	
	
	
	
	def Update_Fourth_Species_Slider(self):
		
		fourth_element_mu4_value_index = self.fourth_element_slider.value()
		fourth_element_mu4_value = self.mu4_value_array[fourth_element_mu4_value_index]
		
		self.PhaseDiagram3D.mu4 = fourth_element_mu4_value
		self.PhaseDiagram3D.Draw_Mu4_Outline()
		
		self.PhaseDiagram2D_TripleView.mu4 = fourth_element_mu4_value
		self.PhaseDiagram2D_TripleView.Plot_PhaseDiagrams()
		
		
		mu4_rounded = round(fourth_element_mu4_value, 4)
		self.fourth_element_slider_value_label.setText("{0:.4f}".format(mu4_rounded))
	"""
	
	
	
	
	
	
	"""
	def Generate_PhaseDiagram3D_Function(self, event):
		
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
		self.fourth_element_slider_label.setText(u"\u0394"+"\u03BC"+"<sub>"+self.fourth_element+"</sub>")
		endpoint_slidebar = self.compounds_info[self.main_compound]["enthalpy"] / self.compounds_info[self.main_compound][self.fourth_element]
		self.mu4_value_array = np.linspace(endpoint_slidebar, -0.0, 1001)
		self.fourth_element_slider.setValue(1000)
		
		self.PhaseDiagram3D.quaternary_phasediagram_3d_plot_axes.clear()
		
		self.PhaseDiagram3D.Set_Elements(element_x = self.first_element, element_y = self.second_element, element_z = self.third_element, dependent_element = self.fourth_element)
		
		self.PhaseDiagram3D.Draw_PhaseDiagram3D()
		self.PhaseDiagram3D.Draw_Mu4_Outline()
		
		self.PhaseDiagram2D_TripleView.Set_Elements(first_element = self.first_element, second_element = self.second_element, third_element = self.third_element, fourth_element = self.fourth_element)
		self.PhaseDiagram2D_TripleView.Update_PhaseDiagram_Object()
		self.PhaseDiagram2D_TripleView.Update_PhaseDiagram_Plot_Axes()
		
		self.PhaseDiagram2D_TripleView.competing_compounds_colorwheel = self.PhaseDiagram3D.competing_compounds_colorwheel
		
		self.PhaseDiagram2D_TripleView.Plot_PhaseDiagrams()
	"""
	
	
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
					self.PhaseDiagram3D.chemicalpotential_phasediagram_plot_figure.savefig(filename, bbox_inches='tight')
			else:
				if figure_type == "Phase Diagram 3D":
					self.PhaseDiagram3D.chemicalpotential_phasediagram_plot_figure.savefig(filename+"."+extension, bbox_inches='tight')






