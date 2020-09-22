
__author__ = 'Michael_Lidia_Jiaxing_Benita_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'



import numpy as np

import PyQt5
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *


class ChemicalPotential_PhaseDiagram_Window(QWidget):
	
	def __init__(self, main_compound: str, PhaseDiagram3D, PhaseDiagram2D_TripleView, type: str):
		
		super().__init__()
		
		# Record type of phase diagram (either ternary or quaternary)
		if type not in ["ternary", "quaternary"]:
			raise ValueError("Argument 'type' can be either 'ternary' or 'quaternary'.")
		self.type = type
		
		
		
		self.tab2 = QWidget()
		self.tab2_layout = QHBoxLayout(self.tab2)
		
		
		
		
		###### "Settings" Widget
		self.phasediagram3d_widget = QWidget()									# One of the main sub-widgets is where the user defines the settings of the plots.
		self.phasediagram3d_widget_layout = QVBoxLayout(self.phasediagram3d_widget)		# The settings should be placed on top of one another, i.e. vertically.
		
		
		
		# (WIDGET) Title of compound
		compound_title_formal = self.Compound_Name_Formal(main_compound)	# Generate Latex-readable version of compound name
		self.compound_title = QLabel(compound_title_formal)					# QLabel is a widget that displays text
		self.compound_title.setAlignment(Qt.AlignCenter)					# Align the text to center
		self.compound_title_font = QFont("sans-serif", 24, QFont.Bold) 		# Declare font
		self.compound_title.setFont(self.compound_title_font)				# Set the font for the QLabel text
		self.phasediagram3d_widget_layout.addWidget(self.compound_title)	# Add the widget to the "main" widget grid layout
		
		
		
		
		
		
		
		
		
		self.PhaseDiagram3D = PhaseDiagram3D
		self.chemicalpotential_phase_diagram_plot = self.PhaseDiagram3D.chemicalpotential_phasediagram_plot_canvas
		self.phasediagram3d_widget_layout.addWidget(self.chemicalpotential_phase_diagram_plot)
		
		
		
		
		
		
		
		self.elements_selector = QWidget()
		self.elements_selector_layout = QHBoxLayout(self.elements_selector)
		self.first_element_selection_box = QComboBox()
		self.elements_selector_layout.addWidget(self.first_element_selection_box)
		self.second_element_selection_box = QComboBox()
		self.elements_selector_layout.addWidget(self.second_element_selection_box)
		self.third_element_selection_box = QComboBox()
		self.elements_selector_layout.addWidget(self.third_element_selection_box)
		self.phasediagram3d_widget_layout.addWidget(self.elements_selector)
		
		
		
		if self.type == "quaternary":
			
			# Create the slider widget for the chemical potential (mu) of the fourth species.
			self.fourth_element_slider_widget = QWidget()
			self.fourth_element_slider_layout = QHBoxLayout(self.fourth_element_slider_widget)
			
			self.fourth_element_selection_box = QComboBox()
			self.fourth_element_slider_layout.addWidget(self.fourth_element_selection_box)		# Add the drop-down menu widget to the slider widget
			
			self.fourth_element_slider_label = QLabel(u"\u0394"+"\u03BC"+"<sub>d</sub>")	# Create a label widget to display the "\Delta\mu_4" text next to the slider (unicode format)
			self.fourth_element_slider_label.setFont(QFont("sans-serif", 12))			# Set the font for the label
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
			self.fourth_element_slider_layout.addWidget(self.fourth_element_slider)	# Add the slider widget to the main slider widget
			
			
			self.fourth_element_slider_value_label = QLabel("{0:.4f}".format(-0.0))				# Create a label widget to display the current value of mu4 up to 4 digits
			self.fourth_element_slider_value_label.setAlignment(Qt.AlignCenter)					# Align the text to center
			self.fourth_element_slider_layout.addWidget(self.fourth_element_slider_value_label)	# Add this label widget to the main slider widget layout
			
			self.phasediagram3d_widget_layout.addWidget(self.fourth_element_slider_widget)
		
		
		
		
		# (WIDGET) Button to generate phase diagram
		self.generate_phasediagram3d_button_widget = QPushButton("Generate 3D Phase Diagram")
		self.generate_phasediagram3d_button_widget.clicked[bool].connect(self.Generate_PhaseDiagram3D_Function)
		self.phasediagram3d_widget_layout.addWidget(self.generate_phasediagram3d_button_widget)
		
		
		# (WIDGET) Save 3D phase diagram as figure
		self.phasediagram3d_savefigure_button = QPushButton("Save 3D Phase Diagram Figure")
		self.phasediagram3d_savefigure_button.clicked[bool].connect(lambda: self.SaveFigure("Phase Diagram 3D"))
		self.phasediagram3d_widget_layout.addWidget(self.phasediagram3d_savefigure_button)
		
		
		
		
		
		self.generate_phasediagram3d_animation_button = QPushButton("360 Rotation Animation")
		self.generate_phasediagram3d_animation_button.clicked[bool].connect(self.PhaseDiagram3D.Make_SpacePotato_Rotation_Animation)
		self.phasediagram3d_widget_layout.addWidget(self.generate_phasediagram3d_animation_button)
		
		
		
		
		
		
		self.tab2_layout.addWidget(self.phasediagram3d_widget)
		
		
		
		self.PhaseDiagram2D_TripleView = PhaseDiagram2D_TripleView
		self.phase_diagram_plot_2d_tripleview = self.PhaseDiagram2D_TripleView.tripleview_phase_diagram_plot_canvas
		self.tab2_layout.addWidget(self.phase_diagram_plot_2d_tripleview)
	
	
	
	
	
	
	
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






