
__author__ = 'Michael_Lidia_Jiaxing_Benita_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'


import numpy as np

import PyQt5
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *



class DefectsDiagram_Window(QWidget):
	
	def __init__(self, main_compound: str, DefectsDiagram_Object):
		
		super().__init__()
		
		# Initialize defects diagram object
		self.DefectsDiagram = DefectsDiagram_Object
		self.DefectsDiagram.Activate_DefectsDiagram_Plot_Axes()
		
		###### Main defects diagram window widget
		self.defectsdiagram_window = QWidget()									# One of the main sub-widgets is where the user defines the settings of the plots.
		self.defectsdiagram_window_layout = QVBoxLayout(self.defectsdiagram_window)		# The settings should be placed on top of one another, i.e. vertically.
		
		# Defects diagram plot
		self.defects_diagram_plot = self.DefectsDiagram.defects_diagram_plot_canvas
		self.defectsdiagram_window_layout.addWidget(self.defects_diagram_plot)
		
		# Y-axis limits for defects diagram
		self.defectsdiagram_viewport = QWidget()
		self.defectsdiagram_viewport_layout = QHBoxLayout(self.defectsdiagram_viewport)
		
		# Y-axis limits for defects diagram
		self.defectsdiagram_Ymin_label = QLabel(u"y"+"<sub>min</sub>")
		self.defectsdiagram_Ymin_label.setAlignment(Qt.AlignRight)
		self.defectsdiagram_viewport_layout.addWidget(self.defectsdiagram_Ymin_label)
		self.defectsdiagram_Ymin_box = QLineEdit("-2.0")
		self.defectsdiagram_Ymin_box.editingFinished.connect(lambda: self.Update_WindowSize("DefectsDiagram", "YMin"))
		self.defectsdiagram_viewport_layout.addWidget(self.defectsdiagram_Ymin_box)
		self.defectsdiagram_Ymax_label = QLabel(u"y"+"<sub>max</sub>")
		self.defectsdiagram_Ymax_label.setAlignment(Qt.AlignRight)
		self.defectsdiagram_viewport_layout.addWidget(self.defectsdiagram_Ymax_label)
		self.defectsdiagram_Ymax_box = QLineEdit("2.0")
		self.defectsdiagram_Ymax_box.editingFinished.connect(lambda: self.Update_WindowSize("DefectsDiagram", "YMax"))
		self.defectsdiagram_viewport_layout.addWidget(self.defectsdiagram_Ymax_box)
		self.defectsdiagram_window_layout.addWidget(self.defectsdiagram_viewport)
		
		# (WIDGET) Button to generate defects diagram
		self.generate_defects_diagram_plot_button_widget = QPushButton("Generate Defects Diagram")
		self.generate_defects_diagram_plot_button_widget.clicked[bool].connect(self.Generate_DefectsDiagram_Plot_Function)
		self.defectsdiagram_window_layout.addWidget(self.generate_defects_diagram_plot_button_widget)
		
		# (WIDGET) Save defects diagram as figure
		self.defects_diagram_savefigure_button = QPushButton("Save Defects Diagram Figure")
		self.defects_diagram_savefigure_button.clicked[bool].connect(lambda: self.SaveFigure("Defects Diagram"))
		self.defectsdiagram_window_layout.addWidget(self.defects_diagram_savefigure_button)
	
	
	
	###############################################################################################
	################################# Generate Defects Diagram ####################################
	###############################################################################################
	
	def Generate_DefectsDiagram_Plot_Function(self, event):
		
		# This function specifies what happens when the user clicks the "Generate Defects Diagram" button.
		
		# Update elements and chemical potentials
		self.DefectsDiagram.mu_elements[self.first_element]["deltamu"] = self.deltamu_values[self.first_element]
		self.DefectsDiagram.mu_elements[self.second_element]["deltamu"] = self.deltamu_values[self.second_element]
		self.DefectsDiagram.mu_elements[self.third_element]["deltamu"] = self.deltamu_values[self.third_element]
		
		# Reset defects diagram
		self.DefectsDiagram.defects_diagram_plot_drawing.remove()
		self.DefectsDiagram.defects_diagram_plot_drawing = self.DefectsDiagram.defects_diagram_plot_figure.add_subplot(111)
		self.DefectsDiagram.Activate_DefectsDiagram_Plot_Axes()
		
		# Calculate defect formation energies
		self.DefectsDiagram.Calculate_DefectFormations()
		
		# Plot defect formation energies
		self.DefectsDiagram.intrinsic_defect_plots = {}
		self.DefectsDiagram.extrinsic_defect_plots = {}
		self.DefectsDiagram.Initialize_Intrinsic_DefectsDiagram_Plot()
	
	
	
	
	
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
				if figure_type == "Defects Diagram":
					self.Compositional_PhaseDiagram.composition_phasediagram_plot_figure.savefig(filename, bbox_inches='tight')
			else:
				if figure_type == "Defects Diagram":
					self.Compositional_PhaseDiagram.composition_phasediagram_plot_figure.savefig(filename+"."+extension, bbox_inches='tight')



