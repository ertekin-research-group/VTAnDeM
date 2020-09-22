
__author__ = 'Michael_Lidia_Jiaxing_Benita_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'


import numpy as np

import PyQt5
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

from vtandem.visualization.compound_name import Compound_Name_Formal


class Compositional_PhaseDiagram_Window(QWidget):
	
	def __init__(self, main_compound: str, Compositional_PhaseDiagram_Object):
		
		super().__init__()
		
		# Initialize compositional phase diagram object
		self.Compositional_PhaseDiagram = Compositional_PhaseDiagram_Object
		
		###### Main compositional phase diagram window widget
		self.compositional_phasediagram_window = QWidget()									# One of the main sub-widgets is where the user defines the settings of the plots.
		self.compositional_phasediagram_window_layout = QVBoxLayout(self.compositional_phasediagram_window)		# The settings should be placed on top of one another, i.e. vertically.
		
		# (WIDGET) Title of compound
		compound_title_formal = Compound_Name_Formal(main_compound, self.Compositional_PhaseDiagram.compounds_info, "unicode")	# Generate Latex-readable version of compound name
		self.compound_title = QLabel(compound_title_formal)					# QLabel is a widget that displays text
		self.compound_title.setAlignment(Qt.AlignCenter)					# Align the text to center
		self.compound_title_font = QFont("sans-serif", 24, QFont.Bold) 		# Declare font
		self.compound_title.setFont(self.compound_title_font)				# Set the font for the QLabel text
		#self.compositional_phasediagram_window_layout.addWidget(self.compound_title)			# Add the widget to the "main" widget grid layout
		
		# (WIDGET) Compositional phase diagram plot object
		self.composition_phase_diagram_plot = self.Compositional_PhaseDiagram.composition_phasediagram_plot_canvas
		#self.compositional_phasediagram_window_layout.addWidget(self.composition_phase_diagram_plot)
		
		# (WIDGET) Save phase diagram as figure
		self.phasediagram_savefigure_button = QPushButton("Save Phase Diagram Figure")
		self.phasediagram_savefigure_button.clicked[bool].connect(lambda: self.SaveFigure_Function())
		#self.compositional_phasediagram_window_layout.addWidget(self.phasediagram_savefigure_button)
	
	
	
	
	###############################################################################################
	###################################### Save Figure ############################################
	###############################################################################################
	
	def SaveFigure_Function(self):
		
		options = QFileDialog.Options()
		options |= QFileDialog.DontUseNativeDialog
		filename, extension_type = QFileDialog.getSaveFileName(filter = "Portable Network Graphics (*.png);;" \
																+"Portable Document Format (*.pdf);;" \
																+"Scalable Vector Graphics (*.svg);;" \
																+"Encapsulated PostScript (*.eps)", options=options)
		self.Compositional_PhaseDiagram.SaveFigure(filename, extension_type)















