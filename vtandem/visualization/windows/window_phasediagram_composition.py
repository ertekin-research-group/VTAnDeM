
__author__ = 'Michael_Lidia_Jiaxing_Benita_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'


import numpy as np

import PyQt5
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

from vtandem.visualization.utils.compound_name import Compound_Name_Formal


class Window_Compositional_PhaseDiagram(QWidget):
	
	def __init__(self, main_compound: str, Compositional_PhaseDiagram_Object):
		
		super().__init__()
		
		# Initialize compositional phase diagram object
		self.Compositional_PhaseDiagram = Compositional_PhaseDiagram_Object
		
		###### Main compositional phase diagram window widget
		self.compositional_phasediagram_window = QWidget()									# One of the main sub-widgets is where the user defines the settings of the plots.
		self.compositional_phasediagram_window_layout = QVBoxLayout(self.compositional_phasediagram_window)		# The settings should be placed on top of one another, i.e. vertically.
		
		"""
		# (WIDGET) Title of compound
		compound_title_formal = Compound_Name_Formal(main_compound, self.Compositional_PhaseDiagram.compounds_info, "unicode")	# Generate Latex-readable version of compound name
		self.compound_title = QLabel(compound_title_formal)					# QLabel is a widget that displays text
		self.compound_title.setAlignment(Qt.AlignCenter)					# Align the text to center
		self.compound_title_font = QFont("sans-serif", 24, QFont.Bold) 		# Declare font
		self.compound_title.setFont(self.compound_title_font)				# Set the font for the QLabel text
		"""
		# (WIDGET) Title
		self.composition_phase_diagram_title = QLabel("Phase \n Diagram \n (Composition)")			# QLabel is a widget that displays text
		self.composition_phase_diagram_title.setAlignment(Qt.AlignCenter)						# Align the text to center
		self.composition_phase_diagram_title_font = QFont("sans-serif", 24, QFont.Bold) 		# Declare font
		self.composition_phase_diagram_title.setFont(self.composition_phase_diagram_title_font)	# Set the font for the QLabel text
		
		# (WIDGET) Compositional phase diagram plot object
		self.composition_phase_diagram_plot = self.Compositional_PhaseDiagram.composition_phasediagram_plot_canvas
		
		# (WIDGET) Save phase diagram as figure
		self.phasediagram_savefigure_button = QPushButton("Save Phase Diagram Figure")
		self.phasediagram_savefigure_button.clicked[bool].connect(lambda: self.Compositional_PhaseDiagram.SaveFigure())















