
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'
__author__ = 'Michael_Lidia_Jiaxing_Elif'

import numpy as np

import PyQt5
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *


from vtandem.visualization.windows.window_defectsdiagram import Window_DefectsDiagram


class Window_DefectsDiagram_Binary(Window_DefectsDiagram):
	
	def __init__(self, show_dopant):
		
		Window_DefectsDiagram.__init__(self, show_dopant=show_dopant)


		# Chemical potential tuners
		self.chemical_potential_displays = QWidget()
		self.chemical_potential_displays_layout = QVBoxLayout(self.chemical_potential_displays)
		self.chemical_potential_first_element = QWidget()
		self.chemical_potential_first_element_layout = QHBoxLayout(self.chemical_potential_first_element)
		self.chemical_potential_first_element_label = QLabel(u"\u0394"+"\u03BC"+"<sub>"+self.first_element+"</sub> =")
		self.chemical_potential_first_element_label.setAlignment(Qt.AlignCenter)
		self.chemical_potential_first_element_layout.addWidget(self.chemical_potential_first_element_label)
		self.chemical_potential_first_element_deltamu = QLabel("-0.0000")
		self.chemical_potential_first_element_layout.addWidget(self.chemical_potential_first_element_deltamu)
		self.chemical_potential_displays_layout.addWidget(self.chemical_potential_first_element)
		self.chemical_potential_second_element = QWidget()
		self.chemical_potential_second_element_layout = QHBoxLayout(self.chemical_potential_second_element)
		self.chemical_potential_second_element_label = QLabel(u"\u0394"+"\u03BC"+"<sub>"+self.second_element+"</sub> =")
		self.chemical_potential_second_element_label.setAlignment(Qt.AlignCenter)
		self.chemical_potential_second_element_layout.addWidget(self.chemical_potential_second_element_label)
		self.chemical_potential_second_element_deltamu = QLabel('{:.4f}'.format(round(self.deltamu_values[self.second_element], 4)))
		self.chemical_potential_second_element_layout.addWidget(self.chemical_potential_second_element_deltamu)
		self.chemical_potential_displays_layout.addWidget(self.chemical_potential_second_element)
		
		# Chemical potential tuners
		self.chemical_potential_sliders = QWidget()
		self.chemical_potential_sliders_layout = QVBoxLayout(self.chemical_potential_sliders)
		self.chemical_potential_first_element_slider = QSlider(Qt.Horizontal)
		self.chemical_potential_first_element_slider.setMinimum(0)
		self.chemical_potential_first_element_slider.setMaximum(1000)
		self.chemical_potential_first_element_slider.setValue(1000)
		self.chemical_potential_first_element_slider.setSingleStep(1)
		self.chemical_potential_first_element_slider.setTickInterval(1)
		self.chemical_potential_first_element_slider.valueChanged.connect(lambda: self.Chemical_Potential_Slider(self.first_element))
		self.chemical_potential_sliders_layout.addWidget(self.chemical_potential_first_element_slider)
		self.chemical_potential_second_element_slider = QSlider(Qt.Horizontal)
		self.chemical_potential_second_element_slider.setMinimum(0)
		self.chemical_potential_second_element_slider.setMaximum(1000)
		self.chemical_potential_second_element_slider.setValue(0)
		self.chemical_potential_second_element_slider.setSingleStep(1)
		self.chemical_potential_second_element_slider.setTickInterval(1)
		self.chemical_potential_second_element_slider.valueChanged.connect(lambda: self.Chemical_Potential_Slider(self.second_element))
		self.chemical_potential_sliders_layout.addWidget(self.chemical_potential_second_element_slider)
		
		# Chemical potentials section
		self.chemical_potentials_section = QWidget()
		self.chemical_potentials_section_layout = QHBoxLayout(self.chemical_potentials_section)
		self.chemical_potentials_section_layout.addWidget(self.chemical_potential_displays)
		self.chemical_potentials_section_layout.addWidget(self.chemical_potential_sliders)
		self.defectsdiagram_window_layout.insertWidget(3, self.chemical_potentials_section)
		

		
		# (WIDGET) Button to generate defects diagram
		self.generate_defects_diagram_plot_button_widget = QPushButton("Generate Defects Diagram")
		self.generate_defects_diagram_plot_button_widget.clicked[bool].connect(self.Generate_DefectsDiagram_Plot_Function_Binary)
		self.defectsdiagram_window_layout.insertWidget(5, self.generate_defects_diagram_plot_button_widget)
		
		


