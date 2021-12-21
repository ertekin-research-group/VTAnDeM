
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'
__author__ = 'Michael_Lidia_Jiaxing_Elif'


import numpy as np

import PyQt5
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *



class Window_DefectsDiagram:
	
	def __init__(self, show_dopant):
		
		self.show_dopant = show_dopant
		
		
		###### Main defects diagram window widget
		self.defectsdiagram_window = QWidget()										# One of the main sub-widgets is where the user defines the settings of the plots.
		self.defectsdiagram_window_layout = QVBoxLayout(self.defectsdiagram_window)	# The settings should be placed on top of one another, i.e. vertically.
		
		# Title
		self.defectsdiagram_title = QLabel("Defect \n Energies")				# QLabel is a widget that displays text
		self.defectsdiagram_title.setAlignment(Qt.AlignCenter)					# Align the text to center
		self.defectsdiagram_title_font = QFont("sans-serif", 24, QFont.Bold) 	# Declare font
		self.defectsdiagram_title.setFont(self.defectsdiagram_title_font)		# Set the font for the QLabel text
		self.defectsdiagram_window_layout.addWidget(self.defectsdiagram_title)
		
		# Defects diagram plot
		self.defects_diagram_plot = self.DefectsDiagram.defects_diagram_plot_canvas
		self.defectsdiagram_window_layout.addWidget(self.defects_diagram_plot)
		
		# Axis limits for defects diagram
		self.defectsdiagram_viewport = QWidget()
		self.defectsdiagram_viewport_layout = QHBoxLayout(self.defectsdiagram_viewport)
		self.defectsdiagram_axislim_boxes = {}

		# X-axis limits for defects diagram
		defectsdiagram_Xmin_label = QLabel(u"x"+"<sub>min</sub>")
		defectsdiagram_Xmin_label.setAlignment(Qt.AlignCenter)
		self.defectsdiagram_viewport_layout.addWidget(defectsdiagram_Xmin_label)
		self.defectsdiagram_Xmin_box = QLineEdit("0.0")
		self.defectsdiagram_axislim_boxes["XMin"] = self.defectsdiagram_Xmin_box
		self.defectsdiagram_Xmin_box.editingFinished.connect(lambda: self.DefectsDiagram.Update_WindowSize("XMin", self.defectsdiagram_axislim_boxes))
		self.defectsdiagram_viewport_layout.addWidget(self.defectsdiagram_Xmin_box)
		defectsdiagram_Xmax_label = QLabel(u"x"+"<sub>max</sub>")
		defectsdiagram_Xmax_label.setAlignment(Qt.AlignCenter)
		self.defectsdiagram_viewport_layout.addWidget(defectsdiagram_Xmax_label)
		self.defectsdiagram_Xmax_box = QLineEdit(str(round(self.DefectsDiagram.ECBM-self.DefectsDiagram.EVBM,4)))
		self.defectsdiagram_axislim_boxes["XMax"] = self.defectsdiagram_Xmax_box
		self.defectsdiagram_Xmax_box.editingFinished.connect(lambda: self.DefectsDiagram.Update_WindowSize("XMax", self.defectsdiagram_axislim_boxes))
		self.defectsdiagram_viewport_layout.addWidget(self.defectsdiagram_Xmax_box)

		# Y-axis limits for defects diagram
		defectsdiagram_Ymin_label = QLabel(u"y"+"<sub>min</sub>")
		defectsdiagram_Ymin_label.setAlignment(Qt.AlignCenter)
		self.defectsdiagram_viewport_layout.addWidget(defectsdiagram_Ymin_label)
		self.defectsdiagram_Ymin_box = QLineEdit("-2.0")
		self.defectsdiagram_axislim_boxes["YMin"] = self.defectsdiagram_Ymin_box
		self.defectsdiagram_Ymin_box.editingFinished.connect(lambda: self.DefectsDiagram.Update_WindowSize("YMin", self.defectsdiagram_axislim_boxes))
		self.defectsdiagram_viewport_layout.addWidget(self.defectsdiagram_Ymin_box)
		defectsdiagram_Ymax_label = QLabel(u"y"+"<sub>max</sub>")
		defectsdiagram_Ymax_label.setAlignment(Qt.AlignCenter)
		self.defectsdiagram_viewport_layout.addWidget(defectsdiagram_Ymax_label)
		self.defectsdiagram_Ymax_box = QLineEdit("2.0")
		self.defectsdiagram_axislim_boxes["YMax"] = self.defectsdiagram_Ymax_box
		self.defectsdiagram_Ymax_box.editingFinished.connect(lambda: self.DefectsDiagram.Update_WindowSize("YMax", self.defectsdiagram_axislim_boxes))
		self.defectsdiagram_viewport_layout.addWidget(self.defectsdiagram_Ymax_box)
		self.defectsdiagram_window_layout.addWidget(self.defectsdiagram_viewport)
		

		# Synthesis temperature
		self.defects_synthesis_temperature = QWidget()
		self.defects_synthesis_temperature_layout = QHBoxLayout(self.defects_synthesis_temperature)
		defects_synthesis_temperature_label = QLabel(u"T<sub>syn</sub> (K) = ")
		defects_synthesis_temperature_label.setAlignment(Qt.AlignCenter)
		self.defects_synthesis_temperature_layout.addWidget(defects_synthesis_temperature_label)
		self.defects_synthesis_temperature_box = QLineEdit("")
		self.defects_synthesis_temperature_box.editingFinished.connect(self.Update_SynthesisTemperature)
		self.defects_synthesis_temperature_box.setEnabled(False)
		self.defects_synthesis_temperature_layout.addWidget(self.defects_synthesis_temperature_box)

		
		if self.show_dopant:
			
			# Extrinsic defect properties
			self.dopant_properties_widget = QWidget()
			self.dopant_properties_widget_layout = QHBoxLayout(self.dopant_properties_widget)
			
			# Extrinsic defect chemical potential
			self.dopant_chemical_potential_label = QLabel(u"\u0394"+"\u03BC"+"<sub>x</sub>")
			self.dopant_chemical_potential_label.setAlignment(Qt.AlignCenter)
			self.dopant_properties_widget_layout.addWidget(self.dopant_chemical_potential_label)
			self.dopant_chemical_potential_deltamu = QLineEdit("-0.0000")
			self.dopant_chemical_potential_deltamu.setMaxLength(7)
			self.dopant_chemical_potential_deltamu.editingFinished.connect(self.Update_ExtrinsicDefect_DeltaMu)
			self.dopant_chemical_potential_deltamu.setEnabled(False)
			self.dopant_properties_widget_layout.addWidget(self.dopant_chemical_potential_deltamu)
			
			# Extrinsic defect selection box
			self.dopant_selection_box = QComboBox()
			self.dopant_selection_box.setEnabled(False)
			self.dopant_selection_box.addItem("None")
			for defect in self.DefectsDiagram.defects_data.keys():
				if "_" not in defect:
					continue
				if self.DefectsDiagram.defects_data[defect]["Extrinsic"] == "Yes":
					### Herein lies the graveyard of various tricks I tried to implement subscripts in QComboBox, and I stand
					###		before you to tell you that: it cannot be done. It's hopeless, trying to write subscripts in QComboBox,
					###		for some reason. For this unknown reason, you (the programmer or future Michael Toriyama) should not
					###		waste any more time with the unforgiving fact that neither unicode nor latex will show up properly
					###		as a text on QComboBox. Do not waste any more of your time trying... Please and thank you.
					#self.extrinsic_defect_selection_box.addItem(u" "+defect.split("_")[0]+"<sub>"+defect.split("_")[-1]+"</sub>")
					#self.extrinsic_defect_selection_box.addItem(u""+defect.split("_")[0]+"<sub>"+defect.split("_")[-1]+"</sub>")
					#self.extrinsic_defect_selection_box.addItem(u""+defect.split("_")[0]+"<sub>"+defect.split("_")[-1]+"</sub>")
					#self.extrinsic_defect_selection_box.addItem(defect.split("_")[0]+"$_{"+defect.split("_")[-1]+"}$")
					#self.extrinsic_defect_selection_box.addItem(defect.split("_")[0]+'\u2083'+defect.split("_")[-1])
					### LOL I don't even need the ones above, I'll just add the names of dopants available in the database.
					dopant = defect.split("_")[0]
					if dopant not in [self.dopant_selection_box.itemText(i) for i in range(self.dopant_selection_box.count())]:
						self.dopant_selection_box.addItem(dopant)
			self.dopant_selection_box.setCurrentIndex(0)
			self.dopant_selection_box.activated.connect(self.Update_ExtrinsicDefect)
			self.dopant_properties_widget_layout.addWidget(self.dopant_selection_box)
			
			
			if self.show_carrier_concentration:
				# Add synthesis temperature button
				self.dopant_properties_widget_layout.addWidget(self.defects_synthesis_temperature)
			
			
			self.defectsdiagram_window_layout.addWidget(self.dopant_properties_widget)
		

		else:
			if self.show_carrier_concentration:
				# Add synthesis temperature button
				self.defectsdiagram_window_layout.addWidget(self.defects_synthesis_temperature)
		
		


		
		# (WIDGET) Save defects diagram as figure
		self.defects_diagram_savefigure_button = QPushButton("Save Defects Diagram Figure")
		self.defects_diagram_savefigure_button.clicked[bool].connect(lambda: self.DefectsDiagram.SaveFigure())
		self.defectsdiagram_window_layout.addWidget(self.defects_diagram_savefigure_button)
	
	
	
	
	
	
	###############################################################################################
	################################# Generate Defects Diagram ####################################
	###############################################################################################
	
	def Generate_DefectsDiagram_Plot_Function(self):
		
		# This function specifies what happens when the user clicks the "Generate Defects Diagram" button.
		
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
		if self.DefectsDiagram.dopant != "None":
			self.DefectsDiagram.Initialize_Extrinsic_DefectsDiagram_Plot()
	
	
	
	
	
	def Update_ExtrinsicDefect_DeltaMu(self):
		
		# Obtain deltamu of dopant
		self.dopant_deltamu = float(self.dopant_chemical_potential_deltamu.text())
		
		# Recalculate defect formation energies
		self.DefectsDiagram.dopant_deltamu = self.dopant_deltamu
		self.DefectsDiagram.Calculate_DefectFormations()
		
		# Redraw defects diagram
		if self.DefectsDiagram.intrinsic_defect_plots != {}:
			self.DefectsDiagram.Update_Intrinsic_DefectsDiagram_Plot()
		if self.DefectsDiagram.extrinsic_defect_plots != {}:
			self.DefectsDiagram.Update_Extrinsic_DefectsDiagram_Plot()
		
		
		
		if self.show_carrier_concentration:
			
			# Recalculate carrier concentrations
			self.CarrierConcentration.dopant_deltamu = self.dopant_deltamu
			
			# Redraw carrier concentration
			self.CarrierConcentration.Update_CarrierConcentration_Plot()

			# Plot the equilibrium Fermi energy
			if (self.DefectsDiagram.intrinsic_defect_plots != {}) and (self.DefectsDiagram.extrinsic_defect_plots != {}):
				self.Update_Equilibrium_Fermi_Energy_Temperature()
	
	
	
	def Update_ExtrinsicDefect(self):
		
		# Obtain selected dopant
		if self.dopant_selection_box.currentText() == "None":
			self.DefectsDiagram.dopant = "None"
			self.dopant_chemical_potential_deltamu.setEnabled(False)
		else:
			dopant_qcombobox_text = self.dopant_selection_box.currentText()
			self.DefectsDiagram.dopant = dopant_qcombobox_text
			self.dopant_chemical_potential_deltamu.setEnabled(True)
		
		# Set intrinsic chemical potential mu0 of dopant
		if self.DefectsDiagram.dopant == "None":
			self.dopant_chemical_potential_label.setText(u"\u0394"+"\u03BC"+"<sub>x</sub>")
			self.DefectsDiagram.dopant_mu0 = 0.0
		else:
			self.dopant_chemical_potential_label.setText(u"\u0394"+"\u03BC"+"<sub>"+self.DefectsDiagram.dopant+"</sub>")
			self.DefectsDiagram.dopant_mu0 = self.compounds_info[self.DefectsDiagram.dopant]["mu0"]
		
		# Reset deltamu of dopant
		self.dopant_chemical_potential_deltamu.setText("-0.0000")
		self.DefectsDiagram.dopant_deltamu = 0.0
		
		
		# Recalculate defect formation energies
		self.DefectsDiagram.Calculate_DefectFormations()
		
		# Draw fresh defects diagram
		self.DefectsDiagram.extrinsic_defect_plots = {}
		self.DefectsDiagram.defects_diagram_plot_drawing.remove()
		self.DefectsDiagram.defects_diagram_plot_drawing = self.DefectsDiagram.defects_diagram_plot_figure.add_subplot(111)
		self.DefectsDiagram.Activate_DefectsDiagram_Plot_Axes()
		self.DefectsDiagram.Initialize_Intrinsic_DefectsDiagram_Plot()
		if self.DefectsDiagram.dopant != "None":
			self.DefectsDiagram.Initialize_Extrinsic_DefectsDiagram_Plot()
		
		
		if self.show_carrier_concentration:
			
			# Recalculate carrier concentrations
			self.CarrierConcentration.dopant = self.DefectsDiagram.dopant
			self.CarrierConcentration.dopant_mu0 = self.DefectsDiagram.dopant_mu0
			self.CarrierConcentration.dopant_deltamu = self.DefectsDiagram.dopant_deltamu
			
			
			# Redraw carrier concentration
			self.CarrierConcentration.Initialize_CarrierConcentration_Plot()

			# Plot the equilibrium Fermi energy
			if self.DefectsDiagram.intrinsic_defect_plots != {}:
				self.Update_Equilibrium_Fermi_Energy_Temperature()

	
	
	def Update_SynthesisTemperature(self):
		
		# Obtain synthesis temperature (K)
		synthesis_temperature = self.defects_synthesis_temperature_box.text()
		
		# Check whether the written synthesis temperature is a possible temperature
		try:
			float(synthesis_temperature) > 0.0
			self.CarrierConcentration.synthesis_temperature = float(synthesis_temperature)
		except:
			self.defects_synthesis_temperature_box.setText("")
			self.CarrierConcentration.synthesis_temperature = None
		
		# Redraw carrier concentration plot with synthesis temperature
		self.CarrierConcentration.Initialize_CarrierConcentration_Plot()

		# Update equilibrium Fermi energy
		# NOTE: The Update_Equilibrium_Fermi_Energy_Temperature function is in window_carrier_concentration
		#		but it can be called because this defects diagram window object is inherited by the tab
		#		objects, which also inherit the carrier concentration window object.
		if self.DefectsDiagram.intrinsic_defect_plots != {}:
			self.Update_Equilibrium_Fermi_Energy_Temperature()




