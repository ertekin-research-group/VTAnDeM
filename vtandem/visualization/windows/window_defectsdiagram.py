
__author__ = 'Michael_Lidia_Jiaxing_Benita_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'


import numpy as np

import PyQt5
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *



class Window_DefectsDiagram(QWidget):
	
	def __init__(self, main_compound: str, DefectsDiagram_Plot, show_carrier_concentration, show_dopant):
		
		super().__init__()
		
		self.show_carrier_concentration = show_carrier_concentration
		self.show_dopant = show_dopant
		
		# Initialize main compound
		self.main_compound = main_compound
		
		# Initialize defects diagram object
		self.DefectsDiagram_Plot = DefectsDiagram_Plot
		self.DefectsDiagram_Plot.Activate_DefectsDiagram_Plot_Axes()
		
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
		self.defects_diagram_plot = self.DefectsDiagram_Plot.defects_diagram_plot_canvas
		self.defectsdiagram_window_layout.addWidget(self.defects_diagram_plot)
		
		# Y-axis limits for defects diagram
		self.defectsdiagram_viewport = QWidget()
		self.defectsdiagram_viewport_layout = QHBoxLayout(self.defectsdiagram_viewport)
		
		# Y-axis limits for defects diagram
		self.defectsdiagram_Ymin_label = QLabel(u"y"+"<sub>min</sub>")
		self.defectsdiagram_Ymin_label.setAlignment(Qt.AlignRight)
		self.defectsdiagram_viewport_layout.addWidget(self.defectsdiagram_Ymin_label)
		self.defectsdiagram_Ymin_box = QLineEdit("-2.0")
		self.defectsdiagram_Ymin_box.editingFinished.connect(lambda: self.DefectsDiagram_Plot.Update_WindowSize("YMin", self.defectsdiagram_Ymin_box))
		self.defectsdiagram_viewport_layout.addWidget(self.defectsdiagram_Ymin_box)
		self.defectsdiagram_Ymax_label = QLabel(u"y"+"<sub>max</sub>")
		self.defectsdiagram_Ymax_label.setAlignment(Qt.AlignRight)
		self.defectsdiagram_viewport_layout.addWidget(self.defectsdiagram_Ymax_label)
		self.defectsdiagram_Ymax_box = QLineEdit("2.0")
		self.defectsdiagram_Ymax_box.editingFinished.connect(lambda: self.DefectsDiagram_Plot.Update_WindowSize("YMax", self.defectsdiagram_Ymax_box))
		self.defectsdiagram_viewport_layout.addWidget(self.defectsdiagram_Ymax_box)
		self.defectsdiagram_window_layout.addWidget(self.defectsdiagram_viewport)
		
		
		
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
			for defect in self.DefectsDiagram_Plot.defects_data.keys():
				if "_" not in defect:
					continue
				if self.DefectsDiagram_Plot.defects_data[defect]["Extrinsic"] == "Yes":
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
				
				# Synthesis temperature
				self.defects_synthesis_temperature = QWidget()
				self.defects_synthesis_temperature_layout = QHBoxLayout(self.defects_synthesis_temperature)
				self.defects_synthesis_temperature_label = QLabel(u"T<sub>syn</sub> (K) = ")
				self.defects_synthesis_temperature_label.setAlignment(Qt.AlignRight)
				self.defects_synthesis_temperature_layout.addWidget(self.defects_synthesis_temperature_label)
				self.defects_synthesis_temperature_box = QLineEdit("")
				self.defects_synthesis_temperature_box.editingFinished.connect(self.Update_SynthesisTemperature)
				self.defects_synthesis_temperature_box.setEnabled(False)
				self.defects_synthesis_temperature_layout.addWidget(self.defects_synthesis_temperature_box)
				self.dopant_properties_widget_layout.addWidget(self.defects_synthesis_temperature)
			
			
			self.defectsdiagram_window_layout.addWidget(self.dopant_properties_widget)
		
		
		
		
		
		# (WIDGET) Save defects diagram as figure
		self.defects_diagram_savefigure_button = QPushButton("Save Defects Diagram Figure")
		#self.defects_diagram_savefigure_button.clicked[bool].connect(lambda: self.SaveFigure("Defects Diagram"))
		self.defects_diagram_savefigure_button.clicked[bool].connect(lambda: self.DefectsDiagram_Plot.SaveFigure("Defects Diagram"))
		self.defectsdiagram_window_layout.addWidget(self.defects_diagram_savefigure_button)
	
	
	
	
	
	
	def Update_ExtrinsicDefect_DeltaMu(self):
		
		# Obtain deltamu of dopant
		self.extrinsic_defect_deltamu = float(self.dopant_chemical_potential_deltamu.text())
		
		# Recalculate defect formation energies
		self.DefectsDiagram_Plot.extrinsic_defect_deltamu = self.extrinsic_defect_deltamu
		self.DefectsDiagram_Plot.Calculate_DefectFormations()
		
		# Redraw defects diagram
		if self.DefectsDiagram_Plot.intrinsic_defect_plots != {}:
			self.DefectsDiagram_Plot.Update_Intrinsic_DefectsDiagram_Plot_Plot()
		if self.DefectsDiagram_Plot.extrinsic_defect_plots != {}:
			self.DefectsDiagram_Plot.Update_Extrinsic_DefectsDiagram_Plot_Plot()
		
		
		
		if self.show_carrier_concentration:
			
			# Recalculate carrier concentrations
			self.CarrierConcentration.extrinsic_defect_deltamu = self.extrinsic_defect_deltamu
			
			# Redraw carrier concentration
			if self.CarrierConcentration.carrier_concentration_intrinsic_defect_hole_plot != None:
				self.CarrierConcentration.Update_HoleConcentration_Plot()
			if self.CarrierConcentration.carrier_concentration_intrinsic_defect_electron_plot != None:
				self.CarrierConcentration.Update_ElectronConcentration_Plot()
			
			# Plot the equilibrium Fermi energy
			if (self.DefectsDiagram_Plot.intrinsic_defect_plots != {}) and (self.DefectsDiagram_Plot.extrinsic_defect_plots != {}):
				self.Update_Equilibrium_Fermi_Energy_Temperature()
	
	
	
	def Update_ExtrinsicDefect(self):
		
		# Obtain selected dopant
		if self.dopant_selection_box.currentText() == "None":
			self.DefectsDiagram_Plot.dopant = "None"
			self.dopant_chemical_potential_deltamu.setEnabled(False)
		else:
			dopant_qcombobox_text = self.dopant_selection_box.currentText()
			#self.DefectsDiagram_Plot.dopant = dopant_qcombobox_text.split("->")[0]+"_"+dopant_qcombobox_text.split("->")[1]
			self.DefectsDiagram_Plot.dopant = dopant_qcombobox_text
			self.dopant_chemical_potential_deltamu.setEnabled(True)
		
		# Set intrinsic chemical potential mu0 of dopant
		if self.DefectsDiagram_Plot.dopant == "None":
			self.dopant_chemical_potential_label.setText(u"\u0394"+"\u03BC"+"<sub>x</sub>")
			self.DefectsDiagram_Plot.dopant_mu0 = 0.0
		else:
			self.dopant_chemical_potential_label.setText(u"\u0394"+"\u03BC"+"<sub>"+self.DefectsDiagram_Plot.dopant+"</sub>")
			self.DefectsDiagram_Plot.dopant_mu0 = self.compounds_info[self.DefectsDiagram_Plot.dopant]["mu0"]
		
		# Reset deltamu of dopant
		self.dopant_chemical_potential_deltamu.setText("-0.0000")
		self.DefectsDiagram_Plot.dopant_deltamu = 0.0
		
		# Recalculate defect formation energies
		self.DefectsDiagram_Plot.dopant = self.extrinsic_defect
		self.DefectsDiagram_Plot.extrinsic_defect_mu0 = self.extrinsic_defect_mu0
		self.DefectsDiagram_Plot.extrinsic_defect_deltamu = self.extrinsic_defect_deltamu
		self.DefectsDiagram_Plot.Calculate_DefectFormations()
		
		# Draw fresh defects diagram
		self.DefectsDiagram_Plot.extrinsic_defect_plots = {}
		self.DefectsDiagram_Plot.defects_diagram_plot_drawing.remove()
		self.DefectsDiagram_Plot.defects_diagram_plot_drawing = self.DefectsDiagram_Plot.defects_diagram_plot_figure.add_subplot(111)
		self.DefectsDiagram_Plot.Activate_DefectsDiagram_Plot_Plot_Axes()
		self.DefectsDiagram_Plot.Initialize_Intrinsic_DefectsDiagram_Plot_Plot()
		if self.DefectsDiagram_Plot.dopant != "None":
			self.DefectsDiagram_Plot.Initialize_Extrinsic_DefectsDiagram_Plot_Plot()
		
		
		
		if self.show_carrier_concentration:
			
			# Recalculate carrier concentrations
			self.CarrierConcentration.extrinsic_defect = self.extrinsic_defect
			self.CarrierConcentration.extrinsic_defect_mu0 = self.extrinsic_defect_mu0
			self.CarrierConcentration.extrinsic_defect_deltamu = self.extrinsic_defect_deltamu
			
			# Redraw carrier concentration
			if self.CarrierConcentration.carrier_concentration_intrinsic_defect_hole_plot != None:
				self.CarrierConcentration.Initialize_HoleConcentration_Plot()
			if self.CarrierConcentration.carrier_concentration_intrinsic_defect_electron_plot != None:
				self.CarrierConcentration.Initialize_ElectronConcentration_Plot()
			
			# Plot the equilibrium Fermi energy
			if self.DefectsDiagram_Plot.intrinsic_defect_plots != {}:
				self.Update_Equilibrium_Fermi_Energy_Temperature()
	
	
	
	
	def Update_SynthesisTemperature(self):
		
		# Obtain synthesis temperature (K)
		synthesis_temperature = self.defects_synthesis_temperature_box.text()
		
		# Check whether the written synthesis temperature is a possible temperature
		try:
			float(synthesis_temperature)
		except:
			self.defects_synthesis_temperature_box.setText("")
			self.CarrierConcentration.synthesis_temperature = None
			if self.CarrierConcentration.carrier_concentration_intrinsic_defect_hole_plot != None:
				self.CarrierConcentration.Initialize_HoleConcentration_Plot()
			if self.CarrierConcentration.carrier_concentration_intrinsic_defect_electron_plot != None:
				self.CarrierConcentration.Initialize_ElectronConcentration_Plot()
			return
		if float(synthesis_temperature) <= 0.0:
			self.defects_synthesis_temperature_box.setText("")
			self.CarrierConcentration.synthesis_temperature = None
			if self.CarrierConcentration.carrier_concentration_intrinsic_defect_hole_plot != None:
				self.CarrierConcentration.Initialize_HoleConcentration_Plot()
			if self.CarrierConcentration.carrier_concentration_intrinsic_defect_electron_plot != None:
				self.CarrierConcentration.Initialize_ElectronConcentration_Plot()
			return
		
		# Update synthesis temperature in CarrierConcentration object
		self.CarrierConcentration.synthesis_temperature = float(synthesis_temperature)
		
		# Redraw carrier concentration plot with synthesis temperature
		if self.CarrierConcentration.carrier_concentration_intrinsic_defect_hole_plot != None:
			self.CarrierConcentration.Initialize_HoleConcentration_Plot()
		if self.CarrierConcentration.carrier_concentration_intrinsic_defect_electron_plot != None:
			self.CarrierConcentration.Initialize_ElectronConcentration_Plot()
	
	
	
	
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
				if figure_type == "Defects Diagram":
					self.DefectsDiagram_Plot.defects_diagram_plot_figure.savefig(filename, bbox_inches='tight')
			else:
				if figure_type == "Defects Diagram":
					self.DefectsDiagram_Plot.defects_diagram_plot_figure.savefig(filename+"."+extension, bbox_inches='tight')
	"""


