
from __future__ import unicode_literals
__author__ = 'Michael_Lidia_Jiaxing_Elif'
__name__ = 'VTanDeM'


###############################################################################################################################
###############################################################################################################################
################################################### Import Libraries ##########################################################
###############################################################################################################################
###############################################################################################################################

import os
import sys
import numpy as np
import pandas as pd
import itertools
import periodictable
import copy
from math import exp, log, floor
from scipy import integrate

import PyQt5
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
from labellines import labelLine, labelLines


###############################################################################################################################
###############################################################################################################################
############################################### Main Application Window #######################################################
###############################################################################################################################
###############################################################################################################################

class Quaternary_ApplicationWindow(QMainWindow):
	
	
	###############################################################################################
	########################### Declare and Store Object Variables ################################
	###############################################################################################
	resized = pyqtSignal()
	def __init__(self, parent = None, main_compound = None, first_species = None, second_species = None, third_species = None, fourth_species = None):	
		
		QMainWindow.__init__(self)
		
		self.font = {'family': 'Arial',
				'color':  'black',
				'weight': 'normal',
				'size': 10 }
		
		self.Setup_Window_Framework()
		
		self.first_species = first_species
		self.second_species = second_species
		self.third_species = third_species
		self.fourth_species = fourth_species
		self.species_list = [self.first_species, self.second_species, self.third_species, self.fourth_species]					
		self.species_list_radiobuttons = [self.first_species, self.second_species, self.third_species, self.fourth_species]		
		
		self.main_compound = main_compound
		
		self.compounds_info = {}		
		self.Obtain_Compounds_Data()
		
		self.df0=pd.read_csv("Defects_Tracker.csv")
		
		self.main_compound_number_first_specie = self.compounds_info[self.main_compound][self.first_species]	
		self.main_compound_number_second_specie = self.compounds_info[self.main_compound][self.second_species]	
		self.main_compound_number_third_specie = self.compounds_info[self.main_compound][self.third_species]	
		self.main_compound_number_fourth_specie = self.compounds_info[self.main_compound][self.fourth_species]	
		self.main_compound_enthalpy = self.compounds_info[self.main_compound]["enthalpy"]	
		

		self.phasediagram_endpoints = min(self.main_compound_enthalpy/self.main_compound_number_first_specie, self.main_compound_enthalpy/self.main_compound_number_second_specie, self.main_compound_enthalpy/self.main_compound_number_third_specie, self.main_compound_enthalpy/self.main_compound_number_fourth_specie)
		

		self.mu4_value_array = np.linspace(-4.0, -0.0, 1001)	
		self.mu4 = -0.0											
		
		self.widgets = QWidget()						

		self.setCentralWidget(self.widgets)				
		self.widgets_grid = QGridLayout(self.widgets)	

		# (WIDGET) Title of compound
		quaternary_compound_title_formal = self.Compound_Name_Formal(self.main_compound)	
		self.quaternary_compound_title = QLabel(quaternary_compound_title_formal)			
		self.quaternary_compound_title.setAlignment(Qt.AlignCenter)							
		self.quaternary_compound_title_font = QFont("serif", 24, QFont.Bold) 				
		self.quaternary_compound_title.setFont(self.quaternary_compound_title_font)			
		self.widgets_grid.addWidget(self.quaternary_compound_title, 0, 0, 1, 4)				
		
		
		# (WIDGET) Buttons to assign species' numbers (defined by function in this script)
		self.Activate_Species_Buttons(button_name = "First\nSpecies", button_x_position = 1, button_y_position = 0)
		self.Activate_Species_Buttons(button_name = "Second\nSpecies", button_x_position = 1, button_y_position = 1)
		self.Activate_Species_Buttons(button_name = "Third\nSpecies", button_x_position = 1, button_y_position = 2)
		self.Activate_Species_Buttons(button_name = "Fourth\nSpecies", button_x_position = 1, button_y_position = 3)
		
		
		# (WIDGET) Display ports for mu values
		self.muvalue_display_widget = QWidget()
		self.muvalue_display_layout = QHBoxLayout(self.muvalue_display_widget)
		self.muvalue_display_layout.setContentsMargins(0, 0, 0, 0)
		self.muvalue_display_layout.setSpacing(0)
		self.widgets_grid.addWidget(self.muvalue_display_widget, 2, 0, 1, 3)
		self.Activate_MuValue_Displays()
		
		
		# (WIDGET) Generate quaternary phase diagram plot button
		self.generate_quaternary_plot_button_widget = QWidget()
		self.generate_quaternary_plot_button_layout = QVBoxLayout(self.generate_quaternary_plot_button_widget)
		self.widgets_grid.addWidget(self.generate_quaternary_plot_button_widget, 3, 0, 1, 3)
		self.Activate_Generate_PhaseDiagram_Plot_Button()
		
		
		# (WIDGET) Generate defects diagram plot button
		self.generate_defects_plot_button_widget = QWidget()
		self.generate_defects_plot_button_layout = QVBoxLayout(self.generate_defects_plot_button_widget)
		self.widgets_grid.addWidget(self.generate_defects_plot_button_widget, 4, 0, 1, 3)
		self.defects_diagram_Yaxis_widget = QWidget()
		self.defects_diagram_Yaxis_layout = QHBoxLayout(self.defects_diagram_Yaxis_widget)
		self.widgets_grid.addWidget(self.defects_diagram_Yaxis_widget, 4, 8, 1, 1)	
		self.defects_diagram_Yaxis_box_widget = QWidget()
		self.defects_diagram_Yaxis_box_layout = QHBoxLayout(self.defects_diagram_Yaxis_box_widget)
		self.widgets_grid.addWidget(self.defects_diagram_Yaxis_box_widget, 4, 9, 1, 2)	
		self.defects_diagram_savefigure_widget = QWidget()
		self.defects_diagram_savefigure_layout = QHBoxLayout(self.defects_diagram_savefigure_widget)		
		self.widgets_grid.addWidget(self.defects_diagram_savefigure_widget, 4, 11, 1, 1)	

		self.ymin = -0.05
		self.ymax = 0.6
		self.EVBM = 0.0
		self.ECBM = 0.0
		self.defect_plots = {}
		self.Activate_Generate_DefectsDiagram_Plot_Tools()
		
		
		# (WIDGET) Quaternary phase diagram plot
		self.quaternary_plot_figure = plt.figure()
		self.quaternary_plot_figure.subplots_adjust(left=0.05, right=0.8)
		self.quaternary_plot_drawing = self.quaternary_plot_figure.add_subplot(111)
		self.quaternary_plot_canvas = FigureCanvas(self.quaternary_plot_figure)
		self.widgets_grid.addWidget(self.quaternary_plot_canvas, 0, 4, 4, 4)
		self.Activate_PhaseDiagram_Plot_Axes()
		
		
		# (WIDGET) Defects diagram plot
		self.defects_diagram_plot_figure = plt.figure()
		self.defects_diagram_plot_drawing = self.defects_diagram_plot_figure.add_subplot(111)
		self.defects_diagram_plot_canvas = FigureCanvas(self.defects_diagram_plot_figure)
		self.widgets_grid.addWidget(self.defects_diagram_plot_canvas, 0, 8, 4, 4)
		self.Activate_DefectsDiagram_Plot_Axes()
		
		
		# (WIDGET) Fourth species slider
		self.fourth_species_slider_widget = QWidget()
		self.fourth_species_slider_layout = QHBoxLayout(self.fourth_species_slider_widget)
		self.widgets_grid.addWidget(self.fourth_species_slider_widget, 4, 4, 1, 4)
		self.Activate_Fourth_Species_Slider()
		
		
		# Initialize necessary objects in app
		self.main_compound_plot = None
		self.competing_compound_plots = {}
		self.competing_compound_colorwheel = {}
		self.PSR_vertices_plot = None
		self.PSR_vertices = []
		self.defect_labels = {}
		self.listdH_def = []
		self.list_carriers = []
		
		# Clicked point in phase diagram
		self.pressed_point = self.quaternary_plot_figure.canvas.mpl_connect('button_press_event', self.Pressed_Point)
		self.pressed_point_plot, = self.quaternary_plot_drawing.plot([], [], color="red", marker="o")
		self.pressed_mu1 = 0.0
		self.pressed_mu2 = 0.0
		self.pressed_mu3 = 0.0
		self.dmu_a = 0
		self.dmu_b = 0
		self.dmu_c = 0
		self.dmu_d = 0

		self.k=8.617e-5
		self.vol=229.36e-24
		self.listdH_def=[]
		self.def_carriers = []

		# (WIDGET) Carrier Concentration Plot
		self.carrier_concentration_plot_figure = plt.figure()
		self.carrier_concentration_plot_drawing = self.carrier_concentration_plot_figure.add_subplot(111)
		self.carrier_concentration_plot_canvas = FigureCanvas(self.carrier_concentration_plot_figure)
		self.widgets_grid.addWidget(self.carrier_concentration_plot_canvas, 0, 12, 4, 4)
		self.Activate_CarrierConcentration_Plot_Axes()

		# (WIDGET) Generate Carrier Concentration Button
		self.generate_carrier_concentration_plot_button_widget = QWidget()
		self.generate_carrier_concentration_plot_button_layout = QVBoxLayout(self.generate_carrier_concentration_plot_button_widget)
		self.widgets_grid.addWidget(self.generate_carrier_concentration_plot_button_widget, 5, 0, 1, 3)

		# (WIDGET) Equilibrium fermi Energy 
		self.Equilibrium_Fermi_Energy_widget = QWidget()
		self.Equilibrium_Fermi_Energy_layout = QHBoxLayout(self.Equilibrium_Fermi_Energy_widget)
		self.widgets_grid.addWidget(self.Equilibrium_Fermi_Energy_widget, 4, 12, 1, 1)	

		self.Equilibrium_Fermi_Energy_display_widget = QWidget()
		self.Equilibrium_Fermi_Energy_display_layout = QHBoxLayout(self.Equilibrium_Fermi_Energy_widget)
		self.widgets_grid.addWidget(self.Equilibrium_Fermi_Energy_display_widget, 4, 14, 1, 1)
		self.Activate_Generate_CarrierConcentration_Plot_Tools()


		self.show()	
	
	
	###############################################################################################
	################################# Main Window Framework #######################################
	###############################################################################################
	
	def Setup_Window_Framework(self):
		

		newappAction = QAction('New', self)
		newappAction.setStatusTip('Launch New Application (Tips by Michael Toriyama TM)')
		

		importMenu = QMenu('Import', self)
		importAction_csv = QAction('&Import CSV', self)
		importAction_other = QAction('&Import other', self)
		importAction_csv.setStatusTip('Import CSV (Tips by Michael Toriyama TM)')
		importAction_other.setStatusTip('Import Other Stuff (Tips by Michael Toriyama TM)')
		importMenu.addAction(importAction_csv)
		importMenu.addAction(importAction_other)
		
		exitAction = QAction('Exit', self)
		exitAction.setStatusTip('Exit Application (Tips by Michael Toriyama TM)')
		exitAction.triggered.connect(qApp.quit)
		
		menubar = self.menuBar()	
		fileMenu = menubar.addMenu('&File')
		
		fileMenu.addAction(newappAction)
		fileMenu.addMenu(importMenu)
		fileMenu.addAction(exitAction)
		
		mainwindow_title = "VTanDeM"
		self.setWindowTitle(mainwindow_title)
		QWidget.setFixedSize(self,1400,600)	
	
	
	###############################################################################################
	############################# Obtain DFT Data of Compounds ####################################
	###############################################################################################
	
	def Obtain_Compounds_Data(self):
		
		compounds_data = pd.read_csv("Compounds_Tracker.csv", header=0, index_col=0)
		header_information = list(compounds_data)
		
		for compound, compound_info in compounds_data.iterrows():						
			if compound in self.compounds_info.keys():									
				print("WARNING: '"+compound+"' is already in the list of compounds! This compound will thus not be recorded.")
				continue
			else:
				self.compounds_info[compound] = {}
				for header_specie, compound_data_info in zip(header_information, list(compound_info)):
					self.compounds_info[compound][header_specie] = compound_data_info	
				for specie in self.species_list:
					if specie not in self.compounds_info[compound].keys():
						self.compounds_info[compound][specie] = 0
	
	
	
	###############################################################################################
	########################## Rewrite Compound Name Latex-Style ##################################
	###############################################################################################
	
	def Compound_Name_Formal(self, compound_name):
		

		compound_species_info = self.compounds_info[compound_name]
		
		compound_name_formal = ""
		for species in self.species_list:				
			if compound_species_info[species] == 0:
				continue								
			elif compound_species_info[species] == 1:
				compound_name_formal += species			
			elif compound_species_info[species] > 1:
				compound_name_formal += species+"<sub>"+str(int(compound_species_info[species]))+"</sub>"	
		
		return compound_name_formal
	
	
	
	###############################################################################################
	#################################### Species Buttons ##########################################
	###############################################################################################
	
	def Activate_Species_Buttons(self, button_name, button_x_position, button_y_position):

		
		species_radiobutton_widget = QWidget()									
																				
		species_radiobutton_layout = QVBoxLayout(species_radiobutton_widget)	
																				
																				
		species_radiobutton_layout.setSpacing(0)

		
		species_radiobutton_label = QLabel(button_name)					
		species_radiobutton_label.setMargin(0)

		species_radiobutton_label.setAlignment(Qt.AlignCenter)			
		species_radiobutton_layout.addWidget(species_radiobutton_label)	
		
		first_species_button = QRadioButton(self.species_list_radiobuttons[0])										
		first_species_button.toggled.connect(lambda:self.Update_Species_Button(first_species_button, button_name))	
		species_radiobutton_layout.addWidget(first_species_button)													
		
		second_species_button = QRadioButton(self.species_list_radiobuttons[1])										
		second_species_button.toggled.connect(lambda:self.Update_Species_Button(second_species_button, button_name))
		species_radiobutton_layout.addWidget(second_species_button)
		
		third_species_button = QRadioButton(self.species_list_radiobuttons[2])										
		third_species_button.toggled.connect(lambda:self.Update_Species_Button(third_species_button, button_name))
		species_radiobutton_layout.addWidget(third_species_button)
		
		fourth_species_button = QRadioButton(self.species_list_radiobuttons[3])										
		fourth_species_button.toggled.connect(lambda:self.Update_Species_Button(fourth_species_button, button_name))
		species_radiobutton_layout.addWidget(fourth_species_button)
		

		if "First" in button_name:
			first_species_button.setChecked(True)
		elif "Second" in button_name:
			second_species_button.setChecked(True)
		elif "Third" in button_name:
			third_species_button.setChecked(True)
		elif "Fourth" in button_name:
			fourth_species_button.setChecked(True)
		
	
		species_radiobutton_widget.setStyleSheet("background-color: rgb(255,255,255)")
		self.widgets_grid.addWidget(species_radiobutton_widget, button_x_position, button_y_position)
	
	
	def Update_Species_Button(self, species_button, button_name):
		
		# Basically changes the species that is called "first species", "second species", etc.
		if species_button.isChecked() == True:
			if button_name == "First Species":
				self.first_species = species_button.text()
			elif button_name == "Second Species":
				self.second_species = species_button.text()
			elif button_name == "Third Species":
				self.third_species = species_button.text()
			elif button_name == "Fourth Species":
				self.fourth_species = species_button.text()
		
		# Reorder the list of the species
		self.species = [self.first_species, self.second_species, self.third_species, self.fourth_species]
		# *Note: 	Clicking a radio button will change the order of the species, but this alone will NOT show any visible changes
		#			in the plots. These changes will be visible in the plots when the user clicks the "Generate Plot" button.
	
	
	
	###############################################################################################
	################################## Mu Values Displays #########################################
	###############################################################################################
	
	def Activate_MuValue_Displays(self):
		
		# Create widgets for the display ports of the mu values. These mu values will change based on where on the quaternary phase
		#	diagram plot the user clicks.
		# *Note: The "sub-main" widget that will contain these mu display ports for each species is already created in __init__ as
		#	the variable self.muvalue_display_widget. It's layout is self.muvalue_display_layout, and this is where each display
		#	port will be placed.
		
		self.display_label_font = QFont("Arial", 18)						# Define the font
		
		self.mu1_display_widget = QWidget()									# Create a widget for the display port
		self.mu1_display_layout = QVBoxLayout(self.mu1_display_widget)		# Create a vertical layout for the widget (should only contain the label and the mu4 number, stacked on top of each other)
		self.mu1_display_label = QLabel(u"\u0394"+"\u03BC"+"<sub>a</sub>")	# Create a label widget to display the "\Delta\mu_1" text above the actual mu1 value (unicode format)
		self.mu1_display_label.setFont(self.display_label_font)			# Set the font
		self.mu1_display_label.setAlignment(Qt.AlignCenter)					# Align the label to center
		self.mu1_display_layout.addWidget(self.mu1_display_label)			# Add the label
		self.mu1_display = QLCDNumber()										# Create the actual port that will display the current mu1 value
		self.mu1_display.setStyleSheet("QLCDNumber {color: black}")			# Color
		self.mu1_display_layout.addWidget(self.mu1_display)					# Add the display port
		
		self.mu2_display_widget = QWidget()
		self.mu2_display_layout = QVBoxLayout(self.mu2_display_widget)
		self.mu2_display_label = QLabel(u"\u0394"+"\u03BC"+"<sub>b</sub>")
		self.mu2_display_label.setFont(self.display_label_font)
		self.mu2_display_label.setAlignment(Qt.AlignCenter)
		self.mu2_display_layout.addWidget(self.mu2_display_label)
		self.mu2_display = QLCDNumber()
		self.mu2_display.setStyleSheet("QLCDNumber {color: black}")
		self.mu2_display_layout.addWidget(self.mu2_display)
		
		self.mu3_display_widget = QWidget()
		self.mu3_display_layout = QVBoxLayout(self.mu3_display_widget)
		self.mu3_display_label = QLabel(u"\u0394"+"\u03BC"+"<sub>c</sub>")
		self.mu3_display_label.setFont(self.display_label_font)
		self.mu3_display_label.setAlignment(Qt.AlignCenter)
		self.mu3_display_layout.addWidget(self.mu3_display_label)
		self.mu3_display = QLCDNumber()
		self.mu3_display.setStyleSheet("QLCDNumber {color: black}")
		self.mu3_display_layout.addWidget(self.mu3_display)
		
		self.muvalue_display_layout.addWidget(self.mu1_display_widget)
		self.muvalue_display_layout.addWidget(self.mu2_display_widget)
		self.muvalue_display_layout.addWidget(self.mu3_display_widget)
	
	
	
	###############################################################################################
	############################## Fourth Species Slide Bar #######################################
	###############################################################################################
	
	def Activate_Fourth_Species_Slider(self):
		
		# Create the slider widget for the chemical potential (mu) of the fourth species.
		# *Note: The widget that will contain the slider, as well as its (horizontal) layout are already defined in __init__.
		#	This function is simply to design and functionalize this "sub-main" widget.
		
		self.fourth_species_slider_label = QLabel(u"\u0394"+"\u03BC"+"<sub>d</sub>")		
		self.fourth_species_slider_label.setAlignment(Qt.AlignCenter)						
		self.fourth_species_slider_layout.addWidget(self.fourth_species_slider_label)		
		
		self.fourth_species_slider = QSlider(Qt.Horizontal)						
		self.fourth_species_slider.setMinimum(0)								
		self.fourth_species_slider.setMaximum(1000)								
		self.fourth_species_slider.setValue(1000)								
		self.fourth_species_slider.setSingleStep(1)								
		self.fourth_species_slider.setTickInterval(1)							 
																				# *Note: The slider widget cannot handle floating values, i.e. you can only define
																				#	integers for the maximum, minimum, etc. The way we get around that is by
																				#	saying that the mu4 value is chosen from an array (self.mu4_value_array), and 
																				#	the integer value of the slider is the INDEX of the array.
		self.fourth_species_slider_layout.addWidget(self.fourth_species_slider)	# Add the slider widget to the main slider widget
		
		self.fourth_species_slider.valueChanged.connect(self.Update_Fourth_Species_Slider)	# How to update the slider in the case that someone moves it
		
		self.fourth_species_slider_value_label = QLabel("{0:.4f}".format(-0.0))				# Create a label widget to display the current value of mu4
		self.fourth_species_slider_value_label.setAlignment(Qt.AlignCenter)					# Align the text to center
		self.fourth_species_slider_layout.addWidget(self.fourth_species_slider_value_label)	# Add this label widget to the main slider widget layout
	
	
	def Update_Fourth_Species_Slider(self):
		
		# This is how we update the properties of the slider given that the user touches it.
		fourth_species_mu4_value_index = self.fourth_species_slider.value()	# Obtain the (integer) value of the slider when the user uses the slider
		self.mu4 = self.mu4_value_array[fourth_species_mu4_value_index]		# Update the mu4 value by 1) using the value of the slider as the index and 
																			#	2) using the index to choose the mu4 value from self.mu4_value_array
		
		mu4_rounded = round(self.mu4, 4)												# Display the updated mu4 value as a float up to four decimal places
		self.fourth_species_slider_value_label.setText("{0:.4f}".format(mu4_rounded))	# Update the text display
		
		if (self.main_compound_plot != None) and (self.competing_compound_plots != {}):		# If the mu4 value changes (i.e. the user uses the slider), update 
																							#	the quaternary phase diagram.
			self.Plot_PhaseDiagram()
		
		if self.defect_plots != {}:
			self.Plot_DefectsDiagram('HSE06')
		
		self.pressed_mu3 = (self.main_compound_enthalpy - self.main_compound_number_first_specie*self.pressed_mu1 - self.main_compound_number_second_specie*self.pressed_mu2 - self.main_compound_number_fourth_specie*self.mu4) / self.main_compound_number_third_specie
		self.mu3_display.display(self.pressed_mu3)		# Every time the mu4 value changes (and the mu1 and mu2 values are held constant), the mu3 value changes
	
	
	
	###############################################################################################
	##################################### Phase Diagram ###########################################
	###############################################################################################
	
	def Activate_PhaseDiagram_Plot_Axes(self):
		
		# This function simply creates the axes of the quaternary phase diagram.
		
		self.quaternary_plot_drawing.set_xlim(self.phasediagram_endpoints, 0.0)
		self.quaternary_plot_drawing.set_ylim(self.phasediagram_endpoints, 0.0)
		self.quaternary_plot_drawing.set_xlabel("$\Delta\mu_{a}$ (eV)",fontdict=self.font,labelpad=12)
		self.quaternary_plot_drawing.set_ylabel("$\Delta\mu_{b}$ (eV)",fontdict=self.font,rotation=270,labelpad=20)
		self.quaternary_plot_drawing.xaxis.tick_top()
		self.quaternary_plot_drawing.yaxis.tick_right()
		self.quaternary_plot_drawing.tick_params(axis='both', labelsize=6)
		self.quaternary_plot_drawing.xaxis.set_label_position("top")
		self.quaternary_plot_drawing.yaxis.set_label_position("right")
		self.quaternary_plot_drawing.spines['left'].set_visible(False)
		self.quaternary_plot_drawing.spines['bottom'].set_visible(False)
		self.quaternary_plot_drawing.set_aspect("equal")
	
	
	def Update_PhaseDiagram_Plot_Axes(self):
		
		# This function changes the axes labels of the quaternary phase diagram when the user selects
		#	the species to be plotted. This will be used in conjunction with the 
		#	Activate_PhaseDiagram_Plot_Axes function, which creates the initial phase diagram plot.
		
		self.quaternary_plot_drawing.set_xlabel("$\Delta\mu_{"+self.first_species+"}$ (eV)",fontdict=self.font,labelpad=12)
		self.quaternary_plot_drawing.set_ylabel("$\Delta\mu_{"+self.second_species+"}$ (eV)",fontdict=self.font,rotation=270,labelpad=20)
	
	
	def Activate_Generate_PhaseDiagram_Plot_Button(self):
		
		# The user is given a button to generate the quaternary phase diagram plot. This is because
		#	when the user changes the species to be plotted, the program does not automatically
		#	plot the new phase diagram.
		
		self.generate_quaternary_plot_button_label = QLabel("Hello. Welcome to [insert name of app].")	
		self.generate_quaternary_plot_button_label.setAlignment(Qt.AlignCenter)
		self.generate_quaternary_plot_button_layout.addWidget(self.generate_quaternary_plot_button_label)
		
		self.generate_quaternary_plot_button = QPushButton("Generate Plot!")							
		self.generate_quaternary_plot_button.clicked[bool].connect(self.Generate_PhaseDiagram_Plot)
		self.generate_quaternary_plot_button_layout.addWidget(self.generate_quaternary_plot_button)
	

	def Generate_PhaseDiagram_Plot(self, event):

		# This function specifies what happens when the user clicks the "Generate Plot" button. The 
		#	only thing it needs to check is whether the elements chosen as the first, second, third, 
		#	and fourth species are unique.
		
		if len(self.species) > len(set(self.species)):	# Every time someone clicks the species' buttons, the self.species list gets updated.
														#	The "set" function checks the self.species list and omits any that are repeated.
														#	If any are repeated, then the chosen species are not unique.
			self.generate_quaternary_plot_button_label.setText("Yo... Pick UNIQUE elements!")	
			self.fourth_species_slider.setEnabled(False)	
		
		else:
			
			self.generate_quaternary_plot_button_label.setText("")
			
			# Update the number of each species in the main compound based on the above change
			self.main_compound_number_first_specie = self.compounds_info[self.main_compound][self.first_species]
			self.main_compound_number_second_specie = self.compounds_info[self.main_compound][self.second_species]
			self.main_compound_number_third_specie = self.compounds_info[self.main_compound][self.third_species]
			self.main_compound_number_fourth_specie = self.compounds_info[self.main_compound][self.fourth_species]
			
			# Reset the slide bar
			self.fourth_species_slider.setEnabled(True)
			self.fourth_species_slider_label.setText(u"\u0394\u03BC<sub>"+self.fourth_species+"</sub>")
			endpoint_slidebar = self.main_compound_enthalpy / self.main_compound_number_fourth_specie
			self.mu4_value_array = np.linspace(endpoint_slidebar, -0.0, 1001)
			self.fourth_species_slider.setValue(1000)
			
			# Reset the mu value displays
			self.mu1_display_label.setText(u"\u0394\u03BC<sub>"+self.first_species+"</sub>")
			self.mu2_display_label.setText(u"\u0394\u03BC<sub>"+self.second_species+"</sub>")
			self.mu3_display_label.setText(u"\u0394\u03BC<sub>"+self.third_species+"</sub>")
			self.mu1_display.display(0.0)
			self.mu2_display.display(0.0)
			self.mu3_display.display(0.0)
			
			# Set up the new plot
			self.quaternary_plot_drawing.remove()
			self.quaternary_plot_drawing = self.quaternary_plot_figure.add_subplot(111)
			
			# Reset the plots of the main and competing compounds
			self.main_compound_plot = None
			self.competing_compound_plots = {}
			
			# Reset the axes of the plot
			self.Activate_PhaseDiagram_Plot_Axes()
			self.Update_PhaseDiagram_Plot_Axes()
			
			# Plot the phase stability diagram of the quaternary compound using the new settings
			self.Plot_PhaseDiagram()
			
			# Reset clicked point
			self.pressed_point_plot, = self.quaternary_plot_drawing.plot([], [], color="red", marker="o")
	
	
	def Plot_PhaseDiagram(self):
		
		# Plotting the phase diagram from DFT data is one of the main features of this app. This function
		#	single-handedly plots the phase diagram, so it's arguably one of the most important block of
		#	code in this script.
		# This function will be used in two forms:
		#	1) when the user generates the phase diagram by clicking the "Generate Plot" button
		#	2) when the user changes the mu4 value
		# Because there are different forms in which this function will be used, the programmer should take
		#	care when editing this block of code.

		main_compound_enthalpy_adjusted = self.main_compound_enthalpy - self.main_compound_number_fourth_specie*self.mu4	# Enthalpy adjusted for mu4 value (main compound)
		
		main_compound_deltamu_first_species = np.linspace(self.main_compound_enthalpy/self.main_compound_number_first_specie, 0, 1000)	# Array of possible mu1 values in the plot
		
		main_compound_stability_limit = (main_compound_enthalpy_adjusted - self.main_compound_number_first_specie*main_compound_deltamu_first_species) / self.main_compound_number_second_specie	# Array that represents the (diagonal) stability limit of the main compound
		
		# Check if the phase diagram has already been plotted
		try:
			self.main_compound_plot.set_ydata(main_compound_stability_limit)	# If the plot's already there, simply update the stability limit of the main compound
		except:
			self.main_compound_plot, = self.quaternary_plot_drawing.plot(main_compound_deltamu_first_species, main_compound_stability_limit)	# If the main compound stabilitiy limit hasn't been plotted yet, plot it
		
		# Find the bounds of the quaternary phase stability region (for shading the stability region of the quaternary compound)
		stability_minimum_bound = []
		stability_maximum_bound = []
		
		vertical_left_values = []
		vertical_right_values = []
		
		# Loop through all compounds in the database
		for competing_compound_index, competing_compound in enumerate(self.compounds_info.keys()):

			# Check if the compound should be considered a "competing compound"
			foreign_compound = False
			for element in periodictable.elements:
				if str(element) in self.species:
					continue
				elif str(element) not in self.compounds_info[competing_compound].keys():
					continue
				else:
					if self.compounds_info[competing_compound][str(element)] != 0.0:
						foreign_compound = True
			
			if foreign_compound:	# Foreign compound (e.g. compound containing at least one element not in the quaternary compound)
				continue
			elif (competing_compound in self.species) or (competing_compound == self.main_compound):	# Single-element compounds and the main compound are not considered "competing compounds" of the main quaternary
				continue
			else:
				#print("Analyzing the competing compound '"+competing_compound+"'...")
				
				competing_compound_number_first_specie = self.compounds_info[competing_compound][self.first_species]
				competing_compound_number_second_specie = self.compounds_info[competing_compound][self.second_species]
				competing_compound_number_third_specie = self.compounds_info[competing_compound][self.third_species]
				competing_compound_number_fourth_specie = self.compounds_info[competing_compound][self.fourth_species]
				
				competing_compound_enthalpy = self.compounds_info[competing_compound]["enthalpy"]
				competing_compound_enthalpy_adjusted = competing_compound_enthalpy - competing_compound_number_fourth_specie*self.mu4		# Enthalpy adjusted for mu4 value (competing compound)
				
				difference_enthalpy_adjusted = competing_compound_enthalpy_adjusted - (competing_compound_number_third_specie / self.main_compound_number_third_specie) * main_compound_enthalpy_adjusted
				
				coefficient_first_specie = competing_compound_number_first_specie - (self.main_compound_number_first_specie * competing_compound_number_third_specie) / self.main_compound_number_third_specie
				coefficient_second_specie = competing_compound_number_second_specie - (self.main_compound_number_second_specie * competing_compound_number_third_specie) / self.main_compound_number_third_specie
				
				if (coefficient_first_specie == 0.0) and (coefficient_second_specie == 0.0):
					continue
					print("ZOMG SOMETHING'S WRONG DUUUUUDE")
					print("Compound may not be stoichiomatrically balanced???")
					
				elif (coefficient_first_specie != 0.0) and (coefficient_second_specie == 0.0):
					constant_deltamu_first_species = difference_enthalpy_adjusted / coefficient_first_specie
					if coefficient_first_specie > 0.0:
						vertical_left_values.append(constant_deltamu_first_species)
					elif coefficient_first_specie < 0.0:
						vertical_right_values.append(constant_deltamu_first_species)
					
					competing_compound_deltamu_first_species = np.ones(len(main_compound_deltamu_first_species)) * constant_deltamu_first_species
					competing_compound_deltamu_second_species = np.linspace(self.main_compound_enthalpy/self.main_compound_number_second_specie, 0, len(competing_compound_deltamu_first_species))
					
					main_compound_stability_limit_vertical = ( main_compound_enthalpy_adjusted - self.main_compound_number_first_specie*constant_deltamu_first_species ) / self.main_compound_number_second_specie
					competing_compound_deltamu_first_species_limit = [competing_compound_deltamu_first_species[i] for i in range(len(competing_compound_deltamu_first_species)) if (main_compound_stability_limit_vertical < competing_compound_deltamu_second_species[i])]
					competing_compound_deltamu_second_species_limit = [competing_compound_deltamu_second_species[i] for i in range(len(competing_compound_deltamu_first_species)) if (main_compound_stability_limit_vertical < competing_compound_deltamu_second_species[i])]
					
				elif (coefficient_first_specie == 0.0) and (coefficient_second_specie != 0.0):
					competing_compound_deltamu_first_species = copy.deepcopy(main_compound_deltamu_first_species)
					constant_deltamu_second_species = difference_enthalpy_adjusted / coefficient_second_specie
					competing_compound_deltamu_second_species = np.ones(len(competing_compound_deltamu_first_species)) * constant_deltamu_second_species
					
					if coefficient_second_specie > 0.0:
						stability_maximum_bound.append(competing_compound_deltamu_second_species)
					elif coefficient_second_specie < 0.0:
						stability_minimum_bound.append(competing_compound_deltamu_second_species)
					
					competing_compound_deltamu_first_species_limit = [competing_compound_deltamu_first_species[i] for i in range(len(competing_compound_deltamu_first_species)) if (main_compound_stability_limit[i] < competing_compound_deltamu_second_species[i])]
					competing_compound_deltamu_second_species_limit = [competing_compound_deltamu_second_species[i] for i in range(len(competing_compound_deltamu_first_species)) if (main_compound_stability_limit[i] < competing_compound_deltamu_second_species[i])]
				
				elif (coefficient_first_specie != 0.0) and (coefficient_second_specie != 0.0):
					competing_compound_deltamu_first_species = copy.deepcopy(main_compound_deltamu_first_species)
					competing_compound_deltamu_second_species = ( difference_enthalpy_adjusted - coefficient_first_specie*competing_compound_deltamu_first_species ) / coefficient_second_specie
					
					if coefficient_second_specie > 0.0:
						stability_maximum_bound.append(competing_compound_deltamu_second_species)
					elif coefficient_second_specie < 0.0:
						stability_minimum_bound.append(competing_compound_deltamu_second_species)
					
					competing_compound_deltamu_first_species_limit = [competing_compound_deltamu_first_species[i] for i in range(len(competing_compound_deltamu_first_species)) if (main_compound_stability_limit[i] < competing_compound_deltamu_second_species[i])]
					competing_compound_deltamu_second_species_limit = [competing_compound_deltamu_second_species[i] for i in range(len(competing_compound_deltamu_first_species)) if (main_compound_stability_limit[i] < competing_compound_deltamu_second_species[i])]
				
				try:
					# See if plot exists so it just needs to be updated
					self.competing_compound_plots[competing_compound].set_data(competing_compound_deltamu_first_species_limit, competing_compound_deltamu_second_species_limit)
				except:
					# If the plot doesn't exist initially, then create it
					self.competing_compound_plots[competing_compound], = self.quaternary_plot_drawing.plot(competing_compound_deltamu_first_species_limit, competing_compound_deltamu_second_species_limit, label=competing_compound)
					self.competing_compound_colorwheel[competing_compound] = self.competing_compound_plots[competing_compound].get_color()
		
		# Check if the legend for the phase diagram has been activated
		try:
			self.quaternary_plot_legend_widget
		except:
			self.Activate_Legend()
		
		stability_minimum_bound.append(main_compound_stability_limit)
		stability_maximum_bound.append(np.zeros(len(main_compound_deltamu_first_species)))
		stability_absolute_minimum = np.fromiter(map(max, zip(*itertools.chain(stability_minimum_bound))), dtype=np.float)
		stability_absolute_maximum = np.fromiter(map(min, zip(*itertools.chain(stability_maximum_bound))), dtype=np.float)
		
		self.main_compound_deltamu_first_species_cutoff = []
		self.stability_minimum_cutoff = []
		self.stability_maximum_cutoff = []
		
		for i in range(len(main_compound_deltamu_first_species)):
			if vertical_left_values != []:
				if (stability_absolute_minimum[i] < stability_absolute_maximum[i]) and (main_compound_deltamu_first_species[i] < min(vertical_left_values)):
					self.main_compound_deltamu_first_species_cutoff.append(main_compound_deltamu_first_species[i])
					self.stability_minimum_cutoff.append(stability_absolute_minimum[i])
					self.stability_maximum_cutoff.append(stability_absolute_maximum[i])
			if vertical_right_values != []:
				if (stability_absolute_minimum[i] < stability_absolute_maximum[i]) and (main_compound_deltamu_first_species[i] > max(vertical_right_values)):
					self.main_compound_deltamu_first_species_cutoff.append(main_compound_deltamu_first_species[i])
					self.stability_minimum_cutoff.append(stability_absolute_minimum[i])
					self.stability_maximum_cutoff.append(stability_absolute_maximum[i])
		
		try:
			self.phase_stability_region.remove()
			self.phase_stability_region = self.quaternary_plot_drawing.fill_between(self.main_compound_deltamu_first_species_cutoff, self.stability_maximum_cutoff, self.stability_minimum_cutoff, facecolor='0.75')
		except:
			self.phase_stability_region = self.quaternary_plot_drawing.fill_between(self.main_compound_deltamu_first_species_cutoff, self.stability_maximum_cutoff, self.stability_minimum_cutoff, facecolor='0.75')
		
		if self.phase_stability_region.get_paths() != []:
			PSR_Vertices = []
			PSR_Bound_Slope_Previous = None
			PSR_Bounding_Point_Previous = None
			tolerance = 1E-6
			for PSR_Bounds in self.phase_stability_region.get_paths()[0].iter_segments():
				PSR_Bounding_Point = PSR_Bounds[0]
				try:
					PSR_Bound_Slope = (PSR_Bounding_Point[1]-PSR_Bounding_Point_Previous[1]) / (PSR_Bounding_Point[0]-PSR_Bounding_Point_Previous[0])
				except:
					PSR_Bounding_Point_Previous = PSR_Bounding_Point
					PSR_Bound_Slope_Previous = 0.0
					continue
				if (PSR_Bound_Slope < PSR_Bound_Slope_Previous - tolerance) or (PSR_Bound_Slope > PSR_Bound_Slope_Previous + tolerance):
					PSR_Vertices.append(PSR_Bounding_Point_Previous)
					PSR_Vertices.append(PSR_Bounding_Point)
				PSR_Bound_Slope_Previous = PSR_Bound_Slope
				PSR_Bounding_Point_Previous = PSR_Bounding_Point
			PSR_Vertices_Omit = []
			for PSR_Vertices_Index in range(len(PSR_Vertices)-1):
				if (np.linalg.norm(PSR_Vertices[PSR_Vertices_Index]-PSR_Vertices[PSR_Vertices_Index+1]) < 0.01):
					PSR_Vertices_Omit.append(PSR_Vertices[PSR_Vertices_Index])
			PSR_Vertices_Unrepeated = [x for x in PSR_Vertices if (not any((x is y for y in PSR_Vertices_Omit)))]
			self.PSR_vertices = PSR_Vertices_Unrepeated
			if PSR_Vertices_Unrepeated == []:
				try:
					self.PSR_vertices_plot.remove()
				except:
					pass
			elif PSR_Vertices_Unrepeated != []:
				try:
					self.PSR_vertices_plot.remove()
					self.PSR_vertices_plot = self.quaternary_plot_drawing.scatter(*zip(*PSR_Vertices_Unrepeated), s=20, c='black')
				except:
					self.PSR_vertices_plot = self.quaternary_plot_drawing.scatter(*zip(*PSR_Vertices_Unrepeated), s=20, c='black')
					pass
		elif self.phase_stability_region.get_paths() == []:
			self.PSR_vertices = []
			try:
				self.PSR_vertices_plot.remove()
			except:
				pass
		
		# Draw the new phase diagram
		self.quaternary_plot_canvas.draw()
	
	
	
	def Pressed_Point(self, event):
		
		# This function reads in the coordinates of the point in the phase stability diagram plot that the user clicks.

		point_x = event.xdata	# x-coordinate
		point_y = event.ydata	# y-coordinate
		
		if (not isinstance(point_x, float)) or (not isinstance(point_y, float)):	# Check that the user presses somewhere on the plot (and not anywhere else)
			pass
		else:

			# Update the clicked mu values
			self.pressed_mu1 = point_x
			self.pressed_mu2 = point_y
			self.pressed_mu3 = (self.main_compound_enthalpy - self.main_compound_number_first_specie*self.pressed_mu1 - self.main_compound_number_second_specie*self.pressed_mu2 - self.main_compound_number_fourth_specie*self.mu4) / self.main_compound_number_third_specie
			
			# Update the display ports for the mu values
			self.mu1_display.display(self.pressed_mu1)
			self.mu2_display.display(self.pressed_mu2)
			self.mu3_display.display(self.pressed_mu3)
			
			# Update the red dot where the user clicked on the phase diagram
			self.pressed_point_plot.set_data([self.pressed_mu1], [self.pressed_mu2])
			self.quaternary_plot_figure.canvas.draw_idle()
			
			if self.defect_plots != {}:
				self.Plot_DefectsDiagram('HSE06')
	
	
	
	
	###############################################################################################
	######################################### Legend ##############################################
	###############################################################################################
	
	def Activate_Legend(self):
		
		self.quaternary_plot_legend_widget = QWidget()
		self.quaternary_plot_legend_layout = QGridLayout(self.quaternary_plot_legend_widget)
		
		quaternary_plot_legend_name = QLabel("Competing\nCompound\nColors")
		quaternary_plot_legend_name.setAlignment(Qt.AlignCenter)
		self.quaternary_plot_legend_layout.addWidget(quaternary_plot_legend_name, 0, 0, 1, 2)
		
		for competing_compound_index, competing_compound_name, competing_compound_color in zip(range(len(self.competing_compound_colorwheel)), self.competing_compound_colorwheel.keys(), self.competing_compound_colorwheel.values()):
			
			competing_compound_name_label = QLabel(self.Compound_Name_Formal(competing_compound_name))
			competing_compound_name_label.setAlignment(Qt.AlignCenter)
			self.quaternary_plot_legend_layout.addWidget(competing_compound_name_label, competing_compound_index+1, 0)
			
			competing_compound_color_label = QWidget()
			#competing_compound_color_label.resize(1, 100)
			competing_compound_color_label.setAutoFillBackground(True)
			competing_compound_color_palette = competing_compound_color_label.palette()
			color = QColor(competing_compound_color)
			competing_compound_color_palette.setColor(competing_compound_color_label.backgroundRole(), color)
			competing_compound_color_label.setPalette(competing_compound_color_palette)
			self.quaternary_plot_legend_layout.addWidget(competing_compound_color_label, competing_compound_index+1, 1)
		
		self.widgets_grid.addWidget(self.quaternary_plot_legend_widget, 2, 3, 2, 1)
	
	
	
	###############################################################################################
	################################### Defects Diagram ###########################################
	###############################################################################################
	
	def Activate_DefectsDiagram_Plot_Axes(self):
		
		self.defects_diagram_plot_drawing.set_xlim(-1.0, 0.0)
		self.defects_diagram_plot_drawing.set_ylim(-0.2, 2.0)
		self.defects_diagram_plot_drawing.set_xlabel("Fermi Energy (eV)", fontdict=self.font, labelpad=12)
		self.defects_diagram_plot_drawing.set_ylabel("$\Delta$H(q,D) (eV)", fontdict=self.font, rotation=90, labelpad=20)
		self.defects_diagram_plot_drawing.xaxis.tick_bottom()
		self.defects_diagram_plot_drawing.yaxis.tick_left()
		self.defects_diagram_plot_drawing.tick_params(axis='both', labelsize=6)
		self.defects_diagram_plot_drawing.xaxis.set_label_position("bottom")
		self.defects_diagram_plot_drawing.yaxis.set_label_position("left")
		#self.defects_diagram_plot_drawing.set_aspect("equal")
		plt.subplots_adjust(left = 0.2, right = 0.85, top = 0.9, bottom = 0.1)
		self.defects_diagram_plot_drawing.xaxis.labelpad = 2.5
		self.defects_diagram_plot_drawing.yaxis.labelpad = 2.0	
		self.defects_diagram_plot_drawing.set_aspect("auto")
	
	
	def Activate_Generate_DefectsDiagram_Plot_Tools(self):
		
		self.generate_defects_plot_button = QPushButton("Defect Diagram")
		self.generate_defects_plot_button.clicked[bool].connect(self.Generate_DefectsDiagram_Plot)
		self.generate_defects_plot_button_layout.addWidget(self.generate_defects_plot_button)
		
		self.defects_diagram_Yaxis_label = QLabel(u"y"+"<sub>min</sub>"+','+u"y"+"<sub>max</sub>")			
		self.defects_diagram_Yaxis_label.setAlignment(Qt.AlignCenter)
		self.defects_diagram_Yaxis_layout.addWidget(self.defects_diagram_Yaxis_label)
		self.defects_diagram_Yaxis_label.setFont(self.display_label_font)	
		
		self.defects_diagram_Yaxis_box = QLineEdit()						
		self.defects_diagram_Yaxis_box_layout.addWidget(self.defects_diagram_Yaxis_box)			
		
		self.defects_diagram_savefigure_button = QPushButton("Save Plot")	
		self.defects_diagram_savefigure_layout.addWidget(self.defects_diagram_savefigure_button)
	
	
	def Generate_DefectsDiagram_Plot(self, event):
		
		self.defects_diagram_plot_drawing.remove()
		self.defects_diagram_plot_drawing = self.defects_diagram_plot_figure.add_subplot(111)
		
		self.defect_plots = {}
		self.Activate_DefectsDiagram_Plot_Axes()
		
		if len(self.defects_diagram_Yaxis_box.text()) == 0:
			self.Plot_DefectsDiagram('HSE06')
			self.defects_diagram_plot_drawing.set_xlim(self.EVBM, self.ECBM)
			self.defects_diagram_plot_canvas.draw()
		else: 
			self.ymin = float(self.defects_diagram_Yaxis_box.text().split(',')[0])
			self.ymax = float(self.defects_diagram_Yaxis_box.text().split(',')[1])
			self.Plot_DefectsDiagram('HSE06')
			self.defects_diagram_plot_drawing.set_xlim(self.EVBM, self.ECBM)
			self.defects_diagram_plot_drawing.set_ylim(self.ymin, self.ymax)
			self.defects_diagram_plot_canvas.draw()
	
	
	def Plot_DefectsDiagram(self, xc):
		
		self.dmu_a = self.pressed_mu1
		self.dmu_b = self.pressed_mu2
		self.dmu_c = self.pressed_mu3
		self.dmu_d = self.mu4

		df = (self.df0[self.df0['Xc']==xc]).copy()
		df['class'] = ['$'+x+'^{'+str(y)+'}$' for x,y in zip(df['Defect'],df['q'])]
		df['En_prist'] = df['En_prist(1x1x1)']
		df['def_form_en'] = df['En(q,D)']-8*(df['En_prist'])-(df['n_a']*(df['mu_a']+float(self.dmu_a)) + df['n_b']*(df['mu_b']+float(self.dmu_b)) + df['n_c']*(df['mu_c']+float(self.dmu_c)) + df['n_d']*(df['mu_d']+float(self.dmu_d)))
		
		for idx in df.index:
			defect_name = df.ix[idx]['Defect'] + str(df.ix[idx]['q'])
			Ef_step = 0.02
			self.EVBM = df.ix[idx]['E_VBM']
			self.ECBM = df.ix[idx]['E_VBM'] + df.ix[idx]['Eg']
			self.Ef_pts = np.arange(self.EVBM, self.ECBM+0.1, Ef_step)
			self.E_enthalpy = (df.ix[idx])['def_form_en'] + ((df.ix[idx])['q']*self.Ef_pts)			

			if 'V' in df.ix[idx]['Defect']:
				if self.first_species in df.ix[idx]['Defect']:
					color = 'blue'
				elif self.second_species in df.ix[idx]['Defect']:
					color = 'purple'
				elif self.third_species in df.ix[idx]['Defect']:
					color = 'black'
				elif self.fourth_species in df.ix[idx]['Defect']:
					color = 'gray'
			else:
				color = 'orange'
			
			if df.ix[idx]['q'] == 0:
				style = '-'
			else:
				style = '--'
			
			# If defect line is there, then update it; otherwise, draw a new defect line
			try:
				self.defect_plots[defect_name].set_data(self.Ef_pts, self.E_enthalpy)
			except:
				self.defect_plots[defect_name], = self.defects_diagram_plot_drawing.plot(self.Ef_pts, self.E_enthalpy, color = color, ls = style, label=str((df.ix[idx])['class']))

		self.defects_diagram_plot_drawing.xaxis.labelpad = 2.5
		self.defects_diagram_plot_drawing.yaxis.labelpad = 2.0
		plt.subplots_adjust(left = 0.2, right = 0.85, top = 0.9, bottom = 0.1)
		plt.xticks((self.EVBM,self.ECBM),('VBM = 0 eV','CBM ='+str(df.ix[idx]['Eg'])+' eV'),fontsize='6',color='slategray')
		
		# Remove labels before redrawing them at new position
		self.defects_diagram_plot_drawing.texts.clear()
		labelLines(self.defects_diagram_plot_drawing.get_lines(),align=False,xvals=[self.EVBM+0.15, self.EVBM+0.10, self.EVBM+0.15, self.EVBM+0.12, self.EVBM+0.15, self.EVBM+0.10, self.EVBM+0.15, self.EVBM+0.02, self.EVBM+0.17, self.EVBM+0.17, self.EVBM+0.17,self.EVBM+0.17, self.EVBM+0.17, self.EVBM+0.16, self.EVBM+0.14,self.EVBM+0.12, self.EVBM+0.15, self.EVBM+0.14, self.EVBM+0.15, self.EVBM+0.12, self.EVBM+0.10, self.EVBM+0.10],fontsize=8,bbox=dict(facecolor='white',alpha=0.8,edgecolor='white',pad=0.5))
		
		try:
			self.defects_stability_region.remove()
			self.defects_stability_region = self.defects_diagram_plot_drawing.fill_between(self.Ef_pts,0,-0.8,facecolor='#614126', interpolate=True, alpha=.1)
		except:
			self.defects_stability_region = self.defects_diagram_plot_drawing.fill_between(self.Ef_pts,0,-0.8,facecolor='#614126', interpolate=True, alpha=.1)
		
		self.defects_diagram_plot_canvas.draw()
		self.defects_diagram_savefigure_button.clicked[bool].connect(self.Defects_Diagram_SaveFigure)
		return self.listdH_def	
	
	def Defects_Diagram_SaveFigure(self):

		plt.savefig(self.main_compound+'_'+str(round(self.pressed_mu1))+'_'+str(round(self.pressed_mu2))+'_'+str(round(self.pressed_mu3))+'_'+str(round(self.mu4))+'.pdf', bbox_inches='tight')

	###############################################################################################
	################################### Carrier Concentration #####################################
	###############################################################################################

	def Activate_CarrierConcentration_Plot_Axes(self):
		
		self.carrier_concentration_plot_drawing.set_xlim(200, 800)
		self.carrier_concentration_plot_drawing.set_ylim((1e16, 1e22))		
		self.carrier_concentration_plot_drawing.set_xlabel("T(K)", fontdict=self.font, labelpad=12)
		self.carrier_concentration_plot_drawing.set_ylabel("p (cm$^{-3}$)", fontdict=self.font, rotation=90, labelpad=20)
		self.carrier_concentration_plot_drawing.xaxis.tick_bottom()
		self.carrier_concentration_plot_drawing.yaxis.tick_left()
		self.carrier_concentration_plot_drawing.tick_params(axis='both', labelsize=6)
		self.carrier_concentration_plot_drawing.xaxis.set_label_position("bottom")
		self.carrier_concentration_plot_drawing.yaxis.set_label_position("left")
		plt.subplots_adjust(left = 0.2, right = 0.85, top = 0.9, bottom = 0.1)
		self.carrier_concentration_plot_drawing.xaxis.labelpad = 2.5
		self.carrier_concentration_plot_drawing.yaxis.labelpad = 2.0	
		self.carrier_concentration_plot_drawing.set_aspect("auto")

	def Activate_Generate_CarrierConcentration_Plot_Tools(self):
		
		self.generate_carrier_concentration_plot_button = QPushButton("Carrier Concentration")
		self.generate_carrier_concentration_plot_button.clicked[bool].connect(self.Generate_CarrierConcentration_Plot)
		self.generate_carrier_concentration_plot_button_layout.addWidget(self.generate_carrier_concentration_plot_button)

		self.Equilibrium_Fermi_Energy_label = QLabel(u"E"+"<sub>f</sub>"+"<sup>eq</sup>")	
		self.Equilibrium_Fermi_Energy_label.setFont(self.display_label_font)	
		self.Equilibrium_Fermi_Energy_label.setAlignment(Qt.AlignCenter)
		self.Equilibrium_Fermi_Energy_layout.addWidget(self.Equilibrium_Fermi_Energy_label)

		self.Equilibrium_Fermi_Energy_display_widget = QWidget()									
		self.Equilibrium_Fermi_Energy_display_layout = QHBoxLayout(self.Equilibrium_Fermi_Energy_widget)						
		self.Equilibrium_Fermi_Energy_display = QLCDNumber()										
		self.Equilibrium_Fermi_Energy_display.setStyleSheet("QLCDNumber {color: black}")			
		self.Equilibrium_Fermi_Energy_display_layout.addWidget(self.Equilibrium_Fermi_Energy_display)		

	def Generate_CarrierConcentration_Plot(self):
		
		self.Equilibrium_Fermi_Energy_display.display(0.0)

		self.carrier_concentration_plot_drawing.remove()
		self.carrier_concentration_plot_drawing = self.carrier_concentration_plot_figure.add_subplot(111)
		
		self.carrier_plots = {}
		self.Activate_CarrierConcentration_Plot_Axes()
		
		self.Plot_Carrier_Concentration()
		self.carrier_concentration_plot_canvas.draw()

	def Defect_formation_energy(self, ef):
		df = (self.df0).copy()

		df['class'] = ['$'+x+'^{'+str(y)+'}$' for x,y in zip(df['Defect'],df['q'])]
		df['En_prist'] = df['En_prist(1x1x1)']
		df['def_form_en'] = df['En(q,D)']-8*(df['En_prist'])-(df['n_a']*(df['mu_a']+float(self.dmu_a)) + df['n_b']*(df['mu_b']+float(self.dmu_b)) + df['n_c']*(df['mu_c']+float(self.dmu_c)) + df['n_d']*(df['mu_d']+float(self.dmu_d)))
		
		self.listdH_def = []
		for idx in df.index:
			self.E_enthalpy = (df.ix[idx])['def_form_en'] + ((df.ix[idx])['q']*float(ef))			
			self.listdH_def.append(self.E_enthalpy)
		return self.listdH_def

	def Calculate_Defect_Concentration(self, ef, T):
		df = (self.df0).copy()
		carriers = []
		self.listdH_def=self.Defect_formation_energy(ef)

		for self.E_enthalpy in self.listdH_def:
			idx = self.listdH_def.index(self.E_enthalpy)
			if df.ix[idx]['q']==0:
				continue
			elif df.ix[idx]['q']>0:
				if self.first_species in df.ix[idx]['Defect']:
					N = self.main_compound_number_first_specie/self.vol
				elif self.second_species in df.ix[idx]['Defect']:
					N = self.main_compound_number_second_specie/self.vol
				elif self.third_species in df.ix[idx]['Defect']:
					N= self.main_compound_number_third_specie/self.vol
				elif self.fourth_species in df.ix[idx]['Defect']:
					N=self.main_compound_number_fourth_specie/self.vol
				cc = N*exp(-self.E_enthalpy/(self.k*T))
				carriers.append(float(cc))

			elif df.ix[idx]['q']<0:
				if self.first_species in df.ix[idx]['Defect']:
					N = self.main_compound_number_first_specie/self.vol
				elif self.second_species in df.ix[idx]['Defect']:
					N = self.main_compound_number_second_specie/self.vol
				elif self.third_species in df.ix[idx]['Defect']:
					N = self.main_compound_number_third_specie/self.vol
				elif elemd in df.ix[idx]['Defect']:
					N = self.main_compound_number_fourth_specie/self.vol
				cc = -N*exp(-self.E_enthalpy/(self.k*T)) 
				carriers.append(float(cc))
		self.def_carriers = sum(carriers)

		return self.def_carriers

	def Calculate_Hole_Concentration(self, ef, T):
		list_gE = []
		list_fE = []
		list_en = []
		doscar=open('DOSCAR','r')
		lines=doscar.readlines()[6:]
		for line in lines:
			gE=float(line.split()[1])*10e24   #  #of states/(en x cm3)
			en=float(line.split()[0])
	
			if en <= self.EVBM:
				fE=gE*(1-1/(1+exp((en-ef)/((self.k)*T))))
				list_gE.append(gE) 
				list_fE.append(fE)
				list_en.append(en)
		p=integrate.simps(list_fE,list_en)
		return p


	def Calculate_Electron_Concentration(self, ef, T):
		list_gE = []
		list_fE = []
		list_en = []
		doscar=open('DOSCAR','r')
		lines=doscar.readlines()[6:]
		for line in lines:
			gE=float(line.split()[1])*10e24   
			en=float(line.split()[0])
			if en >= self.ECBM: 
				fE=gE/(1+exp((en-ef)/((self.k)*T)))
				list_gE.append(gE) 
				list_fE.append(fE)
				list_en.append(en)
		n=-integrate.simps(list_fE,list_en)
		return n

	def Calculate_Total_Charge_Density(self, ef, T):
		p_carriers=self.Calculate_Hole_Concentration(ef, T)
		n_carriers=self.Calculate_Electron_Concentration(ef, T)
		def_carriers=self.Calculate_Defect_Concentration(ef,T)

		self.list_carriers=[]
		self.list_carriers.append(p_carriers)
		self.list_carriers.append(n_carriers)
		self.list_carriers.append(def_carriers)
		tot_charge_density=sum(self.list_carriers)

		return tot_charge_density	


	def Calculate_Equilibrium_Fermi_Energy(self, T):
		list_tot_charge_density=[]
		ef_step=0.02 
		list_ef=[]
		ef_eq = 0.0
		for ef in np.arange(self.EVBM,self.ECBM,ef_step):
			list_ef.append(ef)
			list_tot_charge_density.append(self.Calculate_Total_Charge_Density(ef,T))

		for cdindex in np.arange(0,len(list_tot_charge_density)-1):
			if np.sign(list_tot_charge_density[cdindex])!=np.sign(list_tot_charge_density[cdindex+1]):
				ef_eq = list_ef[cdindex]

		return ef_eq

	def Plot_Carrier_Concentration(self):
		list_T=[]
		list_p=[]
		list_n=[]
		list_pn=[]
		for T in np.arange(200,801,100):
			list_T.append(T)
			ef_eq=self.Calculate_Equilibrium_Fermi_Energy(T)
			p='%10.2e'%(self.Calculate_Hole_Concentration(ef_eq, T))
			n='%10.2e'%(self.Calculate_Electron_Concentration(ef_eq, T))
			pn='%10.2e'%(self.Calculate_Hole_Concentration(ef_eq, T)+self.Calculate_Electron_Concentration(ef_eq, T))
			list_p.append(float(p))
			list_n.append(float(n))
			list_pn.append(float(pn))

		plt.semilogy(list_T,list_pn,'o-',color='orange')
		print(list_p)
		print(list_n)
		print(list_pn)
		print(ef_eq-self.EVBM)
		plt.xlabel('T(K)',fontsize='12')
		plt.ylabel('n$_i$ (cm$^{-3}$)',fontsize='14')

		eq=ef_eq-self.EVBM
		self.Equilibrium_Fermi_Energy_display.display(eq)


app = QApplication([])
aw = Quaternary_ApplicationWindow(main_compound = "Cu2HgGeTe4", first_species = "Cu", second_species = "Hg", third_species = "Ge", fourth_species = "Te")
sys.exit(app.exec_())


















