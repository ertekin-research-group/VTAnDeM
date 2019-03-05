
__author__ = 'Michael_Lidia_Jiaxing_Benita_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'


import os
import sys
import numpy as np
import pandas as pd
import itertools
import periodictable
import copy

import PyQt5
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
import matplotlib.pyplot as plt



class Ternary_Phase_Diagram(object):
	
	def __init__(self, parent = None, main_compound = None, first_species = None, second_species = None, third_species = None):
		
		#QMainWindow.__init__(self)
		
		# Font description for phase stability diagram plot
		self.font = {'family': 'Arial',
				'color':  'black',
				'weight': 'normal',
				'size': 16 }
		
		# Label the first, second, and third species of the atoms in the ternary compound
		self.first_species = first_species
		self.second_species = second_species
		self.third_species = third_species
		self.species_list = [self.first_species, self.second_species, self.third_species]					# Species list (order MAY change)
		self.species_list_radiobuttons = [self.first_species, self.second_species, self.third_species]		# Names of each radio button (order will NOT change)
		
		# Establish a variable for the main ternary compound
		self.main_compound = main_compound
		
		# Go into "Compounds_Tracker.csv" file and extract all DFT results on possible compounds
		self.compounds_info = {}		# *Note: It is important to declare the self.compounds_info variable BEFORE running self.Obtain_Compounds_Data()
		self.Obtain_Compounds_Data()
		
		# Information about main ternary compound
		self.main_compound_number_first_specie = self.compounds_info[self.main_compound][self.first_species]	# Number of first species in ternary compound
		self.main_compound_number_second_specie = self.compounds_info[self.main_compound][self.second_species]	# Number of second species in ternary compound
		self.main_compound_number_third_specie = self.compounds_info[self.main_compound][self.third_species]	# Number of third species in ternary compound
		self.main_compound_enthalpy = self.compounds_info[self.main_compound]["enthalpy"]						# Enthalpy of ternary compound
		
		# Endpoints for ternary phase diagram
		self.phasediagram_endpoints = min(self.main_compound_enthalpy/self.main_compound_number_first_specie, self.main_compound_enthalpy/self.main_compound_number_second_specie, self.main_compound_enthalpy/self.main_compound_number_third_specie)
		
		
		### *Note: The above variables did not utilize any function of PyQT yet.
		###	This is really where we start playing around with the functionalities of PyQT.
		
		# Set up layout of widgets (think of "widgets" as being like objects in the app, like a button or a plot)
		self.widgets = QWidget()						# "Main" widget. This is where all other widgets will be placed.
														# 	The reason why we need a "main" widget is because PyQT only allows the programmer
														#	to declare only ONE central widget, where everything happens. Since we need to pool
														#	into this app other widgets (like buttons and the ternary phase diagram plot),
														#	we place all of our widgets into this "main" widget.
		self.setCentralWidget(self.widgets)				# Declare self.widgets as the "main" (or central) widget.
		self.widgets_grid = QGridLayout(self.widgets)	# In the comment for self.widgets, we say that the programmer needs to place all of our
														#	widgets into this "main" widget. The way to do this is to create a GRID layout
														#	in the "main" widget. Later in the script, when we begin making widgets, we just
														#	place them into this grid using indices.
														# One needs to place widgets into this grid of the "main" widget using the command:
														#	self.widgets_grid.addWidget(widget_name, x-position, y-position, height, width)
														# 	where the x- and y- positions start at the top left.
														# Example: A widget with (x,y,w,h)=(0,0,1,2) in a 3x3 grid would look like:
														#	 _ _ _             _ _ _
														#	|     |           |O O  |
														#	|     |   ===>    |     |
														#	|_ _ _|           |_ _ _|
														# There is no need to specify the size of the grid. It will simply adjust based on how
														#	the programmer places all the widgets in the layout.
		
		
		# (WIDGET) Title of compound
		ternary_compound_title_formal = self.Compound_Name_Formal(self.main_compound)	# Generate Latex-readable version of ternary compound name
		self.ternary_compound_title = QLabel(ternary_compound_title_formal)					# QLabel is a widget that displays text
		self.ternary_compound_title.setAlignment(Qt.AlignCenter)							# Align the text to center
		self.ternary_compound_title_font = QFont("serif", 24, QFont.Bold) 					# Declare font
		self.ternary_compound_title.setFont(self.ternary_compound_title_font)				# Set the font for the QLabel text
		self.widgets_grid.addWidget(self.ternary_compound_title, 0, 0, 1, 4)				# Add the widget to the "main" widget grid layout
		
		
		# (WIDGET) Buttons to assign species' numbers (defined by function in this script)
		self.Activate_Species_Buttons(button_name = "First Species", button_x_position = 1, button_y_position = 0)
		self.Activate_Species_Buttons(button_name = "Second Species", button_x_position = 1, button_y_position = 1)
		self.Activate_Species_Buttons(button_name = "Third Species", button_x_position = 1, button_y_position = 2)
		
		
		# (WIDGET) Display ports for mu values
		self.muvalue_display_widget = QWidget()
		self.muvalue_display_layout = QHBoxLayout(self.muvalue_display_widget)
		self.widgets_grid.addWidget(self.muvalue_display_widget, 2, 0, 1, 3)
		self.Activate_MuValue_Displays()
		
		
		# (WIDGET) Generate plot button
		self.generate_plot_button_widget = QWidget()
		self.generate_plot_button_layout = QVBoxLayout(self.generate_plot_button_widget)
		self.widgets_grid.addWidget(self.generate_plot_button_widget, 3, 0, 1, 3)
		self.Activate_Generate_PhaseDiagram_Plot_Button()
		
		
		# (WIDGET) Ternary phase diagram plot
		self.ternary_plot_figure = plt.figure()
		self.ternary_plot_drawing = self.ternary_plot_figure.add_subplot(111)
		self.ternary_plot_canvas = FigureCanvas(self.ternary_plot_figure)
		self.widgets_grid.addWidget(self.ternary_plot_canvas, 0, 4, 4, 4)
		self.Activate_PhaseDiagram_Plot_Axes()
		
		
		# Initialize necessary objects in app
		self.main_compound_plot = None
		self.competing_compound_plots = {}
		self.competing_compound_colorwheel = {}
		self.PSR_vertices_plot = None
		self.PSR_vertices = []
		
		
		# Clicked point in phase diagram
		self.pressed_point = self.ternary_plot_figure.canvas.mpl_connect('button_press_event', self.Pressed_Point)
		self.pressed_point_plot, = self.ternary_plot_drawing.plot([], [], color="red", marker="o")
		self.pressed_mu1 = 0.0
		self.pressed_mu2 = 0.0
		self.pressed_mu3 = 0.0
		
		
		
		
		self.show()
	
	
	
	###############################################################################################
	############################# Obtain DFT Data of Compounds ####################################
	###############################################################################################
	
	def Obtain_Compounds_Data(self):
		
		# Import compounds info from CSV created by the script "Compounds_DataExtract.py"
		compounds_data = pd.read_csv("Compounds_Tracker.csv", header=0, index_col=0)
		header_information = list(compounds_data)
		
		for compound, compound_info in compounds_data.iterrows():						# Loop through compounds and their rows.
			if compound in self.compounds_info.keys():									# Check if the compound is already in the list.
				print("WARNING: '"+compound+"' is already in the list of compounds! This compound will thus not be recorded.")
				continue
			else:
				self.compounds_info[compound] = {}
				for header_specie, compound_data_info in zip(header_information, list(compound_info)):
					self.compounds_info[compound][header_specie] = compound_data_info	# Record number of each unique element in the compound in the list.
				for specie in self.species_list:
					if specie not in self.compounds_info[compound].keys():
						self.compounds_info[compound][specie] = 0
	
	
	###############################################################################################
	########################## Rewrite Compound Name Latex-Style ##################################
	###############################################################################################
	
	def Compound_Name_Formal(self, compound_name):
		
		# Go into the compounds_info dictionary and obtains the chemistry and stoichiometry of the compound of choice
		compound_species_info = self.compounds_info[compound_name]
		
		compound_name_formal = ""
		for species in self.species_list:				# Loop through the list of possible species that can be contained in the compound
			if compound_species_info[species] == 0:
				continue								# Don't add the species to the name if the compound doesn't contain the species
			elif compound_species_info[species] == 1:
				compound_name_formal += species			# Add the species to the name of the compound
			elif compound_species_info[species] > 1:
				compound_name_formal += species+"<sub>"+str(int(compound_species_info[species]))+"</sub>"	# Add the species to the name of the compound
																											#	with a subscript for the stoichiometry
		
		return compound_name_formal
	
	
	###############################################################################################
	#################################### Species Buttons ##########################################
	###############################################################################################
	
	def Activate_Species_Buttons(self, button_name, button_x_position, button_y_position):
		# *Note: The button name can only be one of: "First_Species", "Second_Species", or "Third_Species".
		# *Note: The button positions (x and y) correspond to the position on the "main" widget's grid layout (self.widgets_grid).
		
		species_radiobutton_widget = QWidget()									# Create a widget for the buttons of the three species. Each button will be a widget
																				#	on its own, so this will be the "sub-main" widget for each species' button.
		species_radiobutton_layout = QVBoxLayout(species_radiobutton_widget)	# Create a layout for this "sub-main" widget. When adding the species' button
																				# 	widgets to this "sub-main" widget, QVBoxLayout automatically places them on
																				#	top of one another.
		
		species_radiobutton_label = QLabel(button_name)					# Create text widget
		species_radiobutton_label.setAlignment(Qt.AlignCenter)			# Align text
		species_radiobutton_layout.addWidget(species_radiobutton_label)	# Add widget to "sub-main" widget layout
		
		first_species_button = QRadioButton(self.species_list_radiobuttons[0])										# Create the first species' button
		first_species_button.toggled.connect(lambda:self.Update_Species_Button(first_species_button, button_name))	# Connect the button to the function Update_Species_Button
		species_radiobutton_layout.addWidget(first_species_button)													# Add widget to "sub-main" widget layout
		
		second_species_button = QRadioButton(self.species_list_radiobuttons[1])										# Create the second species' button
		second_species_button.toggled.connect(lambda:self.Update_Species_Button(second_species_button, button_name))
		species_radiobutton_layout.addWidget(second_species_button)
		
		third_species_button = QRadioButton(self.species_list_radiobuttons[2])										# Create the third species' button
		third_species_button.toggled.connect(lambda:self.Update_Species_Button(third_species_button, button_name))
		species_radiobutton_layout.addWidget(third_species_button)
		
		# Initially (only activated when the app starts), one of the three buttons should be pressed
		if button_name == "First Species":
			first_species_button.setChecked(True)
		elif button_name == "Second Species":
			second_species_button.setChecked(True)
		elif button_name == "Third Species":
			third_species_button.setChecked(True)
		
		# Add the radio button "sub-main" widget to the "main" widget's grid layout
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
		
		# Reorder the list of the species
		self.species = [self.first_species, self.second_species, self.third_species]
		# *Note: 	Clicking a radio button will change the order of the species, but this alone will NOT show any visible changes
		#			in the plots. These changes will be visible in the plots when the user clicks the "Generate Plot" button.
	
	
	
	###############################################################################################
	################################## Mu Values Displays #########################################
	###############################################################################################
	
	def Activate_MuValue_Displays(self):
		
		# Create widgets for the display ports of the mu values. These mu values will change based on where on the ternary phase
		#	diagram plot the user clicks.
		# *Note: The "sub-main" widget that will contain these mu display ports for each species is already created in __init__ as
		#	the variable self.muvalue_display_widget. It's layout is self.muvalue_display_layout, and this is where each display
		#	port will be placed.
		
		self.mu_display_label_font = QFont("Arial", 12)						# Define the font
		
		self.mu1_display_widget = QWidget()									# Create a widget for the display port
		self.mu1_display_layout = QVBoxLayout(self.mu1_display_widget)		# Create a vertical layout for the widget (should only contain the label and the mu4 number, stacked on top of each other)
		self.mu1_display_label = QLabel(u"\u0394"+"\u03BC"+"<sub>a</sub>")	# Create a label widget to display the "\Delta\mu_1" text above the actual mu1 value (unicode format)
		self.mu1_display_label.setFont(self.mu_display_label_font)			# Set the font
		self.mu1_display_label.setAlignment(Qt.AlignCenter)					# Align the label to center
		self.mu1_display_layout.addWidget(self.mu1_display_label)			# Add the label
		self.mu1_display = QLCDNumber()										# Create the actual port that will display the current mu1 value
		self.mu1_display.setStyleSheet("QLCDNumber {color: black}")			# Color
		self.mu1_display_layout.addWidget(self.mu1_display)					# Add the display port
		
		self.mu2_display_widget = QWidget()
		self.mu2_display_layout = QVBoxLayout(self.mu2_display_widget)
		self.mu2_display_label = QLabel(u"\u0394"+"\u03BC"+"<sub>b</sub>")
		self.mu2_display_label.setFont(self.mu_display_label_font)
		self.mu2_display_label.setAlignment(Qt.AlignCenter)
		self.mu2_display_layout.addWidget(self.mu2_display_label)
		self.mu2_display = QLCDNumber()
		self.mu2_display.setStyleSheet("QLCDNumber {color: black}")
		self.mu2_display_layout.addWidget(self.mu2_display)
		
		self.mu3_display_widget = QWidget()
		self.mu3_display_layout = QVBoxLayout(self.mu3_display_widget)
		self.mu3_display_label = QLabel(u"\u0394"+"\u03BC"+"<sub>c</sub>")
		self.mu3_display_label.setFont(self.mu_display_label_font)
		self.mu3_display_label.setAlignment(Qt.AlignCenter)
		self.mu3_display_layout.addWidget(self.mu3_display_label)
		self.mu3_display = QLCDNumber()
		self.mu3_display.setStyleSheet("QLCDNumber {color: black}")
		self.mu3_display_layout.addWidget(self.mu3_display)
		
		self.muvalue_display_layout.addWidget(self.mu1_display_widget)
		self.muvalue_display_layout.addWidget(self.mu2_display_widget)
		self.muvalue_display_layout.addWidget(self.mu3_display_widget)
	
	
	###############################################################################################
	##################################### Phase Diagram ###########################################
	###############################################################################################
	
	def Activate_PhaseDiagram_Plot_Axes(self):
		
		# This function simply creates the axes of the ternary phase diagram.
		
		self.ternary_plot_drawing.set_xlim(self.phasediagram_endpoints, 0.0)
		self.ternary_plot_drawing.set_ylim(self.phasediagram_endpoints, 0.0)
		self.ternary_plot_drawing.set_xlabel("$\Delta\mu_{a}$ (eV)",fontdict=self.font,labelpad=12)
		self.ternary_plot_drawing.set_ylabel("$\Delta\mu_{b}$ (eV)",fontdict=self.font,rotation=270,labelpad=20)
		self.ternary_plot_drawing.xaxis.tick_top()
		self.ternary_plot_drawing.yaxis.tick_right()
		self.ternary_plot_drawing.tick_params(axis='both', labelsize=6)
		self.ternary_plot_drawing.xaxis.set_label_position("top")
		self.ternary_plot_drawing.yaxis.set_label_position("right")
		self.ternary_plot_drawing.spines['left'].set_visible(False)
		self.ternary_plot_drawing.spines['bottom'].set_visible(False)
		self.ternary_plot_drawing.set_aspect("equal")
	
	
	def Reset_PhaseDiagram_Plot_Axes(self):
		
		# This function resets the axes of the ternary phase diagram. This will be used when the 
		#	user changes the species to be plotted, since the phase stability diagram (the triangular
		#	shape) is not always the same between different species.
		
		self.ternary_plot_drawing.set_xlim(self.phasediagram_endpoints, 0.0)
		self.ternary_plot_drawing.set_ylim(self.phasediagram_endpoints, 0.0)
		self.ternary_plot_drawing.set_xlabel("$\Delta\mu_{"+self.first_species+"}$ (eV)",fontdict=self.font,labelpad=12)
		self.ternary_plot_drawing.set_ylabel("$\Delta\mu_{"+self.second_species+"}$ (eV)",fontdict=self.font,rotation=270,labelpad=20)
		self.ternary_plot_drawing.xaxis.tick_top()
		self.ternary_plot_drawing.yaxis.tick_right()
		self.ternary_plot_drawing.tick_params(axis='both', labelsize=6)
		self.ternary_plot_drawing.xaxis.set_label_position("top")
		self.ternary_plot_drawing.yaxis.set_label_position("right")
		self.ternary_plot_drawing.spines['left'].set_visible(False)
		self.ternary_plot_drawing.spines['bottom'].set_visible(False)
		self.ternary_plot_drawing.set_aspect("equal")
	
	
	def Activate_Generate_PhaseDiagram_Plot_Button(self):
		
		# The user is given a button to generate the ternary phase diagram plot. This is because
		#	when the user changes the species to be plotted, the program does not automatically
		#	plot the new phase diagram.
		
		self.generate_plot_button_label = QLabel("Hello. Welcome to [insert name of app].")	# A heart-warming welcome message that is only visible at the start of the app
		self.generate_plot_button_label.setAlignment(Qt.AlignCenter)
		self.generate_plot_button_layout.addWidget(self.generate_plot_button_label)
		
		self.generate_plot_button = QPushButton("Generate Plot!")							# Create the button widget that the user will press to generate a new phase diagram
		self.generate_plot_button.setCheckable(True)
		self.generate_plot_button.clicked[bool].connect(self.Generate_PhaseDiagram_Plot)
		self.generate_plot_button_layout.addWidget(self.generate_plot_button)
	

	def Generate_PhaseDiagram_Plot(self, event):

		# This function specifies what happens when the user clicks the "Generate Plot" button. The 
		#	only thing it needs to check is whether the elements chosen as the first, second, and third 
		#	species are unique.
		
		if len(self.species) > len(set(self.species)):	# Every time someone clicks the species' buttons, the self.species list gets updated.
														#	The "set" function checks the self.species list and omits any that are repeated.
														#	If any are repeated, then the chosen species are not unique.
			self.generate_plot_button_label.setText("Yo... Pick UNIQUE elements!")	# Maybe change the phrasing ;)
		
		else:
			
			self.generate_plot_button_label.setText("")
			
			# Update the number of each species in the main compound based on the above change
			self.main_compound_number_first_specie = self.compounds_info[self.main_compound][self.first_species]
			self.main_compound_number_second_specie = self.compounds_info[self.main_compound][self.second_species]
			self.main_compound_number_third_specie = self.compounds_info[self.main_compound][self.third_species]
			
			# Reset the mu value displays
			self.mu1_display_label.setText(u"\u0394\u03BC<sub>"+self.first_species+"</sub>")
			self.mu2_display_label.setText(u"\u0394\u03BC<sub>"+self.second_species+"</sub>")
			self.mu3_display_label.setText(u"\u0394\u03BC<sub>"+self.third_species+"</sub>")
			self.mu1_display.display(0.0)
			self.mu2_display.display(0.0)
			self.mu3_display.display(0.0)
			
			# Set up the new plot
			self.ternary_plot_drawing.remove()
			self.ternary_plot_drawing = self.ternary_plot_figure.add_subplot(111)
			
			# Reset the plots of the main and competing compounds
			self.main_compound_plot = None
			self.competing_compound_plots = {}
			
			# Reset the axes of the plot
			self.Reset_PhaseDiagram_Plot_Axes()
			
			# Plot the phase stability diagram of the ternary compound using the new settings
			self.Plot_Ternary_PhaseDiagram()
			
			# Reset clicked point
			self.pressed_point_plot, = self.ternary_plot_drawing.plot([], [], color="red", marker="o")
	
	
	def Plot_Ternary_PhaseDiagram(self):
		
		"""
		if len(self.species) != 3:		# Check if there are three species for ternary.
			print "Your species list currently does not support ternary compounds. Either update with Add_Species() or pick a different plotting scheme."
			return
		
		if (first_species not in self.species) or (second_species not in self.species):		# Check if first and second species are in the species list.
			print "One of the species is not a valid species. We recommend using the Add_Species() command to add it."
			return
		
		third_species = [specie for specie in self.species if (specie != first_species) and (specie != second_species)][0]
		"""
		
		main_compound_deltamu_first_species = np.linspace(self.main_compound_enthalpy/self.main_compound_number_first_specie, 0, 1000)
		
		main_compound_stability_limit = (self.main_compound_enthalpy - self.main_compound_number_first_specie * main_compound_deltamu_first_species) / self.main_compound_number_second_specie
		
		
		try:
			self.main_compound_plot.set_ydata(main_compound_stability_limit)
		except:
			self.main_compound_plot, = self.ternary_plot_drawing.plot(main_compound_deltamu_first_species, main_compound_stability_limit)
		
		
		
		stability_minimum_bound = []
		stability_maximum_bound = []
		
		vertical_left_values = []
		vertical_right_values = []
		
		
		for competing_compound_index, competing_compound in enumerate(self.compounds_info.keys()):
			
			foreign_compound = False
			for element in periodictable.elements:
				if str(element) in self.species:
					continue
				elif str(element) not in self.compounds_info[competing_compound].keys():
					continue
				else:
					if self.compounds_info[competing_compound][str(element)] != 0.0:
						foreign_compound = True
			
			
			if foreign_compound:
				continue
			elif (competing_compound in self.species) or (competing_compound == self.main_compound):
				continue
			else:
				#print "Analyzing the competing compound '"+competing_compound+"'..."
				
				competing_compound_number_first_specie = self.compounds_info[competing_compound][self.first_species]
				competing_compound_number_second_specie = self.compounds_info[competing_compound][self.second_species]
				competing_compound_number_third_specie = self.compounds_info[competing_compound][self.third_species]
				
				competing_compound_enthalpy = self.compounds_info[competing_compound]["enthalpy"]
				
				coefficient_first_specie = competing_compound_number_first_specie - float(competing_compound_number_third_specie)*float(self.main_compound_number_first_specie)/float(self.main_compound_number_third_specie)
				coefficient_second_specie = competing_compound_number_second_specie - float(competing_compound_number_third_specie)*float(self.main_compound_number_second_specie)/float(self.main_compound_number_third_specie)
				
				competing_compound_deltamu_first_species = copy.deepcopy(main_compound_deltamu_first_species)
				competing_compound_deltamu_second_species = ( competing_compound_enthalpy - float(competing_compound_number_third_specie)/float(self.main_compound_number_third_specie) * self.main_compound_enthalpy \
															- (coefficient_first_specie * main_compound_deltamu_first_species) ) / coefficient_second_specie
				
				if coefficient_second_specie > 0.0:
					stability_maximum_bound.append(competing_compound_deltamu_second_species)
				elif coefficient_second_specie < 0.0:
					stability_minimum_bound.append(competing_compound_deltamu_second_species)
				
				competing_compound_deltamu_first_species_limit = [competing_compound_deltamu_first_species[i] for i in range(len(main_compound_deltamu_first_species)) if (main_compound_stability_limit[i] < competing_compound_deltamu_second_species[i])]
				competing_compound_deltamu_second_species_limit = [competing_compound_deltamu_second_species[i] for i in range(len(main_compound_deltamu_first_species)) if (main_compound_stability_limit[i] < competing_compound_deltamu_second_species[i])]
				
				try:
					self.competing_compound_plots[competing_compound].set_data(competing_compound_deltamu_first_species_limit, competing_compound_deltamu_second_species_limit)
				except:
					# If the plot doesn't exist initially, then create it
					self.competing_compound_plots[competing_compound], = self.ternary_plot_drawing.plot(competing_compound_deltamu_first_species_limit, competing_compound_deltamu_second_species_limit, label=competing_compound)
					self.competing_compound_colorwheel[competing_compound] = self.competing_compound_plots[competing_compound].get_color()
		
		# Check if the legend for the phase diagram has been activated
		try:
			self.ternary_plot_legend_widget
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
			if stability_absolute_minimum[i] < stability_absolute_maximum[i]:
				self.main_compound_deltamu_first_species_cutoff.append(main_compound_deltamu_first_species[i])
				self.stability_minimum_cutoff.append(stability_absolute_minimum[i])
				self.stability_maximum_cutoff.append(stability_absolute_maximum[i])
		
		try:
			self.phase_stability_region.remove()
			self.phase_stability_region = self.ternary_plot_drawing.fill_between(self.main_compound_deltamu_first_species_cutoff, self.stability_maximum_cutoff, self.stability_minimum_cutoff, facecolor='0.75')
		except:
			self.phase_stability_region = self.ternary_plot_drawing.fill_between(self.main_compound_deltamu_first_species_cutoff, self.stability_maximum_cutoff, self.stability_minimum_cutoff, facecolor='0.75')
		
		
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
					self.PSR_vertices_plot = self.ternary_plot_drawing.scatter(*zip(*PSR_Vertices_Unrepeated), s=20, c='black')
				except:
					self.PSR_vertices_plot = self.ternary_plot_drawing.scatter(*zip(*PSR_Vertices_Unrepeated), s=20, c='black')
					pass
		elif self.phase_stability_region.get_paths() == []:
			self.PSR_vertices = []
			try:
				self.PSR_vertices_plot.remove()
			except:
				pass
		
		
		
		
		self.ternary_plot_figure.canvas.draw_idle()
	
	
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
			self.pressed_mu3 = (self.main_compound_enthalpy - self.main_compound_number_first_specie*self.pressed_mu1 - self.main_compound_number_second_specie*self.pressed_mu2) / self.main_compound_number_third_specie
			
			# Update the display ports for the mu values
			self.mu1_display.display(self.pressed_mu1)
			self.mu2_display.display(self.pressed_mu2)
			self.mu3_display.display(self.pressed_mu3)
			
			# Update the red dot where the user clicked on the phase diagram
			self.pressed_point_plot.set_data([self.pressed_mu1], [self.pressed_mu2])
			self.ternary_plot_figure.canvas.draw_idle()
	
	
	
	
	
	###############################################################################################
	######################################### Legend ##############################################
	###############################################################################################
	
	def Activate_Legend(self):
		
		self.ternary_plot_legend_widget = QWidget()
		self.ternary_plot_legend_layout = QGridLayout(self.ternary_plot_legend_widget)
		
		ternary_plot_legend_name = QLabel("Competing\nCompound\nColors")
		ternary_plot_legend_name.setAlignment(Qt.AlignCenter)
		self.ternary_plot_legend_layout.addWidget(ternary_plot_legend_name, 0, 0, 1, 2)
		
		for competing_compound_index, competing_compound_name, competing_compound_color in zip(range(len(self.competing_compound_colorwheel)), self.competing_compound_colorwheel.keys(), self.competing_compound_colorwheel.values()):
			
			competing_compound_name_label = QLabel(self.Compound_Name_Formal(competing_compound_name))
			competing_compound_name_label.setAlignment(Qt.AlignCenter)
			self.ternary_plot_legend_layout.addWidget(competing_compound_name_label, competing_compound_index+1, 0)
			
			competing_compound_color_label = QWidget()
			#competing_compound_color_label.resize(1, 100)
			competing_compound_color_label.setAutoFillBackground(True)
			competing_compound_color_palette = competing_compound_color_label.palette()
			color = QColor(competing_compound_color)
			competing_compound_color_palette.setColor(competing_compound_color_label.backgroundRole(), color)
			competing_compound_color_label.setPalette(competing_compound_color_palette)
			self.ternary_plot_legend_layout.addWidget(competing_compound_color_label, competing_compound_index+1, 1)
		
		self.widgets_grid.addWidget(self.ternary_plot_legend_widget, 2, 3, 2, 1)






app = QApplication([])
aw = Ternary_ApplicationWindow(main_compound = "CuGaTe2", first_species = "Cu", second_species = "Ga", third_species = "Te")
sys.exit(app.exec_())

















