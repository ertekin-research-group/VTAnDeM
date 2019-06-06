
__author__ = 'Michael_Lidia_Jiaxing_Benita_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'

import numpy as np
import matplotlib.pyplot as plt

import PyQt5
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

from obtain_DFT_data import *
from ternary_phase_diagram import Ternary_Phase_Diagram
from ternary_defects_diagram import Ternary_Defects_Diagram



class Ternary_Main_VTAnDeM_Window(QMainWindow):
	
	def __init__(self, parent = None, main_compound = None, first_species = None, second_species = None, third_species = None):
		
		# Inherit all initial variables from the QMainWindow class
		QMainWindow.__init__(self)
		
		
		# Font description for phase stability diagram plot
		self.font = {'family': 'Arial',
				'color':  'black',
				'weight': 'normal',
				'size': 16 }
		
		
		# Set up the framework of the application window (including file menu, exit function, etc.)
		self.Setup_Window_Framework()
		
		
		# Establish a variable for the main ternary compound
		self.main_compound = main_compound
		
		
		# Label the first, second, and third species of the atoms in the ternary compound
		self.first_species = first_species
		self.second_species = second_species
		self.third_species = third_species
		self.species_list = [self.first_species, self.second_species, self.third_species]					# Species list (order MAY change)
		self.species_list_radiobuttons = [self.first_species, self.second_species, self.third_species]		# Names of each radio button (order will NOT change)
		
		
		# Keep track of mu values of the species in the ternary compound (will be updated as user uses mu4 slider)
		self.mu_values = {}
		self.mu_values[first_species] = 0.0
		self.mu_values[second_species] = 0.0
		self.mu_values[third_species] = 0.0
		
		
		# Obtain DFT data
		self.compounds_info = Obtain_DFT_Data(self)	# Total energies/enthalpies for phase diagram
		self.df0 = Obtain_Defects_Data(self)		# Defect energies for defects diagram
		
		
		# Information about main ternary compound
		self.main_compound_number_first_specie = self.compounds_info[self.main_compound][self.first_species]	# Number of first species in ternary compound
		self.main_compound_number_second_specie = self.compounds_info[self.main_compound][self.second_species]	# Number of second species in ternary compound
		self.main_compound_number_third_specie = self.compounds_info[self.main_compound][self.third_species]	# Number of third species in ternary compound
		self.main_compound_enthalpy = self.compounds_info[self.main_compound]["enthalpy"]						# Enthalpy of ternary compound
		
		
		
		# Set up ternary phase diagram object
		self.PhaseDiagram = Ternary_Phase_Diagram(self, main_compound = main_compound, first_species = first_species, second_species = second_species, third_species = third_species)
		self.PhaseDiagram.compounds_info = self.compounds_info
		self.PhaseDiagram.Activate_PhaseDiagram_Object()
		self.PhaseDiagram.Activate_PhaseDiagram_Plot_Axes()
		
		
		# Set up defects diagram object
		self.DefectsDiagram = Ternary_Defects_Diagram(self, main_compound = main_compound, first_species = first_species, second_species = second_species, third_species = third_species)
		self.DefectsDiagram.df0 = self.df0
		self.DefectsDiagram.Activate_DefectsDiagram_Plot_Axes()
		
		
		#######################################
		####### Define and place widgets ######
		#######################################
		
		# Set up layout of widgets (think of "widgets" as being like objects in the app, like a button or a plot)
		self.widgets = QWidget()						# "Main" widget. This is where all other widgets will be placed.
														# 	The reason why we need a "main" widget is because PyQT only allows the programmer
														#	to declare only ONE central widget, where everything happens. Since we need to pool
														#	into this app other widgets (like buttons and the ternary phase diagram plot),
														#	we place all of our widgets into this "main" widget.
		self.setCentralWidget(self.widgets)				# Declare self.widgets as the "main" (or central) widget.
		self.widgets_grid = QHBoxLayout(self.widgets)	# The layout of the main widget should be such that all widgets placed inside the main
														#	widget (i.e. buttons, plots, etc.) are placed horizontally.
		
		
		###### "Settings" Widget
		self.settings_widget = QWidget()									# One of the main sub-widgets is where the user defines the settings of the plots.
		self.settings_widget_layout = QVBoxLayout(self.settings_widget)		# The settings should be placed on top of one another, i.e. vertically.
		
		# (WIDGET) Title of compound
		ternary_compound_title_formal = self.Compound_Name_Formal(self.main_compound)	# Generate Latex-readable version of ternary compound name
		self.ternary_compound_title = QLabel(ternary_compound_title_formal)			# QLabel is a widget that displays text
		self.ternary_compound_title.setAlignment(Qt.AlignCenter)							# Align the text to center
		self.ternary_compound_title_font = QFont("serif", 24, QFont.Bold) 				# Declare font
		self.ternary_compound_title.setFont(self.ternary_compound_title_font)			# Set the font for the QLabel text
		self.settings_widget_layout.addWidget(self.ternary_compound_title)				# Add the widget to the "main" widget grid layout
		
		# (WIDGET) Buttons to assign species' numbers
		self.species_buttons_widget = QWidget()
		self.species_buttons_widget_layout = QHBoxLayout(self.species_buttons_widget)
		for species_number in range(len(self.species_list_radiobuttons)):
			self.species_buttons_widget_layout.addWidget(self.Generate_Species_Buttons(button_name = "Species "+str(int(species_number)+1)))
		self.settings_widget_layout.addWidget(self.species_buttons_widget)
		
		# (WIDGET) Display ports for mu values
		self.muvalue_display_widget = QWidget()
		self.muvalue_display_layout = QVBoxLayout(self.muvalue_display_widget)
		self.settings_widget_layout.addWidget(self.muvalue_display_widget)
		self.Activate_MuValue_Displays()
		
		# (WIDGET) Button to generate phase diagram
		self.generate_ternary_phase_diagram_plot_button_widget = QPushButton("Generate Phase Diagram")
		self.generate_ternary_phase_diagram_plot_button_widget.clicked[bool].connect(self.Generate_PhaseDiagram_Plot_Function)
		self.settings_widget_layout.addWidget(self.generate_ternary_phase_diagram_plot_button_widget)
		
		# (WIDGET) Button to generate defects diagram
		self.generate_defects_diagram_plot_button_widget = QPushButton("Generate Defect Diagram")
		self.generate_defects_diagram_plot_button_widget.clicked[bool].connect(self.Generate_DefectsDiagram_Plot_Function)
		self.settings_widget_layout.addWidget(self.generate_defects_diagram_plot_button_widget)
		
		# (WIDGET) A heartwarming message
		self.settings_message_widget = QWidget()
		self.settings_message_widget_layout = QHBoxLayout(self.settings_message_widget)
		self.constant_text = QLabel("Message:")
		self.constant_text.setAlignment(Qt.AlignRight)
		self.settings_message_widget_layout.addWidget(self.constant_text)
		self.settings_message = QLabel("Hello. Welcome to VTAnDeM.")
		self.settings_message.setAlignment(Qt.AlignLeft)
		self.settings_message_widget_layout.addWidget(self.settings_message)
		self.settings_widget_layout.addWidget(self.settings_message_widget)
		
		# Add settings widget to the main widget
		self.widgets_grid.addWidget(self.settings_widget)
		
		
		###### "Tabs" Widget
		self.plot_tabs_widget = QTabWidget()	# This is the other main sub-widget, which contains all the viewports for each relevant plot.
		
		# Initialize ternary phase diagram and defects diagram
		self.ternary_phase_diagram_plot = self.PhaseDiagram.ternary_phase_diagram_plot_canvas
		self.defects_diagram_plot = self.DefectsDiagram.ternary_defects_diagram_plot_canvas
		
		
		
		### First tab
		self.tab1 = QWidget()
		self.tab1_layout = QHBoxLayout(self.tab1)
		
		# Phase diagram in first tab
		self.tab1_phasediagram_widget = QWidget()
		self.tab1_phasediagram_widget_layout = QVBoxLayout(self.tab1_phasediagram_widget)
		
		# (WIDGET) Quaternary phase diagram plot
		self.tab1_phasediagram_widget_layout.addWidget(self.ternary_phase_diagram_plot)
		
		# Save phase diagram as figure
		self.phase_diagram_savefigure_button = QPushButton("Save Phase Diagram Figure")
		self.phase_diagram_savefigure_button.clicked[bool].connect(lambda: self.SaveFigure("Phase Diagram"))
		self.tab1_phasediagram_widget_layout.addWidget(self.phase_diagram_savefigure_button)
		
		self.tab1_layout.addWidget(self.tab1_phasediagram_widget)
		
		# Defects diagram in first tab
		self.tab1_defectsdiagram_widget = QWidget()
		self.tab1_defectsdiagram_widget_layout = QVBoxLayout(self.tab1_defectsdiagram_widget)
		
		# (WIDGET) Defects diagram plot
		self.tab1_defectsdiagram_widget_layout.addWidget(self.defects_diagram_plot)
		
		self.defectsdiagram_viewport = QWidget()
		self.defectsdiagram_viewport_layout = QHBoxLayout(self.defectsdiagram_viewport)
		
		### Y-axis limits (label)
		self.defects_diagram_Yaxis_label = QLabel(u"y"+"<sub>min</sub>"+','+u"y"+"<sub>max</sub>")			
		self.defects_diagram_Yaxis_label.setAlignment(Qt.AlignRight)
		self.defectsdiagram_viewport_layout.addWidget(self.defects_diagram_Yaxis_label)
		
		### Y-axis limits (user dialog)
		self.defects_diagram_Yaxis_box = QLineEdit()
		self.defectsdiagram_viewport_layout.addWidget(self.defects_diagram_Yaxis_box)
		
		self.tab1_defectsdiagram_widget_layout.addWidget(self.defectsdiagram_viewport)
		
		
		### Save defects diagram as figure
		self.defects_diagram_savefigure_button = QPushButton("Save Defects Diagram Figure")
		self.defects_diagram_savefigure_button.clicked[bool].connect(lambda: self.SaveFigure("Defects Diagram"))
		self.tab1_defectsdiagram_widget_layout.addWidget(self.defects_diagram_savefigure_button)
		
		self.tab1_layout.addWidget(self.tab1_defectsdiagram_widget)
		
		
		self.plot_tabs_widget.addTab(self.tab1, "Both")
		
		
		
		
		
		
		
		self.widgets_grid.addWidget(self.plot_tabs_widget)
		
		self.pressed_point = self.PhaseDiagram.ternary_phase_diagram_plot_figure.canvas.mpl_connect('button_press_event', self.Pressed_Point)
		self.pressed_point_plot, = self.PhaseDiagram.ternary_phase_diagram_plot_drawing.plot([], [], color="red", marker="o")
		
		
		self.showFullScreen()
	
	
	
	
	
	
	###############################################################################################
	################################# Main Window Framework #######################################
	###############################################################################################
	
	def Setup_Window_Framework(self):
		
		# Set up an option where the user can open a new window (not functionalized yet)
		newappAction = QAction('New', self)
		newappAction.setStatusTip('Launch New Application (Tips by Michael Toriyama TM)')
		
		# Set up an option where the user can import a CSV file (not functionalized yet)
		importMenu = QMenu('Import', self)
		importAction_csv = QAction('&Import CSV', self)
		importAction_other = QAction('&Import other', self)
		importAction_csv.setStatusTip('Import CSV (Tips by Michael Toriyama TM)')
		importAction_other.setStatusTip('Import Other Stuff (Tips by Michael Toriyama TM)')
		importMenu.addAction(importAction_csv)
		importMenu.addAction(importAction_other)
		
		# Set up an option where the user can exit the current window
		exitAction = QAction('Exit', self)
		exitAction.setStatusTip('Exit Application (Tips by Michael Toriyama TM)')
		exitAction.triggered.connect(qApp.quit)
		
		# Create a menu bar
		menubar = self.menuBar()	# This menu bar automatically gets added to the existing main window
		fileMenu = menubar.addMenu('&File')
		
		# Add options (created above) to the menu bar, in order
		fileMenu.addAction(newappAction)
		fileMenu.addMenu(importMenu)
		fileMenu.addAction(exitAction)
		
		# Set the title of the window
		mainwindow_title = "VTAnDeM: Visualization Toolkit for the Analysis of Defects in Materials"
		self.setWindowTitle(mainwindow_title)
	
	
	
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
	
	def Generate_Species_Buttons(self, button_name):
		# *Note: The button name can only be one of: "First_Species", "Second_Species", or "Third_Species".
		# *Note: The button positions (x and y) correspond to the position on the "main" widget's grid layout (self.widgets_grid).
		
		species_radiobutton_widget = QWidget()									# Create a widget for the buttons of the three species. Each button will be a widget
																				#	on its own, so this will be the "sub-main" widget for each species' button.
		species_radiobutton_layout = QVBoxLayout(species_radiobutton_widget)	# Create a layout for this "sub-main" widget. When adding the species' button
																				# 	widgets to this "sub-main" widget, QVBoxLayout automatically places them on
																				#	top of one another.
		species_radiobutton_layout.setSpacing(0)
		
		species_radiobutton_label = QLabel(button_name)					# Create text widget
		species_radiobutton_label.setMargin(0)
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
		if "1" in button_name:
			first_species_button.setChecked(True)
		elif "2" in button_name:
			second_species_button.setChecked(True)
		elif "3" in button_name:
			third_species_button.setChecked(True)
		
		
		# Add the radio button "sub-main" widget to the "main" widget's grid layout
		species_radiobutton_widget.setStyleSheet("background-color: rgb(255,255,255)")
		self.species_buttons_widget_layout.addWidget(species_radiobutton_widget)
		
		return species_radiobutton_widget
	
	
	def Update_Species_Button(self, species_button, button_name):
		
		# Basically changes the species that is called "first species", "second species", etc.
		if species_button.isChecked() == True:
			if button_name == "Species 1":
				self.first_species = species_button.text()
			elif button_name == "Species 2":
				self.second_species = species_button.text()
			elif button_name == "Species 3":
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
		self.mu1_display_layout = QHBoxLayout(self.mu1_display_widget)		# Create a vertical layout for the widget (should only contain the label and the mu4 number, stacked on top of each other)
		self.mu1_display_label = QLabel(u"\u0394"+"\u03BC"+"<sub>a</sub>")	# Create a label widget to display the "\Delta\mu_1" text above the actual mu1 value (unicode format)
		self.mu1_display_label.setFont(self.mu_display_label_font)			# Set the font
		self.mu1_display_label.setAlignment(Qt.AlignCenter)					# Align the label to center
		self.mu1_display_layout.addWidget(self.mu1_display_label)			# Add the label
		self.mu1_display = QLineEdit("0.00")								# Create the actual port that will display the current mu1 value
		self.mu1_display.editingFinished.connect(lambda: self.Update_MuValue_Displays(1))
		self.mu1_display_layout.addWidget(self.mu1_display)					# Add the display port
		
		self.mu2_display_widget = QWidget()
		self.mu2_display_layout = QHBoxLayout(self.mu2_display_widget)
		self.mu2_display_label = QLabel(u"\u0394"+"\u03BC"+"<sub>b</sub>")
		self.mu2_display_label.setFont(self.mu_display_label_font)
		self.mu2_display_label.setAlignment(Qt.AlignCenter)
		self.mu2_display_layout.addWidget(self.mu2_display_label)
		self.mu2_display = QLineEdit("0.00")
		self.mu2_display.editingFinished.connect(lambda: self.Update_MuValue_Displays(2))
		self.mu2_display_layout.addWidget(self.mu2_display)
		
		self.mu3_display_widget = QWidget()
		self.mu3_display_layout = QHBoxLayout(self.mu3_display_widget)
		self.mu3_display_label = QLabel(u"\u0394"+"\u03BC"+"<sub>c</sub>")
		self.mu3_display_label.setFont(self.mu_display_label_font)
		self.mu3_display_label.setAlignment(Qt.AlignCenter)
		self.mu3_display_layout.addWidget(self.mu3_display_label)
		self.mu3_display = QLineEdit("0.00")
		self.mu3_display.setEnabled(False)
		self.mu3_display.setStyleSheet("""QLineEdit { background-color: white; color: black }""")
		self.mu3_display_layout.addWidget(self.mu3_display)
		
		self.muvalue_display_layout.addWidget(self.mu1_display_widget)
		self.muvalue_display_layout.addWidget(self.mu2_display_widget)
		self.muvalue_display_layout.addWidget(self.mu3_display_widget)
	
	
	def Update_MuValue_Displays(self, display_number):
		
		if display_number == 1:
			try:
				float(self.mu1_display.text())
			except:
				self.mu1_display.setText("0.0")
				pass
			if float(self.mu1_display.text()) > 0.0:
				self.mu1_display.setText("0.0")
			self.mu_values[self.first_species] = float(self.mu1_display.text())
		elif display_number == 2:
			try:
				float(self.mu2_display.text())
			except:
				self.mu2_display.setText("0.0")
				pass
			if float(self.mu2_display.text()) > 0.0:
				self.mu2_display.setText("0.0")
			self.mu_values[self.second_species] = float(self.mu2_display.text())
		
		self.mu_values[self.third_species] = (self.main_compound_enthalpy - self.main_compound_number_first_specie*self.mu_values[self.first_species] - self.main_compound_number_second_specie*self.mu_values[self.second_species]) / self.main_compound_number_third_specie
		self.mu3_display.setText(str(self.mu_values[self.third_species]))
		
		self.DefectsDiagram.mu_values[self.first_species] = self.mu_values[self.first_species]
		self.DefectsDiagram.mu_values[self.second_species] = self.mu_values[self.second_species]
		self.DefectsDiagram.mu_values[self.third_species] = self.mu_values[self.third_species]
		
		if (self.PhaseDiagram.main_compound_plot != None) and (self.PhaseDiagram.competing_compound_plots != {}):
			self.pressed_point_plot.set_data([self.mu_values[self.first_species]], [self.mu_values[self.second_species]])
			self.PhaseDiagram.Plot_Ternary_PhaseDiagram()
		
		if self.DefectsDiagram.defect_plots != {}:
			self.DefectsDiagram.Plot_Ternary_DefectsDiagram('HSE06')
	
	
	
	###############################################################################################
	##################################### Phase Diagram ###########################################
	###############################################################################################
	
	def Generate_PhaseDiagram_Plot_Function(self, event):

		# This function specifies what happens when the user clicks the "Generate Plot" button. The 
		#	only thing it needs to check is whether the elements chosen as the first, second, and third 
		#	species are unique.
		
		if len(self.species) > len(set(self.species)):	# Every time someone clicks the species' buttons, the self.species list gets updated.
														#	The "set" function checks the self.species list and omits any that are repeated.
														#	If any are repeated, then the chosen species are not unique.
			self.settings_message.setText("Pick UNIQUE elements!")	# Maybe change the phrasing ;)
		
		else:
			
			self.settings_message.setText("")
			
			# Update the number of each species in the main compound based on the above change
			self.main_compound_number_first_specie = self.compounds_info[self.main_compound][self.first_species]
			self.main_compound_number_second_specie = self.compounds_info[self.main_compound][self.second_species]
			self.main_compound_number_third_specie = self.compounds_info[self.main_compound][self.third_species]
			
			# Reset the mu value displays
			self.mu_values[self.first_species]  = 0.0
			self.mu_values[self.second_species] = 0.0
			self.mu_values[self.third_species]  = (self.main_compound_enthalpy - self.main_compound_number_first_specie*self.mu_values[self.first_species] - self.main_compound_number_second_specie*self.mu_values[self.second_species] ) / self.main_compound_number_third_specie	# Every time the mu4 value changes (and the mu1 and mu2 values are held constant), the mu3 value changes
			self.mu1_display_label.setText(u"\u0394\u03BC<sub>"+self.first_species+"</sub>")
			self.mu2_display_label.setText(u"\u0394\u03BC<sub>"+self.second_species+"</sub>")
			self.mu3_display_label.setText(u"\u0394\u03BC<sub>"+self.third_species+"</sub>")
			self.mu1_display.setText(str(self.mu_values[self.first_species]))
			self.mu2_display.setText(str(self.mu_values[self.second_species]))
			self.mu3_display.setText(str(self.mu_values[self.third_species]))
			
			# Define species
			self.PhaseDiagram.species_list = self.species_list
			self.PhaseDiagram.first_species = self.first_species
			self.PhaseDiagram.second_species = self.second_species
			self.PhaseDiagram.third_species = self.third_species
			self.PhaseDiagram.main_compound_number_first_specie = self.main_compound_number_first_specie
			self.PhaseDiagram.main_compound_number_second_specie = self.main_compound_number_second_specie
			self.PhaseDiagram.main_compound_number_third_specie = self.main_compound_number_third_specie
			
			# Set up the new plot
			self.PhaseDiagram.ternary_phase_diagram_plot_drawing.remove()
			self.PhaseDiagram.ternary_phase_diagram_plot_drawing = self.PhaseDiagram.ternary_phase_diagram_plot_figure.add_subplot(111)
			
			# Reset the plots of the main and competing compounds
			self.PhaseDiagram.main_compound_plot = None
			self.PhaseDiagram.competing_compound_plots = {}
			
			# Reset the axes of the plot
			self.PhaseDiagram.Activate_PhaseDiagram_Plot_Axes()
			self.PhaseDiagram.Update_PhaseDiagram_Plot_Axes()
			
			# Plot the phase stability diagram of the ternary compound using the new settings
			self.PhaseDiagram.Plot_Ternary_PhaseDiagram()
			
			# Reset clicked point
			self.pressed_point_plot, = self.PhaseDiagram.ternary_phase_diagram_plot_drawing.plot([], [], color="red", marker="o")
			
			
			# Reset defects diagram
			self.DefectsDiagram.ternary_defects_diagram_plot_drawing.remove()
			self.DefectsDiagram.ternary_defects_diagram_plot_drawing = self.DefectsDiagram.ternary_defects_diagram_plot_figure.add_subplot(111)
			self.DefectsDiagram.Activate_DefectsDiagram_Plot_Axes()
			self.DefectsDiagram.ternary_defects_diagram_plot_canvas.draw()
	
	
	
	###############################################################################################
	################################# Generate Defects Diagram ####################################
	###############################################################################################
	
	def Generate_DefectsDiagram_Plot_Function(self, event):
		
		self.DefectsDiagram.mu_values[self.first_species]	= self.mu_values[self.first_species]
		self.DefectsDiagram.mu_values[self.second_species]	= self.mu_values[self.second_species]
		self.DefectsDiagram.mu_values[self.third_species]	= self.mu_values[self.third_species]
		self.DefectsDiagram.first_species	= self.first_species
		self.DefectsDiagram.second_species	= self.second_species
		self.DefectsDiagram.third_species	= self.third_species
		
		self.DefectsDiagram.ternary_defects_diagram_plot_drawing.remove()
		self.DefectsDiagram.ternary_defects_diagram_plot_drawing = self.DefectsDiagram.ternary_defects_diagram_plot_figure.add_subplot(111)
		
		self.DefectsDiagram.defect_plots = {}
		self.DefectsDiagram.Activate_DefectsDiagram_Plot_Axes()
		
		if len(self.defects_diagram_Yaxis_box.text()) == 0:
			self.DefectsDiagram.Plot_Ternary_DefectsDiagram('HSE06')
			self.DefectsDiagram.ternary_defects_diagram_plot_drawing.set_xlim(self.DefectsDiagram.EVBM, self.DefectsDiagram.ECBM)
			self.DefectsDiagram.ternary_defects_diagram_plot_canvas.draw()
		else:
			self.DefectsDiagram.ymin = float(self.ternary_defects_diagram_Yaxis_box.text().split(',')[0])
			self.DefectsDiagram.ymax = float(self.ternary_defects_diagram_Yaxis_box.text().split(',')[1])
			self.DefectsDiagram.Plot_Ternary_DefectsDiagram('HSE06')
			self.DefectsDiagram.ternary_defects_diagram_plot_drawing.set_xlim(self.DefectsDiagram.EVBM, self.DefectsDiagram.ECBM)
			self.DefectsDiagram.ternary_defects_diagram_plot_drawing.set_ylim(self.DefectsDiagram.ymin, self.DefectsDiagram.ymax)
			self.DefectsDiagram.ternary_defects_diagram_plot_canvas.draw()
	
	
	
	###############################################################################################
	################################ Clicking on Phase Diagram ####################################
	###############################################################################################
	
	def Pressed_Point(self, event):
		
		# This function reads in the coordinates of the point in the phase stability diagram plot that the user clicks.

		point_x = event.xdata	# x-coordinate
		point_y = event.ydata	# y-coordinate
		
		if (not isinstance(point_x, float)) or (not isinstance(point_y, float)):	# Check that the user presses somewhere on the plot (and not anywhere else)
			pass
		
		else:
			
			# Update the clicked mu values
			self.mu_values[self.first_species] = point_x
			self.mu_values[self.second_species] = point_y
			self.mu_values[self.third_species] = (self.main_compound_enthalpy - self.main_compound_number_first_specie*self.mu_values[self.first_species] - self.main_compound_number_second_specie*self.mu_values[self.second_species]) / self.main_compound_number_third_specie
			
			# Update the display ports for the mu values
			self.mu1_display.setText(str(self.mu_values[self.first_species]))
			self.mu2_display.setText(str(self.mu_values[self.second_species]))
			self.mu3_display.setText(str(self.mu_values[self.third_species]))
			
			# Update the red dot where the user clicked on the phase diagram
			self.pressed_point_plot.set_data([self.mu_values[self.first_species]], [self.mu_values[self.second_species]])
			self.PhaseDiagram.ternary_phase_diagram_plot_figure.canvas.draw_idle()
			
			if self.DefectsDiagram.defect_plots != {}:
				self.DefectsDiagram.mu_values[self.first_species] = self.mu_values[self.first_species]
				self.DefectsDiagram.mu_values[self.second_species] = self.mu_values[self.second_species]
				self.DefectsDiagram.mu_values[self.third_species] = self.mu_values[self.third_species]
				self.DefectsDiagram.Plot_Ternary_DefectsDiagram('HSE06')
				self.DefectsDiagram.ternary_defects_diagram_plot_figure.canvas.draw_idle()
	
	
	
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
				if figure_type == "Phase Diagram":
					self.PhaseDiagram.ternary_phase_diagram_plot_figure.savefig(filename, bbox_inches='tight')
				elif figure_type == "Defects Diagram":
					self.DefectsDiagram.ternary_defects_diagram_plot_figure.savefig(filename, bbox_inches='tight')
			else:
				if figure_type == "Phase Diagram":
					self.PhaseDiagram.ternary_phase_diagram_plot_figure.savefig(filename+"."+extension, bbox_inches='tight')
				elif figure_type == "Defects Diagram":
					self.DefectsDiagram.ternary_defects_diagram_plot_figure.savefig(filename+"."+extension, bbox_inches='tight')
	
	
	###############################################################################################
	###################################### Close Window ###########################################
	###############################################################################################
	
	def closeEvent(self, event):
		msgBox = QMessageBox()
		msgBox.setText("Are you sure you want to quit VTAnDeM?")
		msgBox.setStandardButtons(QMessageBox.Yes | QMessageBox.No);
		msgBox.setDefaultButton(QMessageBox.Yes);
		if msgBox.exec_() == QMessageBox.Yes:
			event.accept()
		else:
			event.ignore()























