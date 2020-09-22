
__author__ = 'Michael_Lidia_Jiaxing_Benita_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'


###############################################################################################################################
###############################################################################################################################
################################################### Import Python Libraries ###################################################
###############################################################################################################################
###############################################################################################################################

import os
import sys
import re
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

from vtandem.dft.import_dft import *
from vtandem.dft.obtain_dft import *


###############################################################################################################################
###############################################################################################################################
################################################### Import VTAnDeM Scripts ####################################################
###############################################################################################################################
###############################################################################################################################

# Ternary scripts
from vtandem.visualization.ternary.tab_ternary_phasediagram_defectsdiagram_carrierconcentration import Tab_Ternary_PhaseDiagram_DefectsDiagram_CarrierConcentration
from vtandem.visualization.ternary.tab_ternary_phasediagram3d import Tab_Ternary_PhaseDiagram3D
from vtandem.visualization.ternary.tab_ternary_phasediagram_composition import Tab_Ternary_Compositional_PhaseDiagram
from vtandem.visualization.ternary.tab_ternary_dopants import Tab_Ternary_Dopants

# Quaternary scripts
from vtandem.visualization.quaternary.tab_quaternary_phasediagram_defectsdiagram_carrierconcentration import Tab_PhaseDiagram_DefectsDiagram_CarrierConcentration
from vtandem.visualization.quaternary.tab_quaternary_phasediagram3d import Tab_PhaseDiagram3D
from vtandem.visualization.quaternary.tab_quaternary_phasediagram_composition import Tab_Quaternary_Compositional_PhaseDiagram3D

# Binary scripts
from vtandem.visualization.binary.tab_binary_defectsdiagram_carrierconcentration import Tab_Binary_DefectsDiagram_CarrierConcentration

script_path = os.path.dirname(__file__)
vtandem_source_path = "/".join(script_path.split("/")[:-1])



###############################################################################################################################
###############################################################################################################################
################################################### Initial VTAnDeM Window ####################################################
###############################################################################################################################
###############################################################################################################################

class Welcome_VTAnDeM_Window(QMainWindow):
	
	def __init__(self):
		
		QMainWindow.__init__(self)
		
		QApplication.setStyle(QStyleFactory.create("Cleanlooks"))
		
		# Set icon as VTAnDeM logo
		self.setWindowIcon(QIcon(vtandem_source_path+"/logo/LogoSmall.png"))
		
		# Set window title
		initial_window_title = "VTAnDeM: Visualization Toolkit for Analyzing Defects in Materials"
		self.setWindowTitle(initial_window_title)
		
		# Set up layout of initial welcome dialog
		self.welcome_dialog_widget = QWidget()
		self.welcome_dialog_widget_layout = QVBoxLayout(self.welcome_dialog_widget)
		self.setCentralWidget(self.welcome_dialog_widget)
		
		# Show VTAnDeM logo
		self.vtandem_logo = QLabel()
		self.vtandem_pixmap = QPixmap(vtandem_source_path+"/logo/LogoLong.png")
		self.vtandem_pixmap_scaled = self.vtandem_pixmap.scaled(512, 512, Qt.KeepAspectRatio)
		self.vtandem_logo.setPixmap( self.vtandem_pixmap_scaled )
		self.welcome_dialog_widget_layout.addWidget(self.vtandem_logo)
		
		# VTAnDeM functionality buttons widget
		self.vtandem_functionality_widget = QWidget()
		self.vtandem_functionality_widget_layout = QHBoxLayout(self.vtandem_functionality_widget)
		
		# Import data button widget
		self.import_data_widget = QWidget()
		self.import_data_widget_layout = QVBoxLayout(self.import_data_widget)
		
		# Import compounds data button
		self.import_compounds_data_button = QPushButton("Import Compound Data")
		self.import_compounds_data_button.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
		self.import_compounds_data_button.clicked[bool].connect(self.Import_Compounds_Data_Function)
		self.import_data_widget_layout.addWidget(self.import_compounds_data_button)
		
		# Import defects data button
		self.import_defects_data_button = QPushButton("Import Defects Data")
		self.import_defects_data_button.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
		self.import_defects_data_button.clicked[bool].connect(self.Import_Defects_Data_Function)
		self.import_data_widget_layout.addWidget(self.import_defects_data_button)
		
		# Import DOS data button
		self.import_dos_data_button = QPushButton("Import DOS Data")
		self.import_dos_data_button.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
		self.import_dos_data_button.clicked[bool].connect(self.Import_DOS_Data_Function)
		self.import_data_widget_layout.addWidget(self.import_dos_data_button)
		
		# Add column of importing data buttons to window
		self.vtandem_functionality_widget_layout.addWidget(self.import_data_widget)
		
		# Visualization button widget
		self.visualize_button = QPushButton("Visualization Toolkit")
		self.visualize_button.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
		self.visualize_button.clicked[bool].connect(self.Visualize_Function)
		self.vtandem_functionality_widget_layout.addWidget(self.visualize_button)
		
		# Add visualization button to window
		self.welcome_dialog_widget_layout.addWidget(self.vtandem_functionality_widget)
		
		self.show()
	
	
	# 
	def Import_Compounds_Data_Function(self):
		
		self.import_compounds_window = Import_Data_Window("Compounds")
		self.import_compounds_window.setWindowIcon(QIcon(vtandem_source_path+"/logo/LogoSmall.png"))
		self.import_compounds_window.show()
	
	
	def Import_Defects_Data_Function(self):
		
		self.import_defects_window = Import_Data_Window("Defects")
		self.import_defects_window.setWindowIcon(QIcon(vtandem_source_path+"/logo/LogoSmall.png"))
		self.import_defects_window.show()
	
	
	def Import_DOS_Data_Function(self):
		
		self.import_dos_window = Import_Data_Window("DOS")
		self.import_dos_window.setWindowIcon(QIcon(vtandem_source_path+"/logo/LogoSmall.png"))
		self.import_dos_window.show()
	
	
	def Visualize_Function(self):
		
		self.hide()
		self.select_material_window = Material_Selection_Window()
		self.select_material_window.setWindowIcon(QIcon(vtandem_source_path+"/logo/LogoSmall.png"))
		self.select_material_window.show()




###############################################################################################################################
###############################################################################################################################
################################################## VTAnDeM Import Functions ###################################################
###############################################################################################################################
###############################################################################################################################


class Import_Data_Window(QMainWindow):
	
	def __init__(self, import_type):
		
		QMainWindow.__init__(self)
		QApplication.setStyle(QStyleFactory.create("Cleanlooks"))
		
		self.import_window = QWidget()
		self.import_window_layout = QVBoxLayout(self.import_window)
		self.setCentralWidget(self.import_window)
		
		if import_type == "Compounds":
			self.Import_Compounds_Window()
		
		elif import_type == "Defects":
			self.Import_Defects_Window()
		
		elif import_type == "DOS":
			self.Import_DOS_Window()
	
	
	def Import_Compounds_Window(self):
		
		# Compound type
		self.compound_type_widget = QWidget()
		self.compound_type_widget_layout = QHBoxLayout(self.compound_type_widget)
		self.compound_type_element = QRadioButton("Element")
		self.compound_type_element.setChecked(True)
		self.compound_type = "Element"
		self.compound_type_element.toggled.connect(lambda: self.Update_Compound_Type_Function("Element"))
		self.compound_type_widget_layout.addWidget(self.compound_type_element)
		self.compound_type_compound = QRadioButton("Compound")
		self.compound_type_compound.toggled.connect(lambda: self.Update_Compound_Type_Function("Compound"))
		self.compound_type_widget_layout.addWidget(self.compound_type_compound)
		self.import_window_layout.addWidget(self.compound_type_widget)
		
		# Compound name
		self.compound_name_widget = QWidget()
		self.compound_name_widget_layout = QHBoxLayout(self.compound_name_widget)
		self.compound_name_label = QLabel("Compound Name: ")
		self.compound_name_widget_layout.addWidget(self.compound_name_label)
		self.compound_name_prompt = QLineEdit()
		self.compound_name_widget_layout.addWidget(self.compound_name_prompt)
		self.question_button = QPushButton()
		self.question_button.setIcon(QIcon(vtandem_source_path+"/icon/QuestionIcon.png"))
		self.question_button.clicked[bool].connect(self.Import_Compounds_Help_Function)
		self.compound_name_widget_layout.addWidget(self.question_button)
		self.import_window_layout.addWidget(self.compound_name_widget)
		
		# Data directory
		self.data_directory_name_widget = QWidget()
		self.data_directory_name_widget_layout = QHBoxLayout(self.data_directory_name_widget)
		self.data_directory_name_label = QLabel("Data Directory: ")
		self.data_directory_name_widget_layout.addWidget(self.data_directory_name_label)
		self.data_directory_name_prompt = QLineEdit()
		self.data_directory_name_widget_layout.addWidget(self.data_directory_name_prompt)
		self.browser_button = QPushButton()
		self.browser_button.setIcon(QIcon(vtandem_source_path+"/icon/FolderBrowserIcon.png"))
		self.browser_button.clicked[bool].connect(self.Data_DirectoryBrowser_Function)
		self.data_directory_name_widget_layout.addWidget(self.browser_button)
		self.import_window_layout.addWidget(self.data_directory_name_widget)
		
		# Submit import request
		self.submit_data_button = QPushButton("Submit")
		self.submit_data_button.clicked[bool].connect(self.Import_Compounds_Function)
		self.import_window_layout.addWidget(self.submit_data_button)
		
		# Set the title of the window
		compounds_import_window_title = "Import Compound Data"
		self.setWindowTitle(compounds_import_window_title)
	
	
	def Update_Compound_Type_Function(self, compound_type):
		if compound_type == "Element":
			self.compound_type = "Element"
		elif compound_type == "Compound":
			self.compound_type = "Compound"
	
	def Import_Compounds_Help_Function(self):
		dialog_instructions = 	"""
Adding your VASP total energy data to the VTAnDeM database: \n \
	- Select material type (element/compound). \n \
	- Enter the name of the compound in 'Compound Name' (e.g. Cu2HgGeTe4, case-sensitive). \n \
	- Browse for the folder where the VASP data is located. \n \
	- Click on the FOLDER and hit the 'Choose' button. \n\n \
Necessary data structure: \n \
	The selected folder should contain the POSCAR and OUTCAR at the very least. \n \
	For ternary and quaternary materials, the vasprun.xml file is also needed.
								"""
		example_file_structure = 	"""
Example data structure: \n\n \
    /path/to/compounds_data 	<------ Select as Data Directory \n \
        |---- OUTCAR \n \
        |---- POSCAR
									"""
		self.message_window = QMainWindow()
		self.message_window.setWindowTitle("Help")
		self.message_window.setWindowIcon(QIcon(vtandem_source_path+"/logo/LogoSmall.png"))
		self.message_widget = QWidget()
		self.message_widget_layout = QVBoxLayout(self.message_widget)
		self.dialog_instructions = QLabel(dialog_instructions)
		self.message_widget_layout.addWidget(self.dialog_instructions)
		self.example_file_structure = QLabel(example_file_structure)
		self.message_widget_layout.addWidget(self.example_file_structure)
		self.message_window.setCentralWidget(self.message_widget)
		self.message_window.show()
	
	
	def Import_Compounds_Function(self):
		
		# Check to see if inputs are given
		self.Check_Inputs()
		
		# Extract relevant phase stability information of compound from directory
		self.compounds_info = Compounds_Import()
		compound_name = self.compound_name_prompt.text()
		data_directory_name = self.data_directory_name_prompt.text()
		if self.compound_type == "Element":
			self.compounds_info.Add_Element(compound_name, data_directory_name)
		elif self.compound_type == "Compound":
			self.compounds_info.Add_Compound(compound_name, data_directory_name)
		self.compounds_info.Update_Compounds_Database()
		
		self.close()
	
	
	
	
	
	
	def Import_Defects_Window(self):
		
		# Compound name
		self.compound_name_widget = QWidget()
		self.compound_name_widget_layout = QHBoxLayout(self.compound_name_widget)
		self.compound_name_label = QLabel("Compound Name: ")
		self.compound_name_widget_layout.addWidget(self.compound_name_label)
		self.compound_name_prompt = QLineEdit()
		self.compound_name_widget_layout.addWidget(self.compound_name_prompt)
		self.question_button = QPushButton()
		self.question_button.setIcon(QIcon(vtandem_source_path+"/icon/QuestionIcon.png"))
		self.question_button.clicked[bool].connect(self.Import_Defects_Help_Function)
		self.compound_name_widget_layout.addWidget(self.question_button)
		self.import_window_layout.addWidget(self.compound_name_widget)
		
		# Data directory
		self.data_directory_name_widget = QWidget()
		self.data_directory_name_widget_layout = QHBoxLayout(self.data_directory_name_widget)
		self.data_directory_name_label = QLabel("Data Directory: ")
		self.data_directory_name_widget_layout.addWidget(self.data_directory_name_label)
		self.data_directory_name_prompt = QLineEdit()
		self.data_directory_name_widget_layout.addWidget(self.data_directory_name_prompt)
		self.browser_button = QPushButton()
		self.browser_button.setIcon(QIcon(vtandem_source_path+"/icon/FolderBrowserIcon.png"))
		self.browser_button.clicked[bool].connect(self.Data_DirectoryBrowser_Function)
		self.data_directory_name_widget_layout.addWidget(self.browser_button)
		self.import_window_layout.addWidget(self.data_directory_name_widget)
		
		# Supercell size
		self.supercell_size_widget = QWidget()
		self.supercell_size_widget_layout = QVBoxLayout(self.supercell_size_widget)
		self.supercell_size_x_widget = QWidget()
		self.supercell_size_x_widget_layout = QHBoxLayout(self.supercell_size_x_widget)
		self.supercell_size_x_label = QLabel("Supercell Size, x: ")
		self.supercell_size_x_widget_layout.addWidget(self.supercell_size_x_label)
		self.supercell_size_x_prompt = QLineEdit()
		self.supercell_size_x_widget_layout.addWidget(self.supercell_size_x_prompt)
		self.supercell_size_widget_layout.addWidget(self.supercell_size_x_widget)
		self.supercell_size_y_widget = QWidget()
		self.supercell_size_y_widget_layout = QHBoxLayout(self.supercell_size_y_widget)
		self.supercell_size_y_label = QLabel("Supercell Size, y: ")
		self.supercell_size_y_widget_layout.addWidget(self.supercell_size_y_label)
		self.supercell_size_y_prompt = QLineEdit()
		self.supercell_size_y_widget_layout.addWidget(self.supercell_size_y_prompt)
		self.supercell_size_widget_layout.addWidget(self.supercell_size_y_widget)
		self.supercell_size_z_widget = QWidget()
		self.supercell_size_z_widget_layout = QHBoxLayout(self.supercell_size_z_widget)
		self.supercell_size_z_label = QLabel("Supercell Size, z: ")
		self.supercell_size_z_widget_layout.addWidget(self.supercell_size_z_label)
		self.supercell_size_z_prompt = QLineEdit()
		self.supercell_size_z_widget_layout.addWidget(self.supercell_size_z_prompt)
		self.supercell_size_widget_layout.addWidget(self.supercell_size_z_widget)
		self.import_window_layout.addWidget(self.supercell_size_widget)
		
		# Submit import request
		self.submit_data_button = QPushButton("Submit")
		self.submit_data_button.clicked[bool].connect(self.Import_Defects_Function)
		self.import_window_layout.addWidget(self.submit_data_button)
		
		# Set the title of the window
		defects_import_window_title = "Import Defects Data"
		self.setWindowTitle(defects_import_window_title)
	
	
	def Import_Defects_Help_Function(self):
		dialog_instructions = 	"""
Adding your VASP defects data to the VTAnDeM database: \n \
	- Enter the name of the compound in 'Compound Name' (e.g. Cu2HgGeTe4, case-sensitive). \n \
	- Browse for the folder where the defects data is located. \n \
	- Click on the FOLDER and hit the 'Choose' button. \n\n \
Necessary defects data structure: \n \
	- The folder should contain a folder for each defect, split by an underscore (e.g. Cu_Hg for Cu antisite on Hg site). \n \
	- Within each defect, a folder for each charge state should exist, named 'q#' where # is the charge state. \n \
	- In each charge state, the OUTCAR or OSZICAR file must be found.
								"""
		example_file_structure = 	"""
Example data structure: \n\n \
    /path/to/defects_data 	<------ Select as Data Directory \n \
        |---- Cu_Hg \n \
                |---- q0 \n \
                        |---- OUTCAR \n \
                |---- q+1 \n \
                        |---- OUTCAR \n \
                |---- q-1 \n \
                        |---- OUTCAR \n \
        |---- V_Cu \n \
                |---- q0 \n \
                        |---- OUTCAR \n \
                |---- q-1 \n \
                        |---- OUTCAR \n \
        |---- V_Te \n \
                |---- q0 \n \
                        |---- OUTCAR \n \
                |---- q+1 \n \
                        |---- OUTCAR \n \
                |---- q+2 \n \
                        |---- OUTCAR
									"""
		self.message_window = QMainWindow()
		self.message_window.setWindowTitle("Help")
		self.message_window.setWindowIcon(QIcon(vtandem_source_path+"/logo/LogoSmall.png"))
		self.message_widget = QWidget()
		self.message_widget_layout = QVBoxLayout(self.message_widget)
		self.dialog_instructions = QLabel(dialog_instructions)
		self.message_widget_layout.addWidget(self.dialog_instructions)
		self.example_file_structure = QLabel(example_file_structure)
		self.message_widget_layout.addWidget(self.example_file_structure)
		self.message_window.setCentralWidget(self.message_widget)
		self.message_window.show()
	
	
	def Import_Defects_Function(self):
		
		# Check to see if inputs are given
		self.Check_Inputs()
		
		# Check to see that supercell sizes are given
		try: 
			supercell_size = float(self.supercell_size_x_prompt.text()) * float(self.supercell_size_y_prompt.text()) * float(self.supercell_size_z_prompt.text())
		except:
			print("Supercell sizes are not properly defined.")
			self.close()
		
		# Extract relevant compounds information from directory
		self.defects_data = Defects_Import()
		compound_name = self.compound_name_prompt.text()
		data_directory_name = self.data_directory_name_prompt.text()
		self.defects_data.Add_Defects(compound_name, data_directory_name, supercell_size)
		self.defects_data.Update_Defects_Database()
		self.close()
	
	
	
	
	
	def Import_DOS_Window(self):
		
		# Compound name
		self.compound_name_widget = QWidget()
		self.compound_name_widget_layout = QHBoxLayout(self.compound_name_widget)
		self.compound_name_label = QLabel("Compound Name: ")
		self.compound_name_widget_layout.addWidget(self.compound_name_label)
		self.compound_name_prompt = QLineEdit()
		self.compound_name_widget_layout.addWidget(self.compound_name_prompt)
		self.question_button = QPushButton()
		self.question_button.setIcon(QIcon(vtandem_source_path+"/icon/QuestionIcon.png"))
		self.question_button.clicked[bool].connect(self.Import_DOS_Help_Function)
		self.compound_name_widget_layout.addWidget(self.question_button)
		self.import_window_layout.addWidget(self.compound_name_widget)
		
		# Data file (DOSCAR)
		self.data_filename_widget = QWidget()
		self.data_filename_widget_layout = QHBoxLayout(self.data_filename_widget)
		self.data_filename_label = QLabel("DOSCAR File: ")
		self.data_filename_widget_layout.addWidget(self.data_filename_label)
		self.data_filename_prompt = QLineEdit()
		self.data_filename_widget_layout.addWidget(self.data_filename_prompt)
		self.browser_button = QPushButton()
		self.browser_button.setIcon(QIcon(vtandem_source_path+"/icon/FolderBrowserIcon.png"))
		self.browser_button.clicked[bool].connect(self.Data_FileBrowser_Function)
		self.data_filename_widget_layout.addWidget(self.browser_button)
		self.import_window_layout.addWidget(self.data_filename_widget)
		
		# Submit import request
		self.submit_data_button = QPushButton("Submit")
		self.submit_data_button.clicked[bool].connect(self.Import_DOS_Function)
		self.import_window_layout.addWidget(self.submit_data_button)
		
		# Set the title of the window
		dos_import_window_title = "Import DOS Data"
		self.setWindowTitle(dos_import_window_title)
	
	
	
	def Import_DOS_Help_Function(self):
		dialog_instructions = 	"""
Adding your VASP Density of States data to the VTAnDeM database: \n \
	- Enter the name of the compound in 'Compound Name' (e.g. Cu2HgGeTe4, case-sensitive). \n \
	- Browse for the DOSCAR file of the compound. \n \
	- Click on the FILE and hit the 'Choose' button.
								"""
		self.message_window = QMainWindow()
		self.message_window.setWindowTitle("Help")
		self.message_window.setWindowIcon(QIcon(vtandem_source_path+"/logo/LogoSmall.png"))
		self.message_widget = QWidget()
		self.message_widget_layout = QVBoxLayout(self.message_widget)
		self.dialog_instructions = QLabel(dialog_instructions)
		self.message_widget_layout.addWidget(self.dialog_instructions)
		self.message_window.setCentralWidget(self.message_widget)
		self.message_window.show()
	
	
	
	def Import_DOS_Function(self):
		
		# Check to see if inputs are given
		self.Check_Inputs()
		
		# Extract DOS information of compound from DOSCAR
		self.dos_data = DOS_Import()
		compound_name = self.compound_name_prompt.text()
		data_filename = self.data_filename_prompt.text()
		self.dos_data.Add_DOS(compound_name, data_filename)
		self.dos_data.Update_DOS_Database()
		
		self.close()
	
	
	
	
	def Check_Inputs(self):
		
		# Check to see that compound name is given
		if self.compound_name_prompt.text() == "":
			print("No compound name given. Exiting...")
			self.close()
		
		# Check to see that name of data is given
		try:
			if self.data_directory_name_prompt.text() == "":
				print("No name for the data directory is given. Exiting...")
				self.close()
		except:
			if self.data_filename_prompt.text() == "":
				print("No name for file is given. Exiting...")
				self.close()
	
	
	def Data_DirectoryBrowser_Function(self):
		filename = QFileDialog.getExistingDirectory(self, "Browse Data Directory")
		#filename = QFileDialog.getOpenFileNames(self, "Browse Data Directory")
		self.data_directory_name_prompt.setText(filename)
	
	
	def Data_FileBrowser_Function(self):
		filename = QFileDialog.getOpenFileNames(self, "Browse Data File")
		self.data_filename_prompt.setText(filename[0][0])








###############################################################################################################################
###############################################################################################################################
##################################################### Material Selection ######################################################
###############################################################################################################################
###############################################################################################################################


class Material_Selection_Window(QMainWindow):
	
	def __init__(self):
		
		QMainWindow.__init__(self)
		
		QApplication.setStyle(QStyleFactory.create("Cleanlooks"))
		
		# Set icon for window as VTAnDeM logo
		self.setWindowIcon(QIcon(vtandem_source_path+"/logo/LogoSmall.png"))
		
		# Set window title
		material_selection_window_title = "Material Selection Hub"
		self.setWindowTitle(material_selection_window_title)
		
		# Set up the layout of the initial dialog
		self.initial_dialog_widgets = QWidget()
		self.setCentralWidget(self.initial_dialog_widgets)
		self.initial_dialog_widgets_layout = QVBoxLayout(self.initial_dialog_widgets)
		
		# Show VTAnDeM logo
		self.vtandem_logo = QLabel()
		self.vtandem_pixmap = QPixmap(vtandem_source_path+"/logo/LogoLong.png")
		self.vtandem_pixmap_scaled = self.vtandem_pixmap.scaled(512, 512, Qt.KeepAspectRatio)
		self.vtandem_logo.setPixmap( self.vtandem_pixmap_scaled )
		self.initial_dialog_widgets_layout.addWidget(self.vtandem_logo)
		
		# Options widget
		self.vtandem_options_widget = QWidget()
		self.vtandem_options_widget_layout = QHBoxLayout(self.vtandem_options_widget)
		
		# Material selection tree widget
		self.compounds_tree = QTreeWidget()
		self.compounds_tree.setHeaderLabels(["Please select material:"])
		
		# Create tree structure for set of available compounds
		self.binary_compounds_set = QTreeWidgetItem(["Binary"])
		self.ternary_compounds_set = QTreeWidgetItem(["Ternary"])
		self.quaternary_compounds_set = QTreeWidgetItem(["Quaternary"])
		
		# Open compounds data
		with open("Compounds_Tracker.json") as CompoundsTracker:
			self.compounds_info = json.load(CompoundsTracker)
		
		# Open defects data
		with open("Defects_Tracker.json") as DefectsTracker:
			self.defects_data = json.load(DefectsTracker)
		
		# Open DOS data
		with open("DOS_Tracker.json") as DOSTracker:
			self.dos_data = json.load(DOSTracker)
		
		# Add compounds to tree
		for compound in self.compounds_info["Compounds"].keys():
			if self.compounds_info["Compounds"][compound]["number_species"] == 2:
				self.binary_compounds_set.addChild(QTreeWidgetItem([compound]))
			elif self.compounds_info["Compounds"][compound]["number_species"] == 3:
				self.ternary_compounds_set.addChild(QTreeWidgetItem([compound]))
			elif self.compounds_info["Compounds"][compound]["number_species"] == 4:
				self.quaternary_compounds_set.addChild(QTreeWidgetItem([compound]))
		
		self.compounds_tree.itemClicked.connect(self.Select_Visualization)
		
		# Add above features to Tree Widget
		self.compounds_tree.addTopLevelItem(self.binary_compounds_set)
		self.compounds_tree.addTopLevelItem(self.ternary_compounds_set)
		self.compounds_tree.addTopLevelItem(self.quaternary_compounds_set)
		self.vtandem_options_widget_layout.addWidget(self.compounds_tree)
		
		# Visualization options widget
		self.vtandem_visualization_options_widget = QWidget()
		self.vtandem_visualization_options_widget_layout = QVBoxLayout(self.vtandem_visualization_options_widget)
		
		# Phase stability visualization option
		self.phase_stability_checkbox = QCheckBox("Phase Stability")
		self.phase_stability_checkbox.setEnabled(False)
		self.phase_stability_checkbox.setStyleSheet("color: gray")
		self.vtandem_visualization_options_widget_layout.addWidget(self.phase_stability_checkbox)
		
		# Defects diagram visualization option
		self.defects_diagram_checkbox = QCheckBox("Defects Diagram")
		self.defects_diagram_checkbox.setEnabled(False)
		self.defects_diagram_checkbox.setStyleSheet("color: gray")
		self.defects_diagram_checkbox.clicked.connect(self.Allow_Visualization)
		self.vtandem_visualization_options_widget_layout.addWidget(self.defects_diagram_checkbox)
		
		# Carrier concentration visualization option
		self.carrier_concentration_checkbox = QCheckBox("Carrier Concentration")
		self.carrier_concentration_checkbox.setEnabled(False)
		self.carrier_concentration_checkbox.setStyleSheet("color: gray")
		self.vtandem_visualization_options_widget_layout.addWidget(self.carrier_concentration_checkbox)
		
		# Button to generate main VTAnDeM window
		self.select_material_button = QPushButton("Visualize!")
		self.select_material_button.clicked[bool].connect(self.Generate_VTAnDeM_Window)
		self.vtandem_visualization_options_widget_layout.addWidget(self.select_material_button)
		
		# Add visualization options widget to options widget
		self.vtandem_options_widget_layout.addWidget(self.vtandem_visualization_options_widget)
		
		# Add options widget to initial dialog widget
		self.initial_dialog_widgets_layout.addWidget(self.vtandem_options_widget)
	
	
	
	
	def Select_Visualization(self):
		
		selected_branch = self.compounds_tree.selectedItems()
		if selected_branch:
			selected_branch_object = selected_branch[0]
			compound_name = selected_branch_object.text(0)
			
			if compound_name not in self.compounds_info["Compounds"].keys():
				self.phase_stability_checkbox.setEnabled(False)
				self.phase_stability_checkbox.setChecked(False)
				self.phase_stability_checkbox.setStyleSheet("color: gray")
				self.defects_diagram_checkbox.setEnabled(False)
				self.defects_diagram_checkbox.setChecked(False)
				self.defects_diagram_checkbox.setStyleSheet("color: gray")
				self.carrier_concentration_checkbox.setEnabled(False)
				self.carrier_concentration_checkbox.setChecked(False)
				self.carrier_concentration_checkbox.setStyleSheet("color: gray")
				return
			
			if self.compounds_info["Compounds"][compound_name]["number_species"] == 2:
				self.phase_stability_checkbox.setChecked(False)
				self.phase_stability_checkbox.setStyleSheet("color: gray")
			else:
				self.phase_stability_checkbox.setChecked(True)
				self.phase_stability_checkbox.setStyleSheet("color: black")
			
			if compound_name in self.defects_data.keys():
				self.defects_diagram_checkbox.setEnabled(True)
				self.defects_diagram_checkbox.setChecked(True)
				self.defects_diagram_checkbox.setStyleSheet("color: black")
				if (compound_name in self.dos_data.keys()) and (self.dos_data[compound_name] != {}):
					self.carrier_concentration_checkbox.setEnabled(True)
					self.carrier_concentration_checkbox.setChecked(True)
					self.carrier_concentration_checkbox.setStyleSheet("color: black")
			else:
				self.defects_diagram_checkbox.setEnabled(False)
				self.defects_diagram_checkbox.setChecked(False)
				self.defects_diagram_checkbox.setStyleSheet("color: gray")
				self.carrier_concentration_checkbox.setEnabled(False)
				self.carrier_concentration_checkbox.setChecked(False)
				self.carrier_concentration_checkbox.setStyleSheet("color: gray")
	
	
	
	def Allow_Visualization(self):
		
		if self.defects_diagram_checkbox.isChecked():
			self.carrier_concentration_checkbox.setEnabled(True)
		else:
			self.carrier_concentration_checkbox.setChecked(False)
			self.carrier_concentration_checkbox.setEnabled(False)
	
	
	
	def Generate_VTAnDeM_Window(self):
		
		selected_branch = self.compounds_tree.selectedItems()
		
		if selected_branch:
			selected_branch_object = selected_branch[0]
			compound_name = selected_branch_object.text(0)
			
			if compound_name not in self.compounds_info["Compounds"].keys():
				return
			
			compound_type = selected_branch_object.parent().text(0)
			
			show_defects_diagram = self.defects_diagram_checkbox.isChecked()
			show_carrier_concentration = self.carrier_concentration_checkbox.isChecked()
			
			self.hide()
			
			if compound_type == "Binary":
				binary_species_list = [ re.sub(r'[0-9]+', '', specie) for specie in re.findall( "[A-Z][^A-Z]*", compound_name ) ]
				self.binary_aw = Binary_Main_VTAnDeM_Window(main_compound = compound_name, first_element = binary_species_list[0], second_element = binary_species_list[1], show_defects_diagram = show_defects_diagram, show_carrier_concentration = show_carrier_concentration)
				self.binary_aw.show()
			
			elif compound_type == "Ternary":
				ternary_species_list = [ re.sub(r'[0-9]+', '', specie) for specie in re.findall( "[A-Z][^A-Z]*", compound_name ) ]
				self.ternary_aw = Ternary_Main_VTAnDeM_Window(main_compound = compound_name, first_element = ternary_species_list[0], second_element = ternary_species_list[1], third_element = ternary_species_list[2], show_defects_diagram = show_defects_diagram, show_carrier_concentration = show_carrier_concentration)
				self.ternary_aw.show()
			
			elif compound_type == "Quaternary":
				quaternary_species_list = [ re.sub(r'[0-9]+', '', specie) for specie in re.findall( "[A-Z][^A-Z]*", compound_name ) ]
				self.quaternary_aw = Quaternary_Main_VTAnDeM_Window(main_compound = compound_name, first_element = quaternary_species_list[0], second_element = quaternary_species_list[1], third_element = quaternary_species_list[2], fourth_element = quaternary_species_list[3], show_defects_diagram = show_defects_diagram, show_carrier_concentration = show_carrier_concentration)
				self.quaternary_aw.show()










###############################################################################################################################
###############################################################################################################################
########################################### VTAnDeM Window for Quaternary Materials ###########################################
###############################################################################################################################
###############################################################################################################################


class Quaternary_Main_VTAnDeM_Window(QMainWindow):
	
	def __init__(self, parent = None, main_compound = None, first_element = None, second_element = None, third_element = None, fourth_element = None, show_defects_diagram = True, show_carrier_concentration = True):	# User specifies the main compound and its constituents
		
		# Inherit all initial variables from the QMainWindow class
		QMainWindow.__init__(self)
		self.setWindowIcon(QIcon(vtandem_source_path+"/logo/LogoSmall.png"))
		
		
		# Font description for phase stability diagram plot
		self.font = {'family': 'sans-serif',
				'color':  'black',
				'weight': 'normal',
				'size': 16 }
		
		
		# Set up the framework of the application window (including file menu, exit function, etc.)
		self.Setup_Window_Framework()
		
		
		
		
		self.elements_list = [first_element, second_element, third_element, fourth_element]					# Species list (order MAY change)
		
		
		
		# Obtain DFT data
		self.compounds_info = Obtain_Compounds_Data(self.elements_list)		# Total energies/enthalpies for phase diagram
		self.defects_data = Obtain_Defects_Data()	# Defect energies for defects diagram
		self.dos_data = Obtain_DOS_Data()
		
		
		
		self.Tab1_PhasesDefectsCarriers_Object = Tab_PhaseDiagram_DefectsDiagram_CarrierConcentration(self, main_compound = main_compound, first_element = first_element, second_element = second_element, third_element = third_element, fourth_element = fourth_element, compounds_info = self.compounds_info, defects_data = self.defects_data, dos_data = self.dos_data, show_defects_diagram = show_defects_diagram, show_carrier_concentration = show_carrier_concentration)
		self.Tab2_PhaseDiagram3D_Object = Tab_PhaseDiagram3D(self, main_compound = main_compound, first_element = first_element, second_element = second_element, third_element = third_element, fourth_element = fourth_element, compounds_info = self.compounds_info)
		self.Tab3_PhaseDiagram3D_Object = Tab_Quaternary_Compositional_PhaseDiagram3D(main_compound = main_compound, first_element = first_element, second_element = second_element, third_element = third_element, fourth_element = fourth_element, compounds_info = self.compounds_info, defects_data = self.defects_data, show_defects_diagram = show_defects_diagram)
		
		
		
		###############################################################################################
		################################## Initialize PyQt5 widgets ###################################
		###############################################################################################
		
		# Set up layout of widgets (think of "widgets" as being like objects in the app, like a button or a plot)
		self.widgets = QWidget()						# "Main" widget. This is where all other widgets will be placed.
														# 	The reason why we need a "main" widget is because PyQT only allows the programmer
														#	to declare only ONE central widget, where everything happens. Since we need to pool
														#	into this app other widgets (like buttons and the quaternary phase diagram plot),
														#	we place all of our widgets into this "main" widget.
		self.setCentralWidget(self.widgets)				# Declare self.widgets as the main (or central) widget. Subsequent widgets (like
														#	the scroll bar or buttons) will be added to this main widget.
		self.widgets_grid = QHBoxLayout(self.widgets)	# The layout of the main widget should be such that all widgets placed inside the main
														#	widget (i.e. buttons, plots, etc.) are placed horizontally. In the case of VTAnDeM,
														#	the widgets will be the phase diagram (and associated widgets), the defects diagram
														#	(with associated widgets), and the carrier concentration (again with associated widgets).
		
		
		
		
		self.plot_tabs_widget = QTabWidget()	# Tabs will be used to organize the layout of VTAnDeM.
		self.plot_tabs_widget.addTab(self.Tab1_PhasesDefectsCarriers_Object.tab1, "Phases and Defects")
		self.plot_tabs_widget.addTab(self.Tab2_PhaseDiagram3D_Object.tab2, "Phase Diagram, Chemical Potential Space")
		self.plot_tabs_widget.addTab(self.Tab3_PhaseDiagram3D_Object.tab3, "Phase Diagram, Composition Space")
		
		
		
		
		self.widgets_grid.addWidget(self.plot_tabs_widget)
		
		
		self.showFullScreen()
	
	
	
	###############################################################################################
	################################# Main Window Framework #######################################
	###############################################################################################
	
	def Setup_Window_Framework(self):
		
		# Create a menu bar
		menubar = self.menuBar()	# This menu bar automatically gets added to the existing main window
		
		# Add options (created above) to the menu bar, in order
		fileMenu = menubar.addMenu("&File")
		
		# Set up an option where the user can open a new window (not functionalized yet)
		newappAction = QAction("New", self)
		newappAction.triggered.connect(self.Open_MaterialSelectionWindow)
		fileMenu.addAction(newappAction)
		
		# Add "About" section
		aboutMenu = menubar.addMenu("&About")
		
		# Set up an option to view a brief description of VTAnDeM
		aboutAction = QAction("About", self)
		aboutAction.triggered.connect(self.VTAnDeM_Description)
		aboutMenu.addAction(aboutAction)
		
		# Set the title of the window
		mainwindow_title = "VTAnDeM: Visualization Toolkit for Analyzing Defects in Materials"
		self.setWindowTitle(mainwindow_title)
	
	
	
	def VTAnDeM_Description(self):
		
		description_window = QMainWindow(self)
		
		### Insert VTAnDeM description
		
		description_window.show()
	
	
	
	
	def Open_MaterialSelectionWindow(self):
		
		self.select_material_window = Material_Selection_Window()
		self.select_material_window.setWindowIcon(QIcon(vtandem_source_path+"/logo/LogoSmall.png"))
		self.select_material_window.show()
		
		self.close()





###############################################################################################################################
###############################################################################################################################
############################################# VTAnDeM Window for Ternary Materials ############################################
###############################################################################################################################
###############################################################################################################################

class Ternary_Main_VTAnDeM_Window(QMainWindow):
	
	def __init__(self, parent = None, main_compound = None, first_element = None, second_element = None, third_element = None, show_defects_diagram = True, show_carrier_concentration = True):
		
		# Inherit all initial variables from the QMainWindow class
		QMainWindow.__init__(self)
		self.setWindowIcon(QIcon(vtandem_source_path+"/logo/LogoSmall.png"))
		
		
		# Set up the framework of the application window (including file menu, exit function, etc.)
		self.Setup_Window_Framework()
		
		# Establish list of elements in compound
		self.elements_list = [first_element, second_element, third_element]					# Species list (order MAY change)
		
		# Obtain DFT data
		self.compounds_info = Obtain_Compounds_Data(self.elements_list)	# Total energies/enthalpies for phase diagram
		self.defects_data = Obtain_Defects_Data()		# Defect energies for defects diagram
		self.dos_data = Obtain_DOS_Data()
		
		
		self.Tab1_PhasesDefectsCarriers_Object = Tab_Ternary_PhaseDiagram_DefectsDiagram_CarrierConcentration(main_compound = main_compound, first_element = first_element, second_element = second_element, third_element = third_element, compounds_info = self.compounds_info, defects_data = self.defects_data, dos_data = self.dos_data, show_defects_diagram = show_defects_diagram, show_carrier_concentration = show_carrier_concentration)
		self.Tab2_PhaseDiagram3D_Object = Tab_Ternary_PhaseDiagram3D(self, main_compound = main_compound, first_element = first_element, second_element = second_element, third_element = third_element, compounds_info = self.compounds_info)
		self.Tab3_PhaseDiagram_Object = Tab_Ternary_Compositional_PhaseDiagram(main_compound = main_compound, first_element = first_element, second_element = second_element, third_element = third_element, compounds_info = self.compounds_info, defects_data = self.defects_data, show_defects_diagram = show_defects_diagram)
		
		if show_defects_diagram:
			self.Tab4_Ternary_Dopants = Tab_Ternary_Dopants(main_compound = main_compound, first_element = first_element, second_element = second_element, third_element = third_element, compounds_info = self.compounds_info, defects_data = self.defects_data)
		
		
		
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
		
		
		self.plot_tabs_widget = QTabWidget()
		
		
		self.plot_tabs_widget.addTab(self.Tab1_PhasesDefectsCarriers_Object.tab1, "Phases and Defects")
		self.plot_tabs_widget.addTab(self.Tab2_PhaseDiagram3D_Object.tab2, "Phase Diagram, Chemical Potential Space")
		self.plot_tabs_widget.addTab(self.Tab3_PhaseDiagram_Object.tab3, "Phase Diagram, Composition Space")
		
		if show_defects_diagram:
			self.plot_tabs_widget.addTab(self.Tab4_Ternary_Dopants.tab4, "Dopants")
		
		
		self.widgets_grid.addWidget(self.plot_tabs_widget)
		
		
		self.showFullScreen()
	
	
	
	
	
	
	###############################################################################################
	################################# Main Window Framework #######################################
	###############################################################################################
	
	def Setup_Window_Framework(self):
		
		# Create a menu bar
		menubar = self.menuBar()	# This menu bar automatically gets added to the existing main window
		
		# Add file options to the menu bar
		fileMenu = menubar.addMenu('&File')
		
		# Set up an option where the user can go back to the material selection window
		newappAction = QAction('New', self)
		newappAction.triggered.connect(self.Open_MaterialSelectionWindow)
		fileMenu.addAction(newappAction)
		
		# Add "About" section
		aboutMenu = menubar.addMenu("&About")
		
		# Set up an option to view a brief description of VTAnDeM
		aboutAction = QAction("About", self)
		aboutAction.triggered.connect(self.VTAnDeM_Description)
		aboutMenu.addAction(aboutAction)
		
		# Set the title of the window
		mainwindow_title = "VTAnDeM: Visualization Toolkit for Analyzing Defects in Materials"
		self.setWindowTitle(mainwindow_title)
	
	
	
	
	def VTAnDeM_Description(self):
		
		description_window = QMainWindow(self)
		
		description_window.show()
	
	
	
	
	
	
	def Open_MaterialSelectionWindow(self):
		
		self.select_material_window = Material_Selection_Window()
		self.select_material_window.setWindowIcon(QIcon(vtandem_source_path+"/logo/LogoSmall.png"))
		self.select_material_window.show()
		
		self.close()







###############################################################################################################################
###############################################################################################################################
############################################# VTAnDeM Window for Binary Materials ############################################
###############################################################################################################################
###############################################################################################################################

class Binary_Main_VTAnDeM_Window(QMainWindow):
	
	def __init__(self, parent = None, main_compound = None, first_element = None, second_element = None, show_defects_diagram = True, show_carrier_concentration = True):
		
		# Inherit all initial variables from the QMainWindow class
		QMainWindow.__init__(self)
		self.setWindowIcon(QIcon(vtandem_source_path+"/logo/LogoSmall.png"))
		
		
		# Set up the framework of the application window (including file menu, exit function, etc.)
		self.Setup_Window_Framework()
		
		# Establish list of elements in compound
		self.elements_list = [first_element, second_element]					# Species list (order MAY change)
		
		# Obtain DFT data
		self.compounds_info = Obtain_Compounds_Data(self.elements_list)	# Total energies/enthalpies for phase diagram
		self.defects_data = Obtain_Defects_Data()		# Defect energies for defects diagram
		self.dos_data = Obtain_DOS_Data()
		
		
		self.Tab1_PhasesDefectsCarriers_Object = Tab_Binary_DefectsDiagram_CarrierConcentration(self, main_compound = main_compound, first_element = first_element, second_element = second_element, compounds_info = self.compounds_info, defects_data = self.defects_data, dos_data = self.dos_data, show_defects_diagram = show_defects_diagram, show_carrier_concentration = show_carrier_concentration)
		
		
		
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
		
		
		self.plot_tabs_widget = QTabWidget()
		
		
		self.plot_tabs_widget.addTab(self.Tab1_PhasesDefectsCarriers_Object.tab1, "Phases and Defects")
		#self.plot_tabs_widget.addTab(self.Tab2_PhaseDiagram3D_Object.tab2, "Phase Diagram, Chemical Potential Space")
		#self.plot_tabs_widget.addTab(self.Tab3_PhaseDiagram_Object.tab3, "Phase Diagram, Composition Space")
		
		
		
		self.widgets_grid.addWidget(self.plot_tabs_widget)
		
		
		self.showFullScreen()
	
	
	
	
	
	
	###############################################################################################
	################################# Main Window Framework #######################################
	###############################################################################################
	
	def Setup_Window_Framework(self):
		
		# Create a menu bar
		menubar = self.menuBar()	# This menu bar automatically gets added to the existing main window
		
		# Add file options to the menu bar
		fileMenu = menubar.addMenu('&File')
		
		# Set up an option where the user can go back to the material selection window
		newappAction = QAction('New', self)
		newappAction.triggered.connect(self.Open_MaterialSelectionWindow)
		fileMenu.addAction(newappAction)
		
		# Add "About" section
		aboutMenu = menubar.addMenu("&About")
		
		# Set up an option to view a brief description of VTAnDeM
		aboutAction = QAction("About", self)
		aboutAction.triggered.connect(self.VTAnDeM_Description)
		aboutMenu.addAction(aboutAction)
		
		# Set the title of the window
		mainwindow_title = "VTAnDeM: Visualization Toolkit for Analyzing Defects in Materials"
		self.setWindowTitle(mainwindow_title)
	
	
	
	
	def VTAnDeM_Description(self):
		
		description_window = QMainWindow(self)
		
		description_window.show()
	
	
	
	
	
	
	def Open_MaterialSelectionWindow(self):
		
		self.select_material_window = Material_Selection_Window()
		self.select_material_window.setWindowIcon(QIcon(vtandem_source_path+"/logo/LogoSmall.png"))
		self.select_material_window.show()
		
		self.close()







def Open_Material_Selection_Window():
	app = QApplication([])
	ex = Material_Selection_Window()
	ex.show()
	sys.exit(app.exec_())


def Open_Welcome_VTAnDeM_Window():
	app = QApplication([])
	ex = Welcome_VTAnDeM_Window()
	sys.exit(app.exec_())





