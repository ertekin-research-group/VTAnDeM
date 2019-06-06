
__author__ = 'Michael_Lidia_Jiaxing_Benita_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'


import sys
import re
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *
from ternary_main_window_vtandem import Ternary_Main_VTAnDeM_Window
from quaternary_main_window_vtandem import Quaternary_Main_VTAnDeM_Window


class Initial_VTAnDeM_Dialog(QMainWindow):
	
	def __init__(self):
		
		QMainWindow.__init__(self)
		
		QApplication.setStyle(QStyleFactory.create("Cleanlooks"))
		
		# Set up the layout of the initial dialog
		self.initial_dialog_widgets = QWidget()
		self.setCentralWidget(self.initial_dialog_widgets)
		self.initial_dialog_widgets_layout = QVBoxLayout(self.initial_dialog_widgets)
		
		# Show VTAnDeM logo
		self.vtandem_logo = QLabel()
		self.vtandem_pixmap = QPixmap("LogoLong.png")
		self.vtandem_pixmap_scaled = self.vtandem_pixmap.scaled(512, 512, Qt.KeepAspectRatio)
		self.vtandem_logo.setPixmap( self.vtandem_pixmap_scaled )
		self.initial_dialog_widgets_layout.addWidget(self.vtandem_logo)
		
		
		###### "Tree" Widget
		self.compounds_tree = QTreeWidget()
		self.compounds_tree.setHeaderLabels(["Please select material:"])
		
		# Create tree structure for set of available compounds
		self.ternary_compounds_set = QTreeWidgetItem(["Ternary"])
		self.quaternary_compounds_set = QTreeWidgetItem(["Quaternary"])
		
		# Add compounds to tree
		child1 = QTreeWidgetItem(["CuGaTe2"])
		child3 = QTreeWidgetItem(["Hg2GeTe4"])
		child4 = QTreeWidgetItem(["Cu2GeTe3"])
		self.ternary_compounds_set.addChild(child1)
		self.ternary_compounds_set.addChild(child3)
		self.ternary_compounds_set.addChild(child4)
		
		child2 = QTreeWidgetItem(["Cu2HgGeTe4"])
		self.quaternary_compounds_set.addChild(child2)
		
		# Add above features to Tree Widget
		self.compounds_tree.addTopLevelItem(self.ternary_compounds_set)
		self.compounds_tree.addTopLevelItem(self.quaternary_compounds_set)
		self.initial_dialog_widgets_layout.addWidget(self.compounds_tree)
		
		
		# Button to generate main VTAnDeM window
		self.select_material_button = QPushButton("Visualize!")
		self.select_material_button.clicked[bool].connect(self.Generate_VTAnDeM_Window)
		self.initial_dialog_widgets_layout.addWidget(self.select_material_button)
		
		self.show()
	
	
	def Generate_VTAnDeM_Window(self):
		selected_branch = self.compounds_tree.selectedItems()
		if selected_branch:
			selected_branch_object = selected_branch[0]
			compound_name = selected_branch_object.text(0)
			compound_type = selected_branch_object.parent().text(0)
			
			self.hide()
			
			if compound_type == "Ternary":
				ternary_species_list = [ re.sub(r'[0-9]+', '', specie) for specie in re.findall( "[A-Z][^A-Z]*", compound_name ) ]
				self.ternary_aw = Ternary_Main_VTAnDeM_Window(main_compound = compound_name, first_species = ternary_species_list[0], second_species = ternary_species_list[1], third_species = ternary_species_list[2])
				self.ternary_aw.show()
				print(self.ternary_aw.isWindow())
				"""
				if self.ternary_aw.closeEvent():
					print("ehllo")
				"""
			
			
			elif compound_type == "Quaternary":
				quaternary_species_list = [ re.sub(r'[0-9]+', '', specie) for specie in re.findall( "[A-Z][^A-Z]*", compound_name ) ]
				self.quaternary_aw = Quaternary_Main_VTAnDeM_Window(main_compound = compound_name, first_species = quaternary_species_list[0], second_species = quaternary_species_list[1], third_species = quaternary_species_list[2], fourth_species = quaternary_species_list[3])
				self.quaternary_aw.show()
			
			
			
			print(compound_name)
			print(compound_type)



app = QApplication([])
ex = Initial_VTAnDeM_Dialog()
sys.exit(app.exec_())




