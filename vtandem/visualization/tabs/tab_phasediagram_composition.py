
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'
__author__ = 'Michael_Lidia_Jiaxing_Elif'

import numpy as np

import PyQt5
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

from vtandem.visualization.windows.window_phasediagram_composition import Window_Compositional_PhaseDiagram
from vtandem.visualization.windows.window_defectsdiagram import Window_DefectsDiagram
from vtandem.visualization.windows.window_carrierconcentration import Window_CarrierConcentration

class Tab_Compositional_PhaseDiagram(Window_DefectsDiagram, Window_CarrierConcentration):
	
	def __init__(self, 	type: str, \
						compounds_info = None, \
						defects_data = None, \
						main_compound_info = None, \
						dos_data = None, \
						show_defects_diagram = True, \
						show_carrier_concentration = True ):	# User specifies the main compound and its constituents
		
		# Check that the 'type' argument is legitimate
		if (type != "ternary") and (type != "quaternary"):
			raise Exception("The argument 'type' must be either 'ternary' or 'quaternary'. Exiting...")
		self.type = type
		
		# Display settings
		self.show_defects_diagram = show_defects_diagram
		self.show_carrier_concentration = show_carrier_concentration

		
		# Compositional phase diagram object is created in tab_quaternary_.../tab_ternary_...
		
		# Defects diagram
		if self.show_defects_diagram:
			
			self.DefectsDiagram.defects_data = defects_data
			self.DefectsDiagram.main_compound_info = main_compound_info
			
			for element in self.elements_list:
				self.DefectsDiagram.mu_elements[element]["mu0"] = compounds_info[element]["mu0"]
				self.DefectsDiagram.mu_elements[element]["deltamu"] = self.Compositional_PhaseDiagram.deltamu_values[element]
			
			self.DefectsDiagram.EVBM = main_compound_info["VBM"]
			self.DefectsDiagram.ECBM = self.DefectsDiagram.EVBM + main_compound_info["BandGap"]
			self.DefectsDiagram.axis_lims["XMin"] = 0.0
			self.DefectsDiagram.axis_lims["XMax"] = main_compound_info["BandGap"]
			self.DefectsDiagram.fermi_energy_array = np.linspace(self.DefectsDiagram.EVBM, self.DefectsDiagram.ECBM, 2000)
			self.DefectsDiagram.Activate_DefectsDiagram_Plot_Axes()
			
			# Update defects diagram in response to clicking on phase diagram
			self.Compositional_PhaseDiagram.composition_phasediagram_plot_figure.canvas.mpl_connect('button_press_event', self.Update_DefectsDiagram_Plot_Function)
		
		
		# Carrier concentration
		if self.show_carrier_concentration:
			
			self.CarrierConcentration.defects_data = defects_data
			self.CarrierConcentration.main_compound_info = main_compound_info
			self.CarrierConcentration.dos_data = dos_data[self.main_compound]
			
			for element in self.elements_list:
				self.CarrierConcentration.mu_elements[element]["mu0"] = compounds_info[element]["mu0"]
				self.CarrierConcentration.mu_elements[element]["deltamu"] = self.Compositional_PhaseDiagram.deltamu_values[element]
			
			self.CarrierConcentration.vol = main_compound_info["Volume"]
			self.CarrierConcentration.EVBM = self.DefectsDiagram.EVBM
			self.CarrierConcentration.ECBM = self.DefectsDiagram.ECBM
			# Fermi energies should sample outside band gap, in case EFeq is not in gap
			self.CarrierConcentration.fermi_energy_array = np.linspace(self.DefectsDiagram.EVBM-1.0, self.DefectsDiagram.ECBM+1.0, 2000)
			self.CarrierConcentration.Activate_CarrierConcentration_Plot_Axes()
			self.CarrierConcentration.Organize_DOS_Data()
			self.CarrierConcentration.Extract_Relevant_Energies_DOSs()
			self.CarrierConcentration.Calculate_Hole_Electron_Concentration_Matrices()
			
			# Update carrier concentration in response to clicking on phase diagram
			self.Compositional_PhaseDiagram.composition_phasediagram_plot_figure.canvas.mpl_connect('button_press_event', self.Update_CarrierConcentration_Plot_Function)
		
		
		###############################################################################################
		###############################################################################################
		#################################### Initialize third tab #####################################
		###############################################################################################
		###############################################################################################
		
		self.tab3 = QWidget()
		self.tab3_layout = QHBoxLayout(self.tab3)
		
		# Add compositional phase diagram window widget to tab3
		self.Compositional_PhaseDiagram_Window = Window_Compositional_PhaseDiagram(self.main_compound, self.Compositional_PhaseDiagram)
		self.Compositional_PhaseDiagram_Window.compositional_phasediagram_window_layout.addWidget(self.Compositional_PhaseDiagram_Window.composition_phase_diagram_title)
		self.Compositional_PhaseDiagram_Window.compositional_phasediagram_window_layout.addWidget(self.Compositional_PhaseDiagram_Window.composition_phase_diagram_plot)
		self.Compositional_PhaseDiagram_Window.compositional_phasediagram_window_layout.addWidget(self.Compositional_PhaseDiagram_Window.phasediagram_savefigure_button)
		self.tab3_layout.addWidget(self.Compositional_PhaseDiagram_Window.compositional_phasediagram_window)
		
		# Add defect formation energy diagram window widget to tab3
		if self.show_defects_diagram:
			Window_DefectsDiagram.__init__(self, show_dopant=False)
			self.tab3_layout.addWidget(self.defectsdiagram_window)
		

		if self.show_carrier_concentration:
			Window_CarrierConcentration.__init__(self)
			self.tab3_layout.addWidget(self.carrierconcentration_window)
	
	
	
	###############################################################################################
	################################# Generate Defects Diagram ####################################
	###############################################################################################
	
	def Update_DefectsDiagram_Plot_Function(self, event):
		print("Defects Updating...")

		# Reset defects diagram
		self.DefectsDiagram.defects_diagram_plot_drawing.remove()
		self.DefectsDiagram.defects_diagram_plot_drawing = self.DefectsDiagram.defects_diagram_plot_figure.add_subplot(111)
		self.DefectsDiagram.Activate_DefectsDiagram_Plot_Axes()
		self.DefectsDiagram.defects_diagram_plot_canvas.draw()

		# Check that phase region has been selected
		if self.Compositional_PhaseDiagram.phaseregion_selected is None:
			return
		
		# Update elements and chemical potentials
		for element in self.elements_list:
			self.DefectsDiagram.mu_elements[element]["deltamu"] = self.Compositional_PhaseDiagram.deltamu_values[element]
		
		# Calculate defect formation energies
		self.DefectsDiagram.Calculate_DefectFormations()
		
		# Plot defect formation energies
		self.DefectsDiagram.intrinsic_defect_plots = {}
		self.DefectsDiagram.extrinsic_defect_plots = {}
		self.DefectsDiagram.Initialize_Intrinsic_DefectsDiagram_Plot()
	
	
	
	###############################################################################################
	############################### Generate Carrier Concentration ################################
	###############################################################################################
	
	def Update_CarrierConcentration_Plot_Function(self, event):
		
		# Reset carrier concentration
		self.CarrierConcentration.carrier_concentration_plot_drawing.remove()
		self.CarrierConcentration.carrier_concentration_plot_drawing = self.CarrierConcentration.carrier_concentration_plot_figure.add_subplot(111)
		self.CarrierConcentration.Activate_CarrierConcentration_Plot_Axes()
		self.CarrierConcentration.carrier_concentration_plot_canvas.draw()
		
		# Check that phase region has been selected
		if self.Compositional_PhaseDiagram.phaseregion_selected is None:
			# NOTE: These objects only exists when self.show_carrier_concentration = True
			self.defects_synthesis_temperature_box.setEnabled(False)
			self.defectconc_stat_box.setEnabled(False)
			self.temperature_selection_box.setEnabled(False)
			self.dos_option_box.setEnabled(False)
			return
		else:
			self.defects_synthesis_temperature_box.setEnabled(True)
			self.defectconc_stat_box.setEnabled(True)
			self.temperature_selection_box.setEnabled(True)
			self.dos_option_box.setEnabled(True)
		
		# Update elements and chemical potentials
		for element in self.elements_list:
			self.CarrierConcentration.mu_elements[element]["deltamu"] = self.Compositional_PhaseDiagram.deltamu_values[element]

		# Plot the carrier concentration (holes and electrons)
		self.CarrierConcentration.carrier_concentration_intrinsic_defect_hole_plot = None
		self.CarrierConcentration.carrier_concentration_intrinsic_defect_electron_plot = None
		self.CarrierConcentration.carrier_concentration_total_hole_plot = None
		self.CarrierConcentration.carrier_concentration_total_electron_plot = None
		self.CarrierConcentration.Initialize_CarrierConcentration_Plot()
		
		# Plot the equilibrium Fermi energy
		if self.DefectsDiagram.intrinsic_defect_plots != {}:
			self.Update_Equilibrium_Fermi_Energy_Temperature()






