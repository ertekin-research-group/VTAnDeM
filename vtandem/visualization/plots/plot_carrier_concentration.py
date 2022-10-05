
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'
__author__ = 'Michael_Lidia_Jiaxing_Elif'

import numpy as np
import os
import matplotlib.pyplot as plt
from scipy import integrate
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

# Import functions for calculating carrier concentration
from vtandem.visualization.utils.carrier_concentration import Calculate_CarrierConcentration, Calculate_FreeHole_FreeElectron_Concentrations

from vtandem.visualization.plots.save_plot import SaveFigure



class Plot_CarrierConcentration(SaveFigure):
	
	def __init__(self, elements_list):
		
		# Font description for defect formation energy diagram
		self.font = { 'color': 'black', 'weight': 'normal', 'size': 12 }
		
		# Plot settings
		self.max_temperature = 1000  # Maximum temperature in Kelvins
		self.temperature_stepsize = 50

		# Store all extracted DFT data
		self.defects_data = None
		self.main_compound_info = None
		self.dos_data = None
		self.vol = 0.0	# Volume of defect supercell (NOT the DOS cell)
		self.EVBM = 0.0
		self.ECBM = 0.0
		self.fermi_energy_array = None
		self.temperature_array = np.arange(200, self.max_temperature+1, self.temperature_stepsize)
		self.synthesis_temperature = None
		
		self.energy = None
		self.gE = None
		self.energies_ValenceBand		= None
		self.gE_ValenceBand 			= None
		self.energies_ConductionBand	= None
		self.gE_ConductionBand 			= None
		
		self.intrinsic_equilibrium_fermi_energy = {}
		self.total_equilibrium_fermi_energy = {}
		for temperature in self.temperature_array:
			self.intrinsic_equilibrium_fermi_energy[temperature] = 0.0
			self.total_equilibrium_fermi_energy[temperature] = 0.0
		
		
		# Initialize all mu values
		self.mu_elements = {}
		for element in elements_list:
			self.mu_elements[element] = {"mu0": 0.0, "deltamu": 0.0}

		
		# Store user-selected dopant
		self.dopant = "None"
		self.dopant_mu0 = 0.0
		self.dopant_deltamu = 0.0
		self.extrinsic_defects = []  # List of extrinsic defects of dopant (e.g. Ge_Bi, Ge_Se, and Ge_O for Ge dopant in Bi2O2Se)
		
		# (WIDGET) Carrier Concentration Plot
		self.carrier_concentration_plot_figure = plt.figure()
		self.carrier_concentration_plot_figure.subplots_adjust(left=0.225)
		self.carrier_concentration_plot_drawing = self.carrier_concentration_plot_figure.add_subplot(111)
		self.carrier_concentration_plot_canvas = FigureCanvas(self.carrier_concentration_plot_figure)
		
		self.carrier_concentration_intrinsic_defect_hole_plot = None
		self.carrier_concentration_intrinsic_defect_electron_plot = None
		self.carrier_concentration_total_hole_plot = None
		self.carrier_concentration_total_electron_plot = None
		
		# Set y-axis minimum and maximum
		self.ymin = 1E16
		self.ymax = 1E23
		
		# Save figure feature
		SaveFigure.__init__(self, self.carrier_concentration_plot_figure)
	
	
	def Activate_CarrierConcentration_Plot_Axes(self):
		
		self.carrier_concentration_plot_drawing.set_xlim(200, self.max_temperature)
		self.carrier_concentration_plot_drawing.set_ylim(self.ymin, self.ymax)
		self.carrier_concentration_plot_drawing.set_xlabel("T(K)", fontdict=self.font)
		self.carrier_concentration_plot_drawing.set_ylabel("Carrier Concentration (cm$^{-3}$)", fontdict=self.font, rotation=90)
		self.carrier_concentration_plot_drawing.set_yscale("log")
		self.carrier_concentration_plot_drawing.xaxis.tick_bottom()
		self.carrier_concentration_plot_drawing.yaxis.tick_left()
		self.carrier_concentration_plot_drawing.tick_params(axis='both', labelsize=self.font['size']-2)
		self.carrier_concentration_plot_drawing.xaxis.set_label_position("bottom")
		self.carrier_concentration_plot_drawing.yaxis.set_label_position("left")
		self.carrier_concentration_plot_drawing.set_aspect("auto")
	
	
	def Update_WindowSize(self, ytype, Ylim_box_object):
		
		# Modify defects diagram y-axis
		if ytype == "YMin":
			self.ymin = float(Ylim_box_object.text())
		if ytype == "YMax":
			self.ymax = float(Ylim_box_object.text())
		self.carrier_concentration_plot_drawing.set_ylim(self.ymin, self.ymax)
		self.carrier_concentration_plot_canvas.draw()

	
	def Organize_DOS_Data(self):
		
		# Initialize data
		energy = []
		gE = []
		
		# Orgnize data
		for energy_dos in sorted(np.asarray([float(i) for i in self.dos_data["DOS"].keys()])):
			energy.append(energy_dos)
			gE.append(float(self.dos_data["DOS"][str(energy_dos)]))
		
		# Store into global variables
		self.energy = np.asarray(energy)
		self.gE = np.asarray(gE)
	

	
	def Extract_Relevant_Energies_DOSs(self):
		
		# Initialize data
		energies_ValenceBand = []
		gE_ValenceBand = []
		energies_ConductionBand = []
		gE_ConductionBand = []
		
		# The DOS band gap may not be the band gap for the defect formation energy
		#	diagram, especially when band gap corrections are applied. To mitigate
		#	this problem, we use a scissor operator where the VBM and CBM in the
		#	DOSCAR file are repositioned to the band gap of the defect formation
		#	energy diagram.
		past_dos_bandgap = False
		for energy, gE in zip(self.energy, self.gE):
			if energy <= 0.0:
				# Get all energies and corresponding DOSs below VBM
				energies_ValenceBand.append(energy)
				gE_ValenceBand.append(gE)
			else:
				# Get all energies and corresponding DOSs above CBM
				if (not past_dos_bandgap) and (gE <= 1E-4):
					continue
				else:
					past_dos_bandgap = True
					energies_ConductionBand.append(energy)
					gE_ConductionBand.append(gE)
		
		# Make data into numpy arrays
		self.energies_ValenceBand = np.asarray(energies_ValenceBand)
		self.gE_ValenceBand = np.asarray(gE_ValenceBand)
		self.energies_ConductionBand = np.asarray(energies_ConductionBand)
		self.gE_ConductionBand = np.asarray(gE_ConductionBand)
		
		# Reposition band edges to corrected values (NOT ZERO-ED)
		self.energies_ValenceBand += self.EVBM
		self.energies_ConductionBand += self.ECBM - np.min(self.energies_ConductionBand)
		
		# Normalize DOS to be per volume
		self.gE_ValenceBand /= self.dos_data["Volume"]
		self.gE_ConductionBand /= self.dos_data["Volume"]
	

	
	def Update_Deltamus(self, deltamu_values):

		# Args:	
		# 	deltamu_values: Dictionary of deltamu values, in element:value pairs
		for element in self.mu_elements.keys():
			self.mu_elements[element]["deltamu"] = deltamu_values[element]




	# Free carrier concentrations are calculated separately from defect concentrations. This is to prevent
	#	having to calculate them repeatedly for different thermodynamic conditions (delta mu values) since
	#	they're the same in each condition.
	def Calculate_Hole_Electron_Concentration_Matrices(self):
		self.hole_concentrations_dict, self.electron_concentrations_dict = Calculate_FreeHole_FreeElectron_Concentrations(	self.temperature_array, \
																															self.fermi_energy_array, \
																															self.gE_ValenceBand, \
																															self.energies_ValenceBand, \
																															self.gE_ConductionBand, \
																															self.energies_ConductionBand )
	
	

	def Initialize_CarrierConcentration_Plot(self):
		
		intrinsic_defect_hole_concentration, intrinsic_defect_electron_concentration, total_hole_concentration, total_electron_concentration, intrinsic_equilibrium_fermi_energy_temperature, total_equilibrium_fermi_energy_temperature = Calculate_CarrierConcentration(	EVBM = self.EVBM, \
																																																																			ECBM = self.ECBM, \
																																																																			energies_ValenceBand = self.energies_ValenceBand, \
																																																																			gE_ValenceBand = self.gE_ValenceBand, \
																																																																			energies_ConductionBand = self.energies_ConductionBand, \
																																																																			gE_ConductionBand = self.gE_ConductionBand, \
																																																																			defects_data = self.defects_data, \
																																																																			main_compound_info = self.main_compound_info, \
																																																																			mu_elements = self.mu_elements, \
																																																																			temperature_array = self.temperature_array, \
																																																																			fermi_energy_array = self.fermi_energy_array, \
																																																																			volume = self.vol, \
																																																																			extrinsic_defects = self.extrinsic_defects, \
																																																																			dopant = self.dopant, \
																																																																			dopant_mu0 = self.dopant_mu0, \
																																																																			dopant_deltamu = self.dopant_deltamu, \
																																																																			hole_concentrations_dict = self.hole_concentrations_dict, \
																																																																			electron_concentrations_dict = self.electron_concentrations_dict, \
																																																																			synthesis_temperature = self.synthesis_temperature )
		
		# Update equilibrium Fermi energy
		self.intrinsic_equilibrium_fermi_energy = intrinsic_equilibrium_fermi_energy_temperature
		self.total_equilibrium_fermi_energy = total_equilibrium_fermi_energy_temperature
		
		try:
			self.carrier_concentration_intrinsic_defect_hole_plot.remove()
			self.carrier_concentration_total_hole_plot.remove()
		except:
			pass

		try:
			self.carrier_concentration_intrinsic_defect_electron_plot.remove()
			self.carrier_concentration_total_electron_plot.remove()
		except:
			pass
		
		self.carrier_concentration_intrinsic_defect_hole_plot, = self.carrier_concentration_plot_drawing.semilogy(self.temperature_array, intrinsic_defect_hole_concentration, 'o-', color='red', label='Hole')
		if self.dopant != "None":
			self.carrier_concentration_total_hole_plot, = self.carrier_concentration_plot_drawing.semilogy(self.temperature_array, total_hole_concentration, 'o-', markerfacecolor='none', markeredgecolor='red', color='red', ls='--', label='Hole (With Dopant)')
		
		self.carrier_concentration_intrinsic_defect_electron_plot, = self.carrier_concentration_plot_drawing.semilogy(self.temperature_array, intrinsic_defect_electron_concentration, 'o-', color='green', label='Electron')
		if self.dopant != "None":
			self.carrier_concentration_total_electron_plot, = self.carrier_concentration_plot_drawing.semilogy(self.temperature_array, total_electron_concentration, 'o-', markerfacecolor='none', markeredgecolor='green', color='green', ls='--', label='Electron (With Dopant)')
		
		self.carrier_concentration_plot_drawing.legend(loc=1, fontsize=self.font['size'])
		self.carrier_concentration_plot_canvas.draw()
	


	def Update_CarrierConcentration_Plot(self):
		
		intrinsic_defect_hole_concentration, intrinsic_defect_electron_concentration, total_hole_concentration, total_electron_concentration, intrinsic_equilibrium_fermi_energy_temperature, total_equilibrium_fermi_energy_temperature = Calculate_CarrierConcentration(	EVBM = self.EVBM, \
																																																																			ECBM = self.ECBM, \
																																																																			energies_ValenceBand = self.energies_ValenceBand, \
																																																																			gE_ValenceBand = self.gE_ValenceBand, \
																																																																			energies_ConductionBand = self.energies_ConductionBand, \
																																																																			gE_ConductionBand = self.gE_ConductionBand, \
																																																																			defects_data = self.defects_data, \
																																																																			main_compound_info = self.main_compound_info, \
																																																																			mu_elements = self.mu_elements, \
																																																																			temperature_array = self.temperature_array, \
																																																																			fermi_energy_array = self.fermi_energy_array, \
																																																																			volume = self.vol, \
																																																																			extrinsic_defects = self.extrinsic_defects, \
																																																																			dopant = self.dopant, \
																																																																			dopant_mu0 = self.dopant_mu0, \
																																																																			dopant_deltamu = self.dopant_deltamu, \
																																																																			hole_concentrations_dict = self.hole_concentrations_dict, \
																																																																			electron_concentrations_dict = self.electron_concentrations_dict, \
																																																																			synthesis_temperature = self.synthesis_temperature )
		
		# Update equilibrium Fermi energy
		self.intrinsic_equilibrium_fermi_energy = intrinsic_equilibrium_fermi_energy_temperature
		self.total_equilibrium_fermi_energy = total_equilibrium_fermi_energy_temperature
		
		self.carrier_concentration_intrinsic_defect_hole_plot.set_ydata(intrinsic_defect_hole_concentration)
		if self.dopant != "None":
			self.carrier_concentration_total_hole_plot.set_ydata(total_hole_concentration)
		
		self.carrier_concentration_intrinsic_defect_electron_plot.set_ydata(intrinsic_defect_electron_concentration)
		if self.dopant != "None":
			self.carrier_concentration_total_electron_plot.set_ydata(total_electron_concentration)
		
		self.carrier_concentration_plot_canvas.draw()





