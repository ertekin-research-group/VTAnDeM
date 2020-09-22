
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'
__author__ = 'Michael_Lidia_Jiaxing_Elif'

import numpy as np
import os
import matplotlib.pyplot as plt
from scipy import integrate
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

# Import functions for calculating carrier concentration
from vtandem.visualization.carrier_concentration import Calculate_CarrierConcentration

class CarrierConcentration_Plot:
	
	def __init__(self):
		
		# Font description for defect formation energy diagram
		self.font = {'family': 'sans-serif',
				'color':  'black',
				'weight': 'normal',
				'size': 14 }
		
		# Store all extracted DFT data
		self.defects_data = None
		self.dos_data = None
		self.vol = 0.0
		self.EVBM = 0.0
		self.ECBM = 0.0
		self.fermi_energy_array = None
		self.temperature_array = np.arange(200, 801, 50)
		self.synthesis_temperature = None
		
		self.energy = None
		self.gE = None
		self.energies_ValenceBand		= None
		self.gE_ValenceBand 			= None
		self.energies_ConductionBand	= None
		self.gE_ConductionBand 			= None
		
		self.hole_concentrations_dict		= {}
		self.electron_concentrations_dict	= {}
		
		self.intrinsic_equilibrium_fermi_energy = {}
		self.total_equilibrium_fermi_energy = {}
		for temperature in self.temperature_array:
			self.intrinsic_equilibrium_fermi_energy[temperature] = 0.0
			self.total_equilibrium_fermi_energy[temperature] = 0.0
		
		
		self.check_outside_bandgap = False
		
		self.k = 8.6173303E-5
		
		
		# Store user-selected extrinsic defect
		self.extrinsic_defect = "None"
		self.extrinsic_defect_mu0 = 0.0
		self.extrinsic_defect_deltamu = 0.0
		
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
	
	
	def Activate_CarrierConcentration_Plot_Axes(self):
		
		self.carrier_concentration_plot_drawing.set_xlim(200, 800)
		self.carrier_concentration_plot_drawing.set_ylim(self.ymin, self.ymax)
		self.carrier_concentration_plot_drawing.set_xlabel("T(K)", fontdict=self.font)
		self.carrier_concentration_plot_drawing.set_ylabel("n$_i$ (cm$^{-3}$)", fontdict=self.font, rotation=90)
		self.carrier_concentration_plot_drawing.set_yscale("log")
		self.carrier_concentration_plot_drawing.xaxis.tick_bottom()
		self.carrier_concentration_plot_drawing.yaxis.tick_left()
		self.carrier_concentration_plot_drawing.tick_params(axis='both', labelsize=9)
		self.carrier_concentration_plot_drawing.xaxis.set_label_position("bottom")
		self.carrier_concentration_plot_drawing.yaxis.set_label_position("left")
		self.carrier_concentration_plot_drawing.set_aspect("auto")
	
	
	def Organize_DOS_Data(self):
		
		# Initialize data
		energy = []
		gE = []
		
		# Orgnize data
		for energy_dos in sorted(np.asarray([float(i) for i in self.dos_data.keys()])):
			energy.append(energy_dos)
			gE.append(float(self.dos_data[str(energy_dos)])) #*1E24)
		
		# Store into global variables
		self.energy = np.asarray(energy)
		self.gE = np.asarray(gE)
	
	
	def Extract_Relevant_Energies_DOSs(self):
		
		# Initialize data
		energies_ValenceBand = []
		gE_ValenceBand = []
		energies_ConductionBand = []
		gE_ConductionBand = []
		
		for energy, gE in zip(self.energy, self.gE):
			# Get all energies and corresponding DOSs below VBM
			if energy <= self.EVBM:
				energies_ValenceBand.append(energy)
				gE_ValenceBand.append(gE)
			# Get all energies and corresponding DOSs above CBM
			if energy >= self.ECBM:
				energies_ConductionBand.append(energy)
				gE_ConductionBand.append(gE)
		
		# Make data into numpy arrays
		self.energies_ValenceBand = np.asarray(energies_ValenceBand)
		self.gE_ValenceBand = np.asarray(gE_ValenceBand)
		self.energies_ConductionBand = np.asarray(energies_ConductionBand)
		self.gE_ConductionBand = np.asarray(gE_ConductionBand)
	
	
	def Calculate_Hole_Electron_Concentration_Matrices(self):
		
		for temperature in self.temperature_array:
			
			self.hole_concentrations_dict[temperature] = np.zeros(len(self.fermi_energy_array))
			self.electron_concentrations_dict[temperature] = np.zeros(len(self.fermi_energy_array))
			
			for ef in self.fermi_energy_array:
				
				# Hole concentration
				fE_holes = self.gE_ValenceBand * (1. - 1./( 1. + np.exp( (self.energies_ValenceBand - ef) / (self.k * temperature) ) ) )
				hole_concentration = integrate.simps(fE_holes, self.energies_ValenceBand)
				self.hole_concentrations_dict[temperature][self.fermi_energy_array.tolist().index(ef)] = hole_concentration
				
				# Electron concentration
				fE_electrons = self.gE_ConductionBand / (1. + np.exp( (self.energies_ConductionBand - ef) / (self.k * temperature) ) )
				electron_concentration = integrate.simps(fE_electrons, self.energies_ConductionBand)
				self.electron_concentrations_dict[temperature][self.fermi_energy_array.tolist().index(ef)] = electron_concentration
	
	
	def Initialize_HoleConcentration_Plot(self):
		
		intrinsic_defect_hole_concentration, intrinsic_defect_electron_concentration, total_hole_concentration, total_electron_concentration, intrinsic_equilibrium_fermi_energy_temperature, total_equilibrium_fermi_energy_temperature = Calculate_CarrierConcentration(	EVBM = self.EVBM, \
																																																																			ECBM = self.ECBM, \
																																																																			energies_ValenceBand = self.energies_ValenceBand, \
																																																																			gE_ValenceBand = self.gE_ValenceBand, \
																																																																			energies_ConductionBand = self.energies_ConductionBand, \
																																																																			gE_ConductionBand = self.gE_ConductionBand, \
																																																																			defects_data = self.defects_data, \
																																																																			main_compound_total_energy = self.main_compound_total_energy, \
																																																																			mu_elements = self.mu_elements, \
																																																																			temperature_array = self.temperature_array, \
																																																																			fermi_energy_array = self.fermi_energy_array, \
																																																																			volume = self.vol, \
																																																																			number_species = self.number_species, \
																																																																			extrinsic_defect = self.extrinsic_defect, \
																																																																			extrinsic_defect_mu0 = self.extrinsic_defect_mu0, \
																																																																			extrinsic_defect_deltamu = self.extrinsic_defect_deltamu, \
																																																																			hole_concentrations_dict = self.hole_concentrations_dict, \
																																																																			electron_concentrations_dict = self.electron_concentrations_dict, \
																																																																			check_outside_bandgap = self.check_outside_bandgap, \
																																																																			synthesis_temperature = self.synthesis_temperature )
		
		# Update equilibrium Fermi energy
		self.intrinsic_equilibrium_fermi_energy = intrinsic_equilibrium_fermi_energy_temperature
		self.total_equilibrium_fermi_energy = total_equilibrium_fermi_energy_temperature
		
		try:
			self.carrier_concentration_intrinsic_defect_hole_plot.remove()
			self.carrier_concentration_total_hole_plot.remove()
		except:
			pass
		
		self.carrier_concentration_intrinsic_defect_hole_plot, = self.carrier_concentration_plot_drawing.semilogy(self.temperature_array, intrinsic_defect_hole_concentration, 'o-', color='red', label='Hole')
		if self.extrinsic_defect != "None":
			self.carrier_concentration_total_hole_plot, = self.carrier_concentration_plot_drawing.semilogy(self.temperature_array, total_hole_concentration, 'o-', markerfacecolor='none', markeredgecolor='red', color='red', ls='--', label='Hole (With Dopant)')
		
		self.carrier_concentration_plot_drawing.legend(loc=1)
		self.carrier_concentration_plot_canvas.draw()
	
	
	
	def Update_HoleConcentration_Plot(self):
		
		intrinsic_defect_hole_concentration, intrinsic_defect_electron_concentration, total_hole_concentration, total_electron_concentration, intrinsic_equilibrium_fermi_energy_temperature, total_equilibrium_fermi_energy_temperature = Calculate_CarrierConcentration(	EVBM = self.EVBM, \
																																																																			ECBM = self.ECBM, \
																																																																			energies_ValenceBand = self.energies_ValenceBand, \
																																																																			gE_ValenceBand = self.gE_ValenceBand, \
																																																																			energies_ConductionBand = self.energies_ConductionBand, \
																																																																			gE_ConductionBand = self.gE_ConductionBand, \
																																																																			defects_data = self.defects_data, \
																																																																			main_compound_total_energy = self.main_compound_total_energy, \
																																																																			mu_elements = self.mu_elements, \
																																																																			temperature_array = self.temperature_array, \
																																																																			fermi_energy_array = self.fermi_energy_array, \
																																																																			volume = self.vol, \
																																																																			number_species = self.number_species, \
																																																																			extrinsic_defect = self.extrinsic_defect, \
																																																																			extrinsic_defect_mu0 = self.extrinsic_defect_mu0, \
																																																																			extrinsic_defect_deltamu = self.extrinsic_defect_deltamu, \
																																																																			hole_concentrations_dict = self.hole_concentrations_dict, \
																																																																			electron_concentrations_dict = self.electron_concentrations_dict, \
																																																																			check_outside_bandgap = self.check_outside_bandgap, \
																																																																			synthesis_temperature = self.synthesis_temperature )
		
		# Update equilibrium Fermi energy
		self.intrinsic_equilibrium_fermi_energy = intrinsic_equilibrium_fermi_energy_temperature
		self.total_equilibrium_fermi_energy = total_equilibrium_fermi_energy_temperature
		
		self.carrier_concentration_intrinsic_defect_hole_plot.set_ydata(intrinsic_defect_hole_concentration)
		if self.extrinsic_defect != "None":
			self.carrier_concentration_total_hole_plot.set_ydata(total_hole_concentration)
		
		self.carrier_concentration_plot_canvas.draw()
	
	
	
	def Initialize_ElectronConcentration_Plot(self):
		
		intrinsic_defect_hole_concentration, intrinsic_defect_electron_concentration, total_hole_concentration, total_electron_concentration, intrinsic_equilibrium_fermi_energy_temperature, total_equilibrium_fermi_energy_temperature = Calculate_CarrierConcentration(	EVBM = self.EVBM, \
																																																																			ECBM = self.ECBM, \
																																																																			energies_ValenceBand = self.energies_ValenceBand, \
																																																																			gE_ValenceBand = self.gE_ValenceBand, \
																																																																			energies_ConductionBand = self.energies_ConductionBand, \
																																																																			gE_ConductionBand = self.gE_ConductionBand, \
																																																																			defects_data = self.defects_data, \
																																																																			main_compound_total_energy = self.main_compound_total_energy, \
																																																																			mu_elements = self.mu_elements, \
																																																																			temperature_array = self.temperature_array, \
																																																																			fermi_energy_array = self.fermi_energy_array, \
																																																																			volume = self.vol, \
																																																																			number_species = self.number_species, \
																																																																			extrinsic_defect = self.extrinsic_defect, \
																																																																			extrinsic_defect_mu0 = self.extrinsic_defect_mu0, \
																																																																			extrinsic_defect_deltamu = self.extrinsic_defect_deltamu, \
																																																																			hole_concentrations_dict = self.hole_concentrations_dict, \
																																																																			electron_concentrations_dict = self.electron_concentrations_dict, \
																																																																			check_outside_bandgap = self.check_outside_bandgap, \
																																																																			synthesis_temperature = self.synthesis_temperature )
		
		# Update equilibrium Fermi energy
		self.intrinsic_equilibrium_fermi_energy = intrinsic_equilibrium_fermi_energy_temperature
		self.total_equilibrium_fermi_energy = total_equilibrium_fermi_energy_temperature
		
		try:
			self.carrier_concentration_intrinsic_defect_electron_plot.remove()
			self.carrier_concentration_total_electron_plot.remove()
		except:
			pass
		
		self.carrier_concentration_intrinsic_defect_electron_plot, = self.carrier_concentration_plot_drawing.semilogy(self.temperature_array, intrinsic_defect_electron_concentration, 'o-', color='green', label='Electron')
		if self.extrinsic_defect != "None":
			self.carrier_concentration_total_electron_plot, = self.carrier_concentration_plot_drawing.semilogy(self.temperature_array, total_electron_concentration, 'o-', markerfacecolor='none', markeredgecolor='green', color='green', ls='--', label='Electron (With Dopant)')
		
		self.carrier_concentration_plot_drawing.legend(loc=1)
		self.carrier_concentration_plot_canvas.draw()
	
	
	
	def Update_ElectronConcentration_Plot(self):
		
		intrinsic_defect_hole_concentration, intrinsic_defect_electron_concentration, total_hole_concentration, total_electron_concentration, intrinsic_equilibrium_fermi_energy_temperature, total_equilibrium_fermi_energy_temperature = Calculate_CarrierConcentration(	EVBM = self.EVBM, \
																																																																			ECBM = self.ECBM, \
																																																																			energies_ValenceBand = self.energies_ValenceBand, \
																																																																			gE_ValenceBand = self.gE_ValenceBand, \
																																																																			energies_ConductionBand = self.energies_ConductionBand, \
																																																																			gE_ConductionBand = self.gE_ConductionBand, \
																																																																			defects_data = self.defects_data, \
																																																																			main_compound_total_energy = self.main_compound_total_energy, \
																																																																			mu_elements = self.mu_elements, \
																																																																			temperature_array = self.temperature_array, \
																																																																			fermi_energy_array = self.fermi_energy_array, \
																																																																			volume = self.vol, \
																																																																			number_species = self.number_species, \
																																																																			extrinsic_defect = self.extrinsic_defect, \
																																																																			extrinsic_defect_mu0 = self.extrinsic_defect_mu0, \
																																																																			extrinsic_defect_deltamu = self.extrinsic_defect_deltamu, \
																																																																			hole_concentrations_dict = self.hole_concentrations_dict, \
																																																																			electron_concentrations_dict = self.electron_concentrations_dict, \
																																																																			check_outside_bandgap = self.check_outside_bandgap,
																																																																			synthesis_temperature = self.synthesis_temperature )
		
		# Update equilibrium Fermi energy
		self.intrinsic_equilibrium_fermi_energy = intrinsic_equilibrium_fermi_energy_temperature
		self.total_equilibrium_fermi_energy = total_equilibrium_fermi_energy_temperature
		
		self.carrier_concentration_intrinsic_defect_electron_plot.set_ydata(intrinsic_defect_electron_concentration)
		if self.extrinsic_defect != "None":
			self.carrier_concentration_total_electron_plot.set_ydata(total_electron_concentration)
		
		self.carrier_concentration_plot_canvas.draw()











