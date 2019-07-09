
__author__ = 'Michael_Lidia_Jiaxing_Benita_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'


import numpy as np
import os
import matplotlib.pyplot as plt
from scipy import integrate
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

import PyQt5
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *


class Ternary_Carrier_Concentration(object):
	
	def __init__(self, parent = None, main_compound = None, first_element = None, second_element = None, third_element = None):
		
		
		# Font description for phase stability diagram plot
		self.font = {'family': 'sans-serif',
				'color':  'black',
				'weight': 'normal',
				'size': 14 }
		
		# Establish the first, second, third, and fourth species of the quaternary compound.
		# Note that this list is subject to change, depending on what the user chooses.
		self.main_compound  = main_compound
		self.first_element	= first_element
		self.second_element	= second_element
		self.third_element	= third_element
		self.elements_list  = [self.first_element, self.second_element, self.third_element]
		
		
		# Keep track of chemical potential values
		self.mu_elements = {self.first_element: {"mu0": 0.0, "deltamu": 0.0},
							self.second_element: {"mu0": 0.0, "deltamu": 0.0},
							self.third_element: {"mu0": 0.0, "deltamu": 0.0} }
		
		"""
		# Establish constants for the species
		self.species_a = first_element
		self.species_b = second_element
		self.species_c = third_element
		"""
		
		# Number of each specie in the main ternary compound
		self.main_compound_number_first_specie  = 0
		self.main_compound_number_second_specie = 0
		self.main_compound_number_third_specie  = 0
		
		
		# Store all extracted DFT data
		self.ternary_defects_data = None
		self.ternary_dos_data = None
		self.main_compound_total_energy = 0.0
		self.first_element_mu0 = 0.0
		self.second_element_mu0 = 0.0
		self.third_element_mu0 = 0.0
		self.EVBM = 0.0
		self.ECBM = 0.0
		self.fermi_energy_array = None
		self.temperature_array = np.arange(200, 801, 50)
		
		"""
		# Store all mu values
		self.mu_values = {}
		self.mu_values[self.first_element]  = 0.0
		self.mu_values[self.second_element] = 0.0
		self.mu_values[self.third_element]  = 0.0
		"""
		
		self.energy = None
		self.gE = None
		self.energy_EVBM	= None
		self.gE_EVBM 		= None
		self.energy_ECBM	= None
		self.gE_ECBM 		= None
		
		
		self.hole_concentrations_dict		= {}
		self.electron_concentrations_dict	= {}
		
		"""
		self.doscar = np.loadtxt("DOSCAR_Hg2GeTe4", skiprows=6)
		self.energy = self.doscar[:,0]
		self.gE = self.doscar[:,1]*1E24
		"""
		
		
		self.intrinsic_equilibrium_fermi_energy = {}
		self.total_equilibrium_fermi_energy = {}
		for temperature in self.temperature_array:
			self.intrinsic_equilibrium_fermi_energy[temperature] = 0.0
			self.total_equilibrium_fermi_energy[temperature] = 0.0
		
		self.check_outside_bandgap = False
		
		self.k = 8.6173303E-5
		self.vol = 230.96781596310623E-24
		
		
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
		for energy_dos in sorted(np.asarray([float(i) for i in self.ternary_dos_data.keys()])):
			energy.append(energy_dos)
			gE.append(float(self.ternary_dos_data[str(energy_dos)])*1E24)
		
		# Store into global variables
		self.energy = np.asarray(energy)
		self.gE = np.asarray(gE)
	
	
	
	def Extract_Relevant_Energies_DOSs(self):
		
		# Initialize data
		energy_EVBM = []
		gE_EVBM = []
		energy_ECBM = []
		gE_ECBM = []
		
		for energy, gE in zip(self.energy, self.gE):
			# Get all energies and corresponding DOSs below VBM
			if energy <= self.EVBM:
				energy_EVBM.append(energy)
				gE_EVBM.append(gE)
			# Get all energies and corresponding DOSs above CBM
			if energy >= self.ECBM:
				energy_ECBM.append(energy)
				gE_ECBM.append(gE)
		
		# Make data into numpy arrays
		self.energy_EVBM = np.asarray(energy_EVBM)
		self.gE_EVBM = np.asarray(gE_EVBM)
		self.energy_ECBM = np.asarray(energy_ECBM)
		self.gE_ECBM = np.asarray(gE_ECBM)
	
	
	def Calculate_Hole_Electron_Concentration_Matrices(self):
		
		for temperature in self.temperature_array:
			
			self.hole_concentrations_dict[temperature] = np.zeros(len(self.fermi_energy_array))
			self.electron_concentrations_dict[temperature] = np.zeros(len(self.fermi_energy_array))
			
			for ef in self.fermi_energy_array:
				
				# Hole concentration
				fE_holes = self.gE_EVBM * (1. - 1./( 1. + np.exp( (self.energy_EVBM - ef) / (self.k * temperature) ) ) )
				hole_concentration = integrate.simps(fE_holes, self.energy_EVBM)
				self.hole_concentrations_dict[temperature][self.fermi_energy_array.tolist().index(ef)] = hole_concentration
				
				# Electron concentration
				fE_electrons = self.gE_ECBM / (1. + np.exp( (self.energy_ECBM - ef) / (self.k * temperature) ) )
				electron_concentration = integrate.simps(fE_electrons, self.energy_ECBM)
				self.electron_concentrations_dict[temperature][self.fermi_energy_array.tolist().index(ef)] = electron_concentration
	
	
	def Calculate_Defect_Carrier_Concentration(self):
		
		"""
		dmu_a = self.mu_values[self.species_a]
		dmu_b = self.mu_values[self.species_b]
		dmu_c = self.mu_values[self.species_c]
		"""
		
		# Initialize
		intrinsic_defect_carrier_concentration_temperature = {}
		extrinsic_defect_carrier_concentration_temperature = {}
		for temperature in self.temperature_array:
			intrinsic_defect_carrier_concentration_temperature[temperature] = np.zeros(len(self.fermi_energy_array))
			extrinsic_defect_carrier_concentration_temperature[temperature] = np.zeros(len(self.fermi_energy_array))
		
		# Obtain intrinsic defect carrier concentration
		for defect in self.ternary_defects_data.keys():
			
			# Check that item is truly a defect
			if "_" not in defect:
				continue
			
			if self.ternary_defects_data[defect]["Extrinsic"] == "No":
				
				for charge in self.ternary_defects_data[defect]["charge"].keys():
					
					# Formation energy
					"""
					defect_formation_enthalpy = self.ternary_defects_data[defect]["charge"][charge]["Energy"] \
												- self.ternary_defects_data["supercellsize"] * self.main_compound_total_energy \
												- self.ternary_defects_data[defect]["n_"+self.species_a] * ( self.first_element_mu0 + float(dmu_a) ) \
												- self.ternary_defects_data[defect]["n_"+self.species_b] * ( self.second_element_mu0 + float(dmu_b) ) \
												- self.ternary_defects_data[defect]["n_"+self.species_c] * ( self.third_element_mu0 + float(dmu_c) ) \
												+ float(charge) * self.fermi_energy_array \
												+ self.ternary_defects_data[defect]["charge"][charge]["ECorr"]
					"""
					defect_formation_enthalpy = self.ternary_defects_data[defect]["charge"][charge]["Energy"] \
												- self.ternary_defects_data["supercellsize"] * self.main_compound_total_energy \
												+ float(charge) * self.fermi_energy_array \
												+ self.ternary_defects_data[defect]["charge"][charge]["ECorr"]
					for element in self.elements_list:
						defect_formation_enthalpy -= self.ternary_defects_data[defect]["n_"+element] * ( self.mu_elements[element]["mu0"] + self.mu_elements[element]["deltamu"] )
					
					# Defect concentration
					if self.first_element == defect.split("_")[-1]:
						N = self.main_compound_number_first_specie/self.vol
					elif self.second_element == defect.split("_")[-1]:
						N = self.main_compound_number_second_specie/self.vol
					elif self.third_element == defect.split("_")[-1]:
						N = self.main_compound_number_third_specie/self.vol
					
					for temperature in self.temperature_array:
						defect_carrier_concentration = float(charge) * N * np.exp(-defect_formation_enthalpy / (self.k * temperature) )
						intrinsic_defect_carrier_concentration_temperature[temperature] += defect_carrier_concentration
		
		# Check that user-selected extrinsic defect is not "None"
		if self.extrinsic_defect == "None":
			return intrinsic_defect_carrier_concentration_temperature, extrinsic_defect_carrier_concentration_temperature
		
		# Obtain extrinsic defect carrier concentration
		for charge in self.ternary_defects_data[self.extrinsic_defect]["charge"].keys():
			
			# Formation energy
			"""
			extrinsic_defect_formation_enthalpy = self.ternary_defects_data[self.extrinsic_defect]["charge"][charge]["Energy"] \
												- self.ternary_defects_data["supercellsize"] * self.main_compound_total_energy \
												- self.ternary_defects_data[self.extrinsic_defect]["n_"+self.species_a] * ( self.first_element_mu0 + float(dmu_a) ) \
												- self.ternary_defects_data[self.extrinsic_defect]["n_"+self.species_b] * ( self.second_element_mu0 + float(dmu_b) ) \
												- self.ternary_defects_data[self.extrinsic_defect]["n_"+self.species_c] * ( self.third_element_mu0 + float(dmu_c) ) \
												- (self.extrinsic_defect_mu0 + self.extrinsic_defect_deltamu) \
												+ float(charge) * self.fermi_energy_array \
												+ self.ternary_defects_data[self.extrinsic_defect]["charge"][charge]["ECorr"]
			"""
			extrinsic_defect_formation_enthalpy = self.ternary_defects_data[self.extrinsic_defect]["charge"][charge]["Energy"] \
												- self.ternary_defects_data["supercellsize"] * self.main_compound_total_energy \
												- (self.extrinsic_defect_mu0 + self.extrinsic_defect_deltamu) \
												+ float(charge) * self.fermi_energy_array \
												+ self.ternary_defects_data[self.extrinsic_defect]["charge"][charge]["ECorr"]
			for element in self.elements_list:
				extrinsic_defect_formation_enthalpy -= self.ternary_defects_data[self.extrinsic_defect]["n_"+element] * ( self.mu_elements[element]["mu0"] + self.mu_elements[element]["deltamu"] )
			
			# Extrinsic defect concentration
			if self.first_element == defect.split("_")[-1]:
				N = self.main_compound_number_first_specie/self.vol
			elif self.second_element == defect.split("_")[-1]:
				N = self.main_compound_number_second_specie/self.vol
			elif self.third_element == defect.split("_")[-1]:
				N = self.main_compound_number_third_specie/self.vol
			
			for temperature in self.temperature_array:
				extrinsic_defect_carrier_concentration = float(charge) * N * np.exp(-extrinsic_defect_formation_enthalpy / (self.k * temperature) )
				extrinsic_defect_carrier_concentration_temperature[temperature] += extrinsic_defect_carrier_concentration
		
		return intrinsic_defect_carrier_concentration_temperature, extrinsic_defect_carrier_concentration_temperature		# Returns a dictionary with temperature as keys and array of defect-induced carrier concentrations (not defect concentrations) as values
	
	
	
	def Calculate_CarrierConcentration(self):
		
		# Carrier concentrations from intrinsic defects only
		intrinsic_defect_hole_concentration = []
		intrinsic_defect_electron_concentration = []
		
		# Carrier concentration from both intrinsic and extrinsic defects
		total_defect_hole_concentration = []
		total_defect_electron_concentration = []
		
		# Calculate defect carrier concentration (for intrinsic defects and extrinsic defects)
		intrinsic_defect_carrier_concentration, extrinsic_defect_carrier_concentration = self.Calculate_Defect_Carrier_Concentration()
		
		for temperature in self.temperature_array:
			
			# Charge density including intrinsic defects only
			intrinsic_defect_charge_density_array = intrinsic_defect_carrier_concentration[temperature] + self.hole_concentrations_dict[temperature] - self.electron_concentrations_dict[temperature]
			intrinsic_equilibrium_fermi_energy = 0.0
			intrinsic_equilibrium_fermi_energy_index = 0
			
			# Charge density including both intrinsic and extrinsic defects
			total_charge_density_array = intrinsic_defect_carrier_concentration[temperature] + extrinsic_defect_carrier_concentration[temperature] + self.hole_concentrations_dict[temperature] - self.electron_concentrations_dict[temperature]
			total_equilibrium_fermi_energy = 0.0
			total_equilibrium_fermi_energy_index = 0
			
			
			
			# Check if equilibrium Fermi energy can be found within band gap (for charge density including intrinsic defects only)
			if np.sign(intrinsic_defect_charge_density_array[0]) == np.sign(intrinsic_defect_charge_density_array[-1]):
				
				if self.check_outside_bandgap:
					
					intrinsic_equilibrium_fermi_energy, intrinsic_hole_concentration, intrinsic_electron_concentration = self.Check_Outside_BandGap(intrinsic_defect_charge_density_array, temperature, calculation_type="intrinsic")
					
					self.intrinsic_equilibrium_fermi_energy[temperature] = intrinsic_equilibrium_fermi_energy - self.EVBM
					intrinsic_defect_hole_concentration.append(intrinsic_hole_concentration)
					intrinsic_defect_electron_concentration.append(intrinsic_electron_concentration)
				
				else:
					if ( (intrinsic_defect_charge_density_array[0] < intrinsic_defect_charge_density_array[-1]) and (np.sign(intrinsic_defect_charge_density_array[0]) == 1) ) or ( (intrinsic_defect_charge_density_array[0] > intrinsic_defect_charge_density_array[-1]) and (np.sign(intrinsic_defect_charge_density_array[0]) == -1) ):
						self.intrinsic_equilibrium_fermi_energy[temperature] = "< EVBM"
					elif ( (intrinsic_defect_charge_density_array[0] > intrinsic_defect_charge_density_array[-1]) and (np.sign(intrinsic_defect_charge_density_array[0]) == 1) ) or ( (intrinsic_defect_charge_density_array[0] < intrinsic_defect_charge_density_array[-1]) and (np.sign(intrinsic_defect_charge_density_array[0]) == -1) ):
						self.intrinsic_equilibrium_fermi_energy[temperature] = "> ECBM"
					intrinsic_defect_hole_concentration.append(0.0)
					intrinsic_defect_electron_concentration.append(0.0)
			
			
			elif np.sign(intrinsic_defect_charge_density_array[0]) != np.sign(intrinsic_defect_charge_density_array[-1]):
			
				# Search for equilibrium Fermi energy within band gap of material
				for total_charge_density_index in range(len(self.fermi_energy_array)-1):
					
					# Find equilibrium Fermi energy from charge density including intrinsic defects only
					if np.sign(intrinsic_defect_charge_density_array[total_charge_density_index]) != np.sign(intrinsic_defect_charge_density_array[total_charge_density_index+1]):
						intrinsic_equilibrium_fermi_energy = self.fermi_energy_array[total_charge_density_index]
						intrinsic_equilibrium_fermi_energy_index = total_charge_density_index
					
				self.intrinsic_equilibrium_fermi_energy[temperature] = intrinsic_equilibrium_fermi_energy - self.EVBM
				intrinsic_defect_hole_concentration.append(self.hole_concentrations_dict[temperature][intrinsic_equilibrium_fermi_energy_index])
				intrinsic_defect_electron_concentration.append(self.electron_concentrations_dict[temperature][intrinsic_equilibrium_fermi_energy_index])
			
			
			
			# Check if equilibrium Fermi energy can be found within band gap (for charge density including both intrinsic and extrinsic defects)
			if np.sign(total_charge_density_array[0]) == np.sign(total_charge_density_array[-1]):
				
				if self.check_outside_bandgap:
					
					total_equilibrium_fermi_energy, total_hole_concentration, total_electron_concentration = self.Check_Outside_BandGap(total_charge_density_array, temperature, calculation_type="total")
					
					self.total_equilibrium_fermi_energy[temperature] = total_equilibrium_fermi_energy - self.EVBM
					total_defect_hole_concentration.append(total_hole_concentration)
					total_defect_electron_concentration.append(total_electron_concentration)
				
				else:
					if ( (total_charge_density_array[0] < total_charge_density_array[-1]) and (np.sign(total_charge_density_array[0]) == 1) ) or ( (total_charge_density_array[0] > total_charge_density_array[-1]) and (np.sign(total_charge_density_array[0]) == -1) ):
						self.total_equilibrium_fermi_energy[temperature] = "< EVBM"
					elif ( (total_charge_density_array[0] > total_charge_density_array[-1]) and (np.sign(total_charge_density_array[0]) == 1) ) or ( (total_charge_density_array[0] < total_charge_density_array[-1]) and (np.sign(total_charge_density_array[0]) == -1) ):
						self.total_equilibrium_fermi_energy[temperature] = "> ECBM"
					total_defect_hole_concentration.append(0.0)
					total_defect_electron_concentration.append(0.0)
			
			
			
			elif np.sign(total_charge_density_array[0]) != np.sign(total_charge_density_array[-1]):
				
				for total_charge_density_index in range(len(self.fermi_energy_array)-1):
					
					# Find equilibrium Fermi energy from charge density including both intrinsic defects and user-selected extrinsic defect
					if np.sign(total_charge_density_array[total_charge_density_index]) != np.sign(total_charge_density_array[total_charge_density_index+1]):
						total_equilibrium_fermi_energy = self.fermi_energy_array[total_charge_density_index]
						total_equilibrium_fermi_energy_index = total_charge_density_index
				
				self.total_equilibrium_fermi_energy[temperature] = total_equilibrium_fermi_energy - self.EVBM
				total_defect_hole_concentration.append(self.hole_concentrations_dict[temperature][total_equilibrium_fermi_energy_index])
				total_defect_electron_concentration.append(self.electron_concentrations_dict[temperature][total_equilibrium_fermi_energy_index])
		
		return intrinsic_defect_hole_concentration, intrinsic_defect_electron_concentration, total_defect_hole_concentration, total_defect_electron_concentration
	
	
	
	def Initialize_HoleConcentration_Plot(self):
		
		intrinsic_defect_hole_concentration, intrinsic_defect_electron_concentration, total_defect_hole_concentration, total_defect_electron_concentration = self.Calculate_CarrierConcentration()
		
		try:
			self.carrier_concentration_intrinsic_defect_hole_plot.remove()
			self.carrier_concentration_total_hole_plot.remove()
		except:
			pass
		
		self.carrier_concentration_intrinsic_defect_hole_plot, = self.carrier_concentration_plot_drawing.semilogy(self.temperature_array, intrinsic_defect_hole_concentration, 'o-', color='red', label='Hole')
		if self.extrinsic_defect != "None":
			self.carrier_concentration_total_hole_plot, = self.carrier_concentration_plot_drawing.semilogy(self.temperature_array, total_defect_hole_concentration, 'o-', markerfacecolor='none', markeredgecolor='red', color='red', ls='--', label='Hole (With Dopant)')
		
		self.carrier_concentration_plot_drawing.legend(loc=1)
		self.carrier_concentration_plot_canvas.draw()
	
	
	
	def Update_HoleConcentration_Plot(self):
		
		intrinsic_defect_hole_concentration, intrinsic_defect_electron_concentration, total_defect_hole_concentration, total_defect_electron_concentration = self.Calculate_CarrierConcentration()
		self.carrier_concentration_intrinsic_defect_hole_plot.set_ydata(intrinsic_defect_hole_concentration)
		if self.extrinsic_defect != "None":
			self.carrier_concentration_total_hole_plot.set_ydata(total_defect_hole_concentration)
		
		self.carrier_concentration_plot_canvas.draw()
	
	
	
	def Initialize_ElectronConcentration_Plot(self):
		
		intrinsic_defect_hole_concentration, intrinsic_defect_electron_concentration, total_defect_hole_concentration, total_defect_electron_concentration = self.Calculate_CarrierConcentration()
		
		try:
			self.carrier_concentration_intrinsic_defect_electron_plot.remove()
			self.carrier_concentration_total_electron_plot.remove()
		except:
			pass
		
		self.carrier_concentration_intrinsic_defect_electron_plot, = self.carrier_concentration_plot_drawing.semilogy(self.temperature_array, intrinsic_defect_electron_concentration, 'o-', color='green', label='Electron')
		if self.extrinsic_defect != "None":
			self.carrier_concentration_total_electron_plot, = self.carrier_concentration_plot_drawing.semilogy(self.temperature_array, total_defect_electron_concentration, 'o-', markerfacecolor='none', markeredgecolor='green', color='green', ls='--', label='Electron (With Dopant)')
		
		self.carrier_concentration_plot_drawing.legend(loc=1)
		self.carrier_concentration_plot_canvas.draw()
	
	
	
	def Update_ElectronConcentration_Plot(self):
		
		intrinsic_defect_hole_concentration, intrinsic_defect_electron_concentration, total_defect_hole_concentration, total_defect_electron_concentration = self.Calculate_CarrierConcentration()
		self.carrier_concentration_intrinsic_defect_electron_plot.set_ydata(intrinsic_defect_electron_concentration)
		if self.extrinsic_defect != "None":
			self.carrier_concentration_total_electron_plot.set_ydata(total_defect_electron_concentration)
		
		self.carrier_concentration_plot_canvas.draw()
	
	
	
	
	
	
	
	
	
	
	def Check_Outside_BandGap(self, charge_density_array, temperature, calculation_type="intrinsic"):
		
		equilibrium_fermi_energy = 0.0
		
		# If not, determine whether it's in the valence or conduction band
		if ( (charge_density_array[0] < charge_density_array[-1]) and (np.sign(charge_density_array[0]) == 1) ) or ( (charge_density_array[0] > charge_density_array[-1]) and (np.sign(charge_density_array[0]) == -1) ):
			
			# Equilibrium Fermi energy is in the VALENCE band
			fermi_energy_increment = self.EVBM
			while equilibrium_fermi_energy == 0.0:
				
				# Fermi energy 1
				fermi_energy1 = fermi_energy_increment
				
				# Hole concentration
				fE_holes1 = self.gE_EVBM * (1. - 1./( 1. + np.exp( (self.energy_EVBM - fermi_energy1) / (self.k * temperature) ) ) )
				hole_concentration1 = integrate.simps(fE_holes1, self.energy_EVBM)
				
				# Electron concentration
				fE_electrons1 = self.gE_ECBM / (1. + np.exp( (self.energy_ECBM - fermi_energy1) / (self.k * temperature) ) )
				electron_concentration1 = integrate.simps(fE_electrons1, self.energy_ECBM)
				
				# Defects concentration (mirroring self.Calculate_Defect_Carrier_Concentration() function)
				defects_carrier_concentration1 = 0
				"""
				dmu_a = self.mu_values[self.species_a]
				dmu_b = self.mu_values[self.species_b]
				dmu_c = self.mu_values[self.species_c]
				"""
				for defect in self.ternary_defects_data.keys():
					if "_" not in defect:
						continue
					if self.ternary_defects_data[defect]["Extrinsic"] == "No":
						for charge in self.ternary_defects_data[defect]["charge"].keys():
							"""
							intrinsic_defect_formation_enthalpy = self.ternary_defects_data[defect]["charge"][charge]["Energy"] \
																- self.ternary_defects_data["supercellsize"] * self.main_compound_total_energy \
																- self.ternary_defects_data[defect]["n_"+self.species_a] * ( self.first_element_mu0 + float(dmu_a) ) \
																- self.ternary_defects_data[defect]["n_"+self.species_b] * ( self.second_element_mu0 + float(dmu_b) ) \
																- self.ternary_defects_data[defect]["n_"+self.species_c] * ( self.third_element_mu0 + float(dmu_c) ) \
																+ float(charge) * fermi_energy1 \
																+ self.ternary_defects_data[defect]["charge"][charge]["ECorr"]
							"""
							intrinsic_defect_formation_enthalpy = self.ternary_defects_data[defect]["charge"][charge]["Energy"] \
																- self.ternary_defects_data["supercellsize"] * self.main_compound_total_energy \
																+ float(charge) * fermi_energy1 \
																+ self.ternary_defects_data[defect]["charge"][charge]["ECorr"]
							for element in self.elements_list:
								intrinsic_defect_formation_enthalpy -= self.ternary_defects_data[defect]["n_"+element] * ( self.mu_elements[element]["mu0"] + self.mu_elements[element]["deltamu"] )
							if self.first_element == defect.split("_")[-1]:
								N = self.main_compound_number_first_specie/self.vol
							elif self.second_element == defect.split("_")[-1]:
								N = self.main_compound_number_second_specie/self.vol
							elif self.third_element == defect.split("_")[-1]:
								N = self.main_compound_number_third_specie/self.vol
							intrinsic_defect_charge_concentration = float(charge) * N * np.exp(-intrinsic_defect_formation_enthalpy / (self.k * temperature) )
							defects_carrier_concentration1 += intrinsic_defect_charge_concentration
				if calculation_type == "total":
					if self.extrinsic_defect != "None":
						for charge in self.ternary_defects_data[self.extrinsic_defect]["charge"].keys():
							"""
							extrinsic_defect_formation_enthalpy = self.ternary_defects_data[self.extrinsic_defect]["charge"][charge]["Energy"] \
																- self.ternary_defects_data["supercellsize"] * self.main_compound_total_energy \
																- self.ternary_defects_data[self.extrinsic_defect]["n_"+self.species_a] * ( self.first_element_mu0 + float(dmu_a) ) \
																- self.ternary_defects_data[self.extrinsic_defect]["n_"+self.species_b] * ( self.second_element_mu0 + float(dmu_b) ) \
																- self.ternary_defects_data[self.extrinsic_defect]["n_"+self.species_c] * ( self.third_element_mu0 + float(dmu_c) ) \
																- (self.extrinsic_defect_mu0 + self.extrinsic_defect_deltamu) \
																+ float(charge) * fermi_energy1 \
																+ self.ternary_defects_data[self.extrinsic_defect]["charge"][charge]["ECorr"]
							"""
							extrinsic_defect_formation_enthalpy = self.ternary_defects_data[self.extrinsic_defect]["charge"][charge]["Energy"] \
																- self.ternary_defects_data["supercellsize"] * self.main_compound_total_energy \
																- (self.extrinsic_defect_mu0 + self.extrinsic_defect_deltamu) \
																+ float(charge) * fermi_energy1 \
																+ self.ternary_defects_data[self.extrinsic_defect]["charge"][charge]["ECorr"]
							for element in self.elements_list:
								extrinsic_defect_formation_enthalpy -= self.ternary_defects_data[self.extrinsic_defect]["n_"+element] * ( self.mu_elements[element]["mu0"] + self.mu_elements[element]["deltamu"] )
							if self.first_element == defect.split("_")[-1]:
								N = self.main_compound_number_first_specie/self.vol
							elif self.second_element == defect.split("_")[-1]:
								N = self.main_compound_number_second_specie/self.vol
							elif self.third_element == defect.split("_")[-1]:
								N = self.main_compound_number_third_specie/self.vol
							extrinsic_defect_charge_concentration = float(charge) * N * np.exp(-extrinsic_defect_formation_enthalpy / (self.k * temperature) )
							defects_carrier_concentration1 += extrinsic_defect_charge_concentration
				
				# Total charge density at Fermi energy 1
				defect_charge_density1 = hole_concentration1 + electron_concentration1 + defects_carrier_concentration1
				
				# Increment Fermi energy (since we are searching the valence band, decrease the energy to which to search)
				fermi_energy_increment -= (self.ECBM-self.EVBM)/100.
				
				# Fermi energy 2 (incremented)
				fermi_energy2 = fermi_energy_increment
				
				# Hole concentration
				fE_holes2 = self.gE_EVBM * (1. - 1./( 1. + np.exp( (self.energy_EVBM - fermi_energy2) / (self.k * temperature) ) ) )
				hole_concentration2 = integrate.simps(fE_holes2, self.energy_EVBM)
				
				# Electron concentration
				fE_electrons2 = self.gE_ECBM / (1. + np.exp( (self.energy_ECBM - fermi_energy2) / (self.k * temperature) ) )
				electron_concentration2 = integrate.simps(fE_electrons2, self.energy_ECBM)
				
				# Defects concentration (mirroring self.Calculate_Defect_Carrier_Concentration() function)
				defects_carrier_concentration2 = 0
				"""
				dmu_a = self.mu_values[self.species_a]
				dmu_b = self.mu_values[self.species_b]
				dmu_c = self.mu_values[self.species_c]
				"""
				for defect in self.ternary_defects_data.keys():
					if "_" not in defect:
						continue
					if self.ternary_defects_data[defect]["Extrinsic"] == "No":
						for charge in self.ternary_defects_data[defect]["charge"].keys():
							"""
							intrinsic_defect_formation_enthalpy = self.ternary_defects_data[defect]["charge"][charge]["Energy"] \
																- self.ternary_defects_data["supercellsize"] * self.main_compound_total_energy \
																- self.ternary_defects_data[defect]["n_"+self.species_a] * ( self.first_element_mu0 + float(dmu_a) ) \
																- self.ternary_defects_data[defect]["n_"+self.species_b] * ( self.second_element_mu0 + float(dmu_b) ) \
																- self.ternary_defects_data[defect]["n_"+self.species_c] * ( self.third_element_mu0 + float(dmu_c) ) \
																+ float(charge) * fermi_energy2 \
																+ self.ternary_defects_data[defect]["charge"][charge]["ECorr"]
							"""
							intrinsic_defect_formation_enthalpy = self.ternary_defects_data[defect]["charge"][charge]["Energy"] \
																- self.ternary_defects_data["supercellsize"] * self.main_compound_total_energy \
																+ float(charge) * fermi_energy2 \
																+ self.ternary_defects_data[defect]["charge"][charge]["ECorr"]
							for element in self.elements_list:
								intrinsic_defect_formation_enthalpy -= self.ternary_defects_data[defect]["n_"+element] * ( self.mu_elements[element]["mu0"] + self.mu_elements[element]["deltamu"] )
							if self.first_element == defect.split("_")[-1]:
								N = self.main_compound_number_first_specie/self.vol
							elif self.second_element == defect.split("_")[-1]:
								N = self.main_compound_number_second_specie/self.vol
							elif self.third_element == defect.split("_")[-1]:
								N = self.main_compound_number_third_specie/self.vol
							intrinsic_defect_charge_concentration = float(charge) * N * np.exp(-intrinsic_defect_formation_enthalpy / (self.k * temperature) )
							defects_carrier_concentration2 += intrinsic_defect_charge_concentration
				if calculation_type == "total":
					if self.extrinsic_defect != "None":
						for charge in self.ternary_defects_data[self.extrinsic_defect]["charge"].keys():
							"""
							extrinsic_defect_formation_enthalpy = self.ternary_defects_data[self.extrinsic_defect]["charge"][charge]["Energy"] \
																- self.ternary_defects_data["supercellsize"] * self.main_compound_total_energy \
																- self.ternary_defects_data[self.extrinsic_defect]["n_"+self.species_a] * ( self.first_element_mu0 + float(dmu_a) ) \
																- self.ternary_defects_data[self.extrinsic_defect]["n_"+self.species_b] * ( self.second_element_mu0 + float(dmu_b) ) \
																- self.ternary_defects_data[self.extrinsic_defect]["n_"+self.species_c] * ( self.third_element_mu0 + float(dmu_c) ) \
																- (self.extrinsic_defect_mu0 + self.extrinsic_defect_deltamu) \
																+ float(charge) * fermi_energy2 \
																+ self.ternary_defects_data[self.extrinsic_defect]["charge"][charge]["ECorr"]
							"""
							extrinsic_defect_formation_enthalpy = self.ternary_defects_data[self.extrinsic_defect]["charge"][charge]["Energy"] \
																- self.ternary_defects_data["supercellsize"] * self.main_compound_total_energy \
																- (self.extrinsic_defect_mu0 + self.extrinsic_defect_deltamu) \
																+ float(charge) * fermi_energy2 \
																+ self.ternary_defects_data[self.extrinsic_defect]["charge"][charge]["ECorr"]
							for element in self.elements_list:
								extrinsic_defect_formation_enthalpy -= self.ternary_defects_data[self.extrinsic_defect]["n_"+element] * ( self.mu_elements[element]["mu0"] + self.mu_elements[element]["deltamu"] )
							if self.first_element == defect.split("_")[-1]:
								N = self.main_compound_number_first_specie/self.vol
							elif self.second_element == defect.split("_")[-1]:
								N = self.main_compound_number_second_specie/self.vol
							elif self.third_element == defect.split("_")[-1]:
								N = self.main_compound_number_third_specie/self.vol
							extrinsic_defect_charge_concentration = float(charge) * N * np.exp(-extrinsic_defect_formation_enthalpy / (self.k * temperature) )
							defects_carrier_concentration2 += extrinsic_defect_charge_concentration
				
				# Total charge density at Fermi energy 2
				defect_charge_density2 = hole_concentration2 + electron_concentration2 + defects_carrier_concentration2
				
				# Check if Fermi energy admits equilibrium conditions
				if np.sign(defect_charge_density1) != np.sign(defect_charge_density2):
					equilibrium_fermi_energy = fermi_energy1
					equilibrium_hole_concentration = hole_concentration1
					equilibrium_electron_concentration = electron_concentration1
		
		
		
		elif ( (charge_density_array[0] > charge_density_array[-1]) and (np.sign(charge_density_array[0]) == 1) ) or ( (charge_density_array[0] < charge_density_array[-1]) and (np.sign(charge_density_array[0]) == -1) ):
			
			# Equilibrium Fermi energy is in the CONDUCTION band
			fermi_energy_increment = self.ECBM
			while equilibrium_fermi_energy == 0.0:
				
				# Fermi energy 1
				fermi_energy1 = fermi_energy_increment
				
				# Hole concentration
				fE_holes1 = self.gE_EVBM * (1. - 1./( 1. + np.exp( (self.energy_EVBM - fermi_energy1) / (self.k * temperature) ) ) )
				hole_concentration1 = integrate.simps(fE_holes1, self.energy_EVBM)
				
				# Electron concentration
				fE_electrons1 = self.gE_ECBM / (1. + np.exp( (self.energy_ECBM - fermi_energy1) / (self.k * temperature) ) )
				electron_concentration1 = integrate.simps(fE_electrons1, self.energy_ECBM)
				
				# Defects concentration (mirroring self.Calculate_Defect_Carrier_Concentration() function)
				defects_carrier_concentration1 = 0
				"""
				dmu_a = self.mu_values[self.species_a]
				dmu_b = self.mu_values[self.species_b]
				dmu_c = self.mu_values[self.species_c]
				"""
				for defect in self.ternary_defects_data.keys():
					if "_" not in defect:
						continue
					if self.ternary_defects_data[defect]["Extrinsic"] == "No":
						for charge in self.ternary_defects_data[defect]["charge"].keys():
							"""
							intrinsic_defect_formation_enthalpy = self.ternary_defects_data[defect]["charge"][charge]["Energy"] \
																- self.ternary_defects_data["supercellsize"] * self.main_compound_total_energy \
																- self.ternary_defects_data[defect]["n_"+self.species_a] * ( self.first_element_mu0 + float(dmu_a) ) \
																- self.ternary_defects_data[defect]["n_"+self.species_b] * ( self.second_element_mu0 + float(dmu_b) ) \
																- self.ternary_defects_data[defect]["n_"+self.species_c] * ( self.third_element_mu0 + float(dmu_c) ) \
																+ float(charge) * fermi_energy1 \
																+ self.ternary_defects_data[defect]["charge"][charge]["ECorr"]
							"""
							intrinsic_defect_formation_enthalpy = self.ternary_defects_data[defect]["charge"][charge]["Energy"] \
																- self.ternary_defects_data["supercellsize"] * self.main_compound_total_energy \
																+ float(charge) * fermi_energy1 \
																+ self.ternary_defects_data[defect]["charge"][charge]["ECorr"]
							for element in self.elements_list:
								intrinsic_defect_formation_enthalpy -= self.ternary_defects_data[defect]["n_"+element] * ( self.mu_elements[element]["mu0"] + self.mu_elements[element]["deltamu"] )
							if self.first_element == defect.split("_")[-1]:
								N = self.main_compound_number_first_specie/self.vol
							elif self.second_element == defect.split("_")[-1]:
								N = self.main_compound_number_second_specie/self.vol
							elif self.third_element == defect.split("_")[-1]:
								N = self.main_compound_number_third_specie/self.vol
							intrinsic_defect_charge_concentration = float(charge) * N * np.exp(-intrinsic_defect_formation_enthalpy / (self.k * temperature) )
							defects_carrier_concentration1 += intrinsic_defect_charge_concentration
				if calculation_type == "total":
					if self.extrinsic_defect != "None":
						for charge in self.ternary_defects_data[self.extrinsic_defect]["charge"].keys():
							"""
							extrinsic_defect_formation_enthalpy = self.ternary_defects_data[self.extrinsic_defect]["charge"][charge]["Energy"] \
																- self.ternary_defects_data["supercellsize"] * self.main_compound_total_energy \
																- self.ternary_defects_data[self.extrinsic_defect]["n_"+self.species_a] * ( self.first_element_mu0 + float(dmu_a) ) \
																- self.ternary_defects_data[self.extrinsic_defect]["n_"+self.species_b] * ( self.second_element_mu0 + float(dmu_b) ) \
																- self.ternary_defects_data[self.extrinsic_defect]["n_"+self.species_c] * ( self.third_element_mu0 + float(dmu_c) ) \
																- (self.extrinsic_defect_mu0 + self.extrinsic_defect_deltamu) \
																+ float(charge) * fermi_energy1 \
																+ self.ternary_defects_data[self.extrinsic_defect]["charge"][charge]["ECorr"]
							"""
							extrinsic_defect_formation_enthalpy = self.ternary_defects_data[self.extrinsic_defect]["charge"][charge]["Energy"] \
																- self.ternary_defects_data["supercellsize"] * self.main_compound_total_energy \
																- (self.extrinsic_defect_mu0 + self.extrinsic_defect_deltamu) \
																+ float(charge) * fermi_energy1 \
																+ self.ternary_defects_data[self.extrinsic_defect]["charge"][charge]["ECorr"]
							for element in self.elements_list:
								extrinsic_defect_formation_enthalpy -= self.ternary_defects_data[self.extrinsic_defect]["n_"+element] * ( self.mu_elements[element]["mu0"] + self.mu_elements[element]["deltamu"] )
							if self.first_element == defect.split("_")[-1]:
								N = self.main_compound_number_first_specie/self.vol
							elif self.second_element == defect.split("_")[-1]:
								N = self.main_compound_number_second_specie/self.vol
							elif self.third_element == defect.split("_")[-1]:
								N = self.main_compound_number_third_specie/self.vol
							extrinsic_defect_charge_concentration = float(charge) * N * np.exp(-extrinsic_defect_formation_enthalpy / (self.k * temperature) )
							defects_carrier_concentration1 += extrinsic_defect_charge_concentration
				
				# Total charge density at Fermi energy 1
				defect_charge_density1 = hole_concentration1 + electron_concentration1 + defects_carrier_concentration1
				
				# Increment Fermi energy (since we are searching the conduction band, inrease the energy to which to search)
				fermi_energy_increment += (self.ECBM-self.EVBM)/100.
				
				# Fermi energy 2 (incremented)
				fermi_energy2 = fermi_energy_increment
				
				# Hole concentration
				fE_holes2 = self.gE_EVBM * (1. - 1./( 1. + np.exp( (self.energy_EVBM - fermi_energy2) / (self.k * temperature) ) ) )
				hole_concentration2 = integrate.simps(fE_holes2, self.energy_EVBM)
				
				# Electron concentration
				fE_electrons2 = self.gE_ECBM / (1. + np.exp( (self.energy_ECBM - fermi_energy2) / (self.k * temperature) ) )
				electron_concentration2 = integrate.simps(fE_electrons2, self.energy_ECBM)
				
				# Defects concentration (mirroring self.Calculate_Defect_Carrier_Concentration() function)
				defects_carrier_concentration2 = 0
				"""
				dmu_a = self.mu_values[self.species_a]
				dmu_b = self.mu_values[self.species_b]
				dmu_c = self.mu_values[self.species_c]
				"""
				for defect in self.ternary_defects_data.keys():
					if "_" not in defect:
						continue
					if self.ternary_defects_data[defect]["Extrinsic"] == "No":
						for charge in self.ternary_defects_data[defect]["charge"].keys():
							"""
							intrinsic_defect_formation_enthalpy = self.ternary_defects_data[defect]["charge"][charge]["Energy"] \
																- self.ternary_defects_data["supercellsize"] * self.main_compound_total_energy \
																- self.ternary_defects_data[defect]["n_"+self.species_a] * ( self.first_element_mu0 + float(dmu_a) ) \
																- self.ternary_defects_data[defect]["n_"+self.species_b] * ( self.second_element_mu0 + float(dmu_b) ) \
																- self.ternary_defects_data[defect]["n_"+self.species_c] * ( self.third_element_mu0 + float(dmu_c) ) \
																+ float(charge) * fermi_energy2 \
																+ self.ternary_defects_data[defect]["charge"][charge]["ECorr"]
							"""
							intrinsic_defect_formation_enthalpy = self.ternary_defects_data[defect]["charge"][charge]["Energy"] \
																- self.ternary_defects_data["supercellsize"] * self.main_compound_total_energy \
																+ float(charge) * fermi_energy2 \
																+ self.ternary_defects_data[defect]["charge"][charge]["ECorr"]
							for element in self.elements_list:
								intrinsic_defect_formation_enthalpy -= self.ternary_defects_data[defect]["n_"+element] * ( self.mu_elements[element]["mu0"] + self.mu_elements[element]["deltamu"] )
							if self.first_element == defect.split("_")[-1]:
								N = self.main_compound_number_first_specie/self.vol
							elif self.second_element == defect.split("_")[-1]:
								N = self.main_compound_number_second_specie/self.vol
							elif self.third_element == defect.split("_")[-1]:
								N = self.main_compound_number_third_specie/self.vol
							intrinsic_defect_charge_concentration = float(charge) * N * np.exp(-intrinsic_defect_formation_enthalpy / (self.k * temperature) )
							defects_carrier_concentration2 += intrinsic_defect_charge_concentration
				if calculation_type == "total":
					if self.extrinsic_defect != "None":
						for charge in self.ternary_defects_data[self.extrinsic_defect]["charge"].keys():
							"""
							extrinsic_defect_formation_enthalpy = self.ternary_defects_data[self.extrinsic_defect]["charge"][charge]["Energy"] \
																- self.ternary_defects_data["supercellsize"] * self.main_compound_total_energy \
																- self.ternary_defects_data[self.extrinsic_defect]["n_"+self.species_a] * (self.first_element_mu0 + float(dmu_a) ) \
																- self.ternary_defects_data[self.extrinsic_defect]["n_"+self.species_b] * ( self.second_element_mu0 + float(dmu_b) ) \
																- self.ternary_defects_data[self.extrinsic_defect]["n_"+self.species_c] * ( self.third_element_mu0 + float(dmu_c) ) \
																- (self.extrinsic_defect_mu0 + self.extrinsic_defect_deltamu) \
																+ float(charge) * fermi_energy2 \
																+ self.ternary_defects_data[self.extrinsic_defect]["charge"][charge]["ECorr"]
							"""
							extrinsic_defect_formation_enthalpy = self.ternary_defects_data[self.extrinsic_defect]["charge"][charge]["Energy"] \
																- self.ternary_defects_data["supercellsize"] * self.main_compound_total_energy \
																- (self.extrinsic_defect_mu0 + self.extrinsic_defect_deltamu) \
																+ float(charge) * fermi_energy2 \
																+ self.ternary_defects_data[self.extrinsic_defect]["charge"][charge]["ECorr"]
							for element in self.elements_list:
								extrinsic_defect_formation_enthalpy -= self.ternary_defects_data[self.extrinsic_defect]["n_"+element] * (self.mu_elements[element]["mu0"] + self.mu_elements[element]["deltamu"] )
							if self.first_element == defect.split("_")[-1]:
								N = self.main_compound_number_first_specie/self.vol
							elif self.second_element == defect.split("_")[-1]:
								N = self.main_compound_number_second_specie/self.vol
							elif self.third_element == defect.split("_")[-1]:
								N = self.main_compound_number_third_specie/self.vol
							extrinsic_defect_charge_concentration = float(charge) * N * np.exp(-extrinsic_defect_formation_enthalpy / (self.k * temperature) )
							defects_carrier_concentration2 += extrinsic_defect_charge_concentration
				
				# Total charge density at Fermi energy 2
				defect_charge_density2 = hole_concentration2 + electron_concentration2 + defects_carrier_concentration2
				
				# Check if Fermi energy admits equilibrium conditions
				if np.sign(defect_charge_density1) != np.sign(defect_charge_density2):
					equilibrium_fermi_energy = fermi_energy1
					equilibrium_hole_concentration = hole_concentration1
					equilibrium_electron_concentration = electron_concentration1
		
		return equilibrium_fermi_energy, equilibrium_hole_concentration, equilibrium_electron_concentration













