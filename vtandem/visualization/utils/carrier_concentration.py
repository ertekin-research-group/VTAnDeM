
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'
__author__ = 'Michael_Lidia_Jiaxing_Elif'


import numpy as np
import os
from scipy import integrate

from vtandem.visualization.utils.defect_formation_energy import *



def Calculate_FreeHole_FreeElectron_Concentrations(	temperature_array, \
													fermi_energy_array, \
													gE_ValenceBand, \
													energies_ValenceBand, \
													gE_ConductionBand, \
													energies_ConductionBand ):
	
	k = 8.6173303E-5

	hole_concentrations_dict = {}
	electron_concentrations_dict = {}

	for temperature in temperature_array:
		
		hole_concentrations_dict[temperature] = []
		electron_concentrations_dict[temperature] = []
		
		for ef in fermi_energy_array:
			
			# Hole concentration
			fE_holes = gE_ValenceBand * (1. - 1./( 1. + np.exp( (energies_ValenceBand - ef) / (k * temperature) ) ) )
			hole_concentration = integrate.simps(fE_holes, energies_ValenceBand)
			hole_concentrations_dict[temperature].append(hole_concentration)	# In units of cm^-3
			
			# Electron concentration
			fE_electrons = gE_ConductionBand / (1. + np.exp( (energies_ConductionBand - ef) / (k * temperature) ) )
			electron_concentration = integrate.simps(fE_electrons, energies_ConductionBand)
			electron_concentrations_dict[temperature].append(electron_concentration)	# In units of cm^-3
		
		hole_concentrations_dict[temperature] = np.asarray(hole_concentrations_dict[temperature])
		electron_concentrations_dict[temperature] = np.asarray(electron_concentrations_dict[temperature])

	return hole_concentrations_dict, electron_concentrations_dict



def Calculate_Defect_Carrier_Concentration(	defects_data, \
											main_compound_info, \
											mu_elements, \
											temperature_array, \
											fermi_energy_array, \
											volume, \
											extrinsic_defects, \
											dopant, \
											dopant_mu0, \
											dopant_deltamu, \
											synthesis_temperature = None ):
	
	k = 8.6173303E-5
	
	# Initialize
	intrinsic_defect_carrier_concentration_temperature = {}
	extrinsic_defect_carrier_concentration_temperature = {}
	for temperature in temperature_array:
		intrinsic_defect_carrier_concentration_temperature[temperature] = np.zeros(len(fermi_energy_array))
		extrinsic_defect_carrier_concentration_temperature[temperature] = np.zeros(len(fermi_energy_array))

	# Obtain formation enthalpies of intrinsic defects
	intrinsic_defects_enthalpy_data = Calculate_IntrinsicDefectFormationEnthalpies(	defects_data, \
																					main_compound_info, \
																					fermi_energy_array, \
																					mu_elements )
	
	# Calculate intrinsic defect carrier concentration
	# Loop through intrinsic defects
	for intrinsic_defect in intrinsic_defects_enthalpy_data.keys():
		
		# Loop through charge states
		for charge in intrinsic_defects_enthalpy_data[intrinsic_defect].keys():
			
			# Prefactor
			N = defects_data[intrinsic_defect]["site_multiplicity"] / volume
			
			# Defect concentration
			for temperature in temperature_array:
				if synthesis_temperature is None:
					defect_carrier_concentration = float(charge) * N * np.exp( -intrinsic_defects_enthalpy_data[intrinsic_defect][charge] / (k * temperature) )
				elif synthesis_temperature is not None:
					defect_carrier_concentration = float(charge) * N * np.exp( -intrinsic_defects_enthalpy_data[intrinsic_defect][charge] / (k * synthesis_temperature) )
				
				"""
				if (intrinsic_defect == "V_Se") and (charge == "+2") and (temperature == 1000):
					print(intrinsic_defect, charge, temperature, defect_carrier_concentration, intrinsic_defects_enthalpy_data[intrinsic_defect][charge], defects_data[intrinsic_defect]["site_multiplicity"], volume)
				"""
				
				intrinsic_defect_carrier_concentration_temperature[temperature] += defect_carrier_concentration
	
	# Check if the user-selected dopant is "None"
	if dopant == "None":
		return intrinsic_defect_carrier_concentration_temperature, extrinsic_defect_carrier_concentration_temperature
	
	# Obtain the formation enthalpy of dopant on different sites (i.e. extrinsic defects)
	extrinsic_defects_enthalpy_data = Calculate_ExtrinsicDefectFormationEnthalpies(	defects_data, \
																					main_compound_info, \
																					fermi_energy_array, \
																					mu_elements, \
																					extrinsic_defects, \
																					dopant, \
																					dopant_mu0, \
																					dopant_deltamu )
	
	print(dopant, extrinsic_defects)

	# Loop through extrinsic defects
	for extrinsic_defect in extrinsic_defects:

		# Carrier concentration prefactor
		N_extrinsic = defects_data[extrinsic_defect]["site_multiplicity"] / volume
		
		# Calculate extrinsic defect carrier concentration
		# Loop through charge states
		for charge in extrinsic_defects_enthalpy_data[extrinsic_defect].keys():
			
			# Defect concentration
			for temperature in temperature_array:
				if synthesis_temperature is None:
					extrinsic_defect_carrier_concentration = float(charge) * N_extrinsic * np.exp( -extrinsic_defects_enthalpy_data[extrinsic_defect][charge] / (k * temperature) )
				elif synthesis_temperature is not None:
					extrinsic_defect_carrier_concentration = float(charge) * N_extrinsic * np.exp( -extrinsic_defects_enthalpy_data[extrinsic_defect][charge] / (k * synthesis_temperature) )
				extrinsic_defect_carrier_concentration_temperature[temperature] += extrinsic_defect_carrier_concentration
		
	# Returns a dictionary with temperature as keys and array of defect-induced carrier concentrations (not defect concentrations) as values
	return intrinsic_defect_carrier_concentration_temperature, extrinsic_defect_carrier_concentration_temperature



def Calculate_CarrierConcentration(	EVBM, \
									ECBM, \
									energies_ValenceBand, \
									gE_ValenceBand, \
									energies_ConductionBand, \
									gE_ConductionBand, \
									defects_data, \
									main_compound_info, \
									mu_elements, \
									temperature_array, \
									fermi_energy_array, \
									volume, \
									extrinsic_defects, \
									dopant, \
									dopant_mu0, \
									dopant_deltamu, \
									hole_concentrations_dict, \
									electron_concentrations_dict, \
									synthesis_temperature = None ):
	
	# Calculate defect carrier concentration (for intrinsic defects and extrinsic defects)
	intrinsic_defect_carrier_concentration_temperature, extrinsic_defect_carrier_concentration_temperature = Calculate_Defect_Carrier_Concentration(	defects_data = defects_data, \
																																						main_compound_info = main_compound_info, \
																																						mu_elements = mu_elements, \
																																						temperature_array = temperature_array, \
																																						fermi_energy_array = fermi_energy_array, \
																																						volume = volume, \
																																						extrinsic_defects = extrinsic_defects, \
																																						dopant = dopant, \
																																						dopant_mu0 = dopant_mu0, \
																																						dopant_deltamu = dopant_deltamu, \
																																						synthesis_temperature = synthesis_temperature )
	
	# Carrier concentrations from intrinsic defects only
	intrinsic_defect_hole_concentration = []
	intrinsic_defect_electron_concentration = []
	
	# Carrier concentration from both intrinsic and extrinsic defects
	total_hole_concentration = []
	total_electron_concentration = []
	
	# Track equilibrium Fermi energy for both intrinsic defects only and total at each temperature
	intrinsic_equilibrium_fermi_energy_temperature = {}
	total_equilibrium_fermi_energy_temperature = {}
	
	for temperature in temperature_array:
		
		# Charge density including intrinsic defects only
		intrinsic_defect_charge_density_array = intrinsic_defect_carrier_concentration_temperature[temperature] + hole_concentrations_dict[temperature] - electron_concentrations_dict[temperature]
		intrinsic_equilibrium_fermi_energy = 0.0
		intrinsic_equilibrium_fermi_energy_index = 0
		
		# Charge density including both intrinsic and extrinsic defects
		total_charge_density_array = intrinsic_defect_carrier_concentration_temperature[temperature] + extrinsic_defect_carrier_concentration_temperature[temperature] + hole_concentrations_dict[temperature] - electron_concentrations_dict[temperature]
		total_equilibrium_fermi_energy = 0.0
		total_equilibrium_fermi_energy_index = 0
		

		# Search for equilibrium Fermi energy within band gap of material (for only intrinsic defects)
		for intrinsic_defect_charge_density_index in range(len(fermi_energy_array)-1):
			
			# Find equilibrium Fermi energy from charge density including intrinsic defects only
			if np.sign(intrinsic_defect_charge_density_array[intrinsic_defect_charge_density_index]) != np.sign(intrinsic_defect_charge_density_array[intrinsic_defect_charge_density_index+1]):
				intrinsic_equilibrium_fermi_energy = fermi_energy_array[intrinsic_defect_charge_density_index]
				intrinsic_equilibrium_fermi_energy_index = intrinsic_defect_charge_density_index
				break
		
		intrinsic_equilibrium_fermi_energy_temperature[temperature] = intrinsic_equilibrium_fermi_energy - EVBM
		intrinsic_defect_hole_concentration.append(hole_concentrations_dict[temperature][intrinsic_equilibrium_fermi_energy_index])
		intrinsic_defect_electron_concentration.append(electron_concentrations_dict[temperature][intrinsic_equilibrium_fermi_energy_index])


		
		# Search for equilibrium Fermi energy within band gap of material
		for total_charge_density_index in range(len(fermi_energy_array)-1):
			
			# Find equilibrium Fermi energy from charge density including both intrinsic defects and user-selected extrinsic defect
			if np.sign(total_charge_density_array[total_charge_density_index]) != np.sign(total_charge_density_array[total_charge_density_index+1]):
				total_equilibrium_fermi_energy = fermi_energy_array[total_charge_density_index]
				total_equilibrium_fermi_energy_index = total_charge_density_index
				break
		
		total_equilibrium_fermi_energy_temperature[temperature] = total_equilibrium_fermi_energy - EVBM
		total_hole_concentration.append(hole_concentrations_dict[temperature][total_equilibrium_fermi_energy_index])
		total_electron_concentration.append(electron_concentrations_dict[temperature][total_equilibrium_fermi_energy_index])




	return intrinsic_defect_hole_concentration, intrinsic_defect_electron_concentration, total_hole_concentration, total_electron_concentration, intrinsic_equilibrium_fermi_energy_temperature, total_equilibrium_fermi_energy_temperature



