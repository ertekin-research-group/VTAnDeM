
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'
__author__ = 'Michael_Lidia_Jiaxing_Elif'


import numpy as np
import os
from scipy import integrate

from vtandem.visualization.defect_formation_energy import Calculate_IntrinsicDefectFormationEnthalpies
from vtandem.visualization.defect_formation_energy import Calculate_ExtrinsicDefectFormationEnthalpies
from vtandem.visualization.defect_formation_energy import Find_SiteMultiplicity



def Calculate_Defect_Carrier_Concentration(	defects_data, \
											main_compound_total_energy, \
											mu_elements, \
											temperature_array, \
											fermi_energy_array, \
											volume, \
											number_species, \
											extrinsic_defect, \
											extrinsic_defect_mu0, \
											extrinsic_defect_deltamu, \
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
																					main_compound_total_energy, \
																					fermi_energy_array, \
																					mu_elements )
	
	# Calculate intrinsic defect carrier concentration
	# Loop through intrinsic defects
	for intrinsic_defect in intrinsic_defects_enthalpy_data.keys():
		
		# Prefactor
		N = Find_SiteMultiplicity(intrinsic_defect, len(number_species), number_species, volume)
		
		# Loop through charge states
		for charge in intrinsic_defects_enthalpy_data[intrinsic_defect].keys():
			
			# Defect concentration
			for temperature in temperature_array:
				if synthesis_temperature is None:
					defect_carrier_concentration = float(charge) * N * np.exp( -intrinsic_defects_enthalpy_data[intrinsic_defect][charge] / (k * temperature) )
				elif synthesis_temperature is not None:
					defect_carrier_concentration = float(charge) * N * np.exp( -intrinsic_defects_enthalpy_data[intrinsic_defect][charge] / (k * synthesis_temperature) )
				intrinsic_defect_carrier_concentration_temperature[temperature] += defect_carrier_concentration
	
	# Check that user-selected extrinsic defect is not "None"
	if extrinsic_defect == "None":
		return intrinsic_defect_carrier_concentration_temperature, extrinsic_defect_carrier_concentration_temperature
	
	# Obtain the formation enthalpy of extrinsic defect
	extrinsic_defects_enthalpy_data = Calculate_ExtrinsicDefectFormationEnthalpies(	defects_data, \
																					main_compound_total_energy, \
																					fermi_energy_array, \
																					mu_elements, \
																					extrinsic_defect, \
																					extrinsic_defect_mu0, \
																					extrinsic_defect_deltamu )
	
	# Carrier concentration prefactor
	N_extrinsic = Find_SiteMultiplicity(extrinsic_defect, len(number_species), number_species, volume)
	
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
									main_compound_total_energy, \
									mu_elements, \
									temperature_array, \
									fermi_energy_array, \
									volume, \
									number_species, \
									extrinsic_defect, \
									extrinsic_defect_mu0, \
									extrinsic_defect_deltamu, \
									hole_concentrations_dict, \
									electron_concentrations_dict, \
									check_outside_bandgap, \
									synthesis_temperature = None ):
	
	# Carrier concentrations from intrinsic defects only
	intrinsic_defect_hole_concentration = []
	intrinsic_defect_electron_concentration = []
	
	# Carrier concentration from both intrinsic and extrinsic defects
	total_hole_concentration = []
	total_electron_concentration = []
	
	# Track equilibrium Fermi energy for both intrinsic defects only and total at each temperature
	intrinsic_equilibrium_fermi_energy_temperature = {}
	total_equilibrium_fermi_energy_temperature = {}
	
	# Calculate defect carrier concentration (for intrinsic defects and extrinsic defects)
	intrinsic_defect_carrier_concentration_temperature, extrinsic_defect_carrier_concentration_temperature = Calculate_Defect_Carrier_Concentration(	defects_data, \
																																						main_compound_total_energy, \
																																						mu_elements, \
																																						temperature_array, \
																																						fermi_energy_array, \
																																						volume, \
																																						number_species, \
																																						extrinsic_defect, \
																																						extrinsic_defect_mu0, \
																																						extrinsic_defect_deltamu, \
																																						synthesis_temperature )
	
	for temperature in temperature_array:
		
		# Charge density including intrinsic defects only
		intrinsic_defect_charge_density_array = intrinsic_defect_carrier_concentration_temperature[temperature] + hole_concentrations_dict[temperature] - electron_concentrations_dict[temperature]
		intrinsic_equilibrium_fermi_energy = 0.0
		intrinsic_equilibrium_fermi_energy_index = 0
		
		# Charge density including both intrinsic and extrinsic defects
		total_charge_density_array = intrinsic_defect_carrier_concentration_temperature[temperature] + extrinsic_defect_carrier_concentration_temperature[temperature] + hole_concentrations_dict[temperature] - electron_concentrations_dict[temperature]
		total_equilibrium_fermi_energy = 0.0
		total_equilibrium_fermi_energy_index = 0
		
		# Check if equilibrium Fermi energy can be found within band gap (for charge density including INTRINSIC DEFECTS ONLY)
		if np.sign(intrinsic_defect_charge_density_array[0]) == np.sign(intrinsic_defect_charge_density_array[-1]):
			
			if check_outside_bandgap:
				# Find equilibrium Fermi energy and carrier concentration when EF is outside the band gap
				intrinsic_equilibrium_fermi_energy, intrinsic_hole_concentration, intrinsic_electron_concentration = Check_Outside_BandGap(	intrinsic_defect_charge_density_array, \
																																			EVBM, \
																																			ECBM, \
																																			energies_ValenceBand, \
																																			gE_ValenceBand, \
																																			energies_ConductionBand, \
																																			gE_ConductionBand, \
																																			temperature, \
																																			defects_data, \
																																			main_compound_total_energy, \
																																			mu_elements, \
																																			volume, \
																																			number_species, \
																																			extrinsic_defect, \
																																			extrinsic_defect_mu0, \
																																			extrinsic_defect_deltamu, \
																																			calculation_type = "intrinsic", \
																																			synthesis_temperature = synthesis_temperature	)
				intrinsic_equilibrium_fermi_energy_temperature[temperature] = intrinsic_equilibrium_fermi_energy - EVBM
				intrinsic_defect_hole_concentration.append(intrinsic_hole_concentration)
				intrinsic_defect_electron_concentration.append(intrinsic_electron_concentration)
			
			else:
				# Don't search for equilibrium Fermi energy, record carrier concentrations as 0.0
				if ( (intrinsic_defect_charge_density_array[0] < intrinsic_defect_charge_density_array[-1]) and (np.sign(intrinsic_defect_charge_density_array[0]) == 1) ) or ( (intrinsic_defect_charge_density_array[0] > intrinsic_defect_charge_density_array[-1]) and (np.sign(intrinsic_defect_charge_density_array[0]) == -1) ):
					intrinsic_equilibrium_fermi_energy_temperature[temperature] = "< EVBM"
				elif ( (intrinsic_defect_charge_density_array[0] > intrinsic_defect_charge_density_array[-1]) and (np.sign(intrinsic_defect_charge_density_array[0]) == 1) ) or ( (intrinsic_defect_charge_density_array[0] < intrinsic_defect_charge_density_array[-1]) and (np.sign(intrinsic_defect_charge_density_array[0]) == -1) ):
					intrinsic_equilibrium_fermi_energy_temperature[temperature] = "> ECBM"
				intrinsic_defect_hole_concentration.append(0.0)
				intrinsic_defect_electron_concentration.append(0.0)
		
		elif np.sign(intrinsic_defect_charge_density_array[0]) != np.sign(intrinsic_defect_charge_density_array[-1]):
			
			# Search for equilibrium Fermi energy within band gap of material
			for intrinsic_defect_charge_density_index in range(len(fermi_energy_array)-1):
				
				# Find equilibrium Fermi energy from charge density including intrinsic defects only
				if np.sign(intrinsic_defect_charge_density_array[intrinsic_defect_charge_density_index]) != np.sign(intrinsic_defect_charge_density_array[intrinsic_defect_charge_density_index+1]):
					intrinsic_equilibrium_fermi_energy = fermi_energy_array[intrinsic_defect_charge_density_index]
					intrinsic_equilibrium_fermi_energy_index = intrinsic_defect_charge_density_index
			
			intrinsic_equilibrium_fermi_energy_temperature[temperature] = intrinsic_equilibrium_fermi_energy - EVBM
			intrinsic_defect_hole_concentration.append(hole_concentrations_dict[temperature][intrinsic_equilibrium_fermi_energy_index])
			intrinsic_defect_electron_concentration.append(electron_concentrations_dict[temperature][intrinsic_equilibrium_fermi_energy_index])
		
		
		# Check if equilibrium Fermi energy can be found within band gap (for charge density including BOTH intrinsic and extrinsic defects)
		if np.sign(total_charge_density_array[0]) == np.sign(total_charge_density_array[-1]):
			
			if check_outside_bandgap:
				# Find equilibrium Fermi energy and carrier concentration when EF is outside the band gap
				total_equilibrium_fermi_energy, total_equilibrium_hole_concentration, total_equilibrium_electron_concentration = Check_Outside_BandGap(	total_charge_density_array, \
																																						EVBM, \
																																						ECBM, \
																																						energies_ValenceBand, \
																																						gE_ValenceBand, \
																																						energies_ConductionBand, \
																																						gE_ConductionBand, \
																																						temperature, \
																																						defects_data, \
																																						main_compound_total_energy, \
																																						mu_elements, \
																																						volume, \
																																						number_species, \
																																						extrinsic_defect, \
																																						extrinsic_defect_mu0, \
																																						extrinsic_defect_deltamu, \
																																						calculation_type = "total", \
																																						synthesis_temperature = synthesis_temperature	)
				total_equilibrium_fermi_energy_temperature[temperature] = total_equilibrium_fermi_energy - EVBM
				total_hole_concentration.append(total_equilibrium_hole_concentration)
				total_electron_concentration.append(total_equilibrium_electron_concentration)
			
			else:
				# Don't search for equilibrium Fermi energy, record carrier concentrations as 0.0
				if ( (total_charge_density_array[0] < total_charge_density_array[-1]) and (np.sign(total_charge_density_array[0]) == 1) ) or ( (total_charge_density_array[0] > total_charge_density_array[-1]) and (np.sign(total_charge_density_array[0]) == -1) ):
					total_equilibrium_fermi_energy_temperature[temperature] = "< EVBM"
				elif ( (total_charge_density_array[0] > total_charge_density_array[-1]) and (np.sign(total_charge_density_array[0]) == 1) ) or ( (total_charge_density_array[0] < total_charge_density_array[-1]) and (np.sign(total_charge_density_array[0]) == -1) ):
					total_equilibrium_fermi_energy_temperature[temperature] = "> ECBM"
				total_hole_concentration.append(0.0)
				total_electron_concentration.append(0.0)
		
		elif np.sign(total_charge_density_array[0]) != np.sign(total_charge_density_array[-1]):
			
			# Search for equilibrium Fermi energy within band gap of material
			for total_charge_density_index in range(len(fermi_energy_array)-1):
				
				# Find equilibrium Fermi energy from charge density including both intrinsic defects and user-selected extrinsic defect
				if np.sign(total_charge_density_array[total_charge_density_index]) != np.sign(total_charge_density_array[total_charge_density_index+1]):
					total_equilibrium_fermi_energy = fermi_energy_array[total_charge_density_index]
					total_equilibrium_fermi_energy_index = total_charge_density_index
			
			total_equilibrium_fermi_energy_temperature[temperature] = total_equilibrium_fermi_energy - EVBM
			total_hole_concentration.append(hole_concentrations_dict[temperature][total_equilibrium_fermi_energy_index])
			total_electron_concentration.append(electron_concentrations_dict[temperature][total_equilibrium_fermi_energy_index])
	
	return intrinsic_defect_hole_concentration, intrinsic_defect_electron_concentration, total_hole_concentration, total_electron_concentration, intrinsic_equilibrium_fermi_energy_temperature, total_equilibrium_fermi_energy_temperature



def Check_Outside_BandGap( 	charge_density_array, \
							EVBM, \
							ECBM, \
							energies_ValenceBand, \
							gE_ValenceBand, \
							energies_ConductionBand, \
							gE_ConductionBand, \
							temperature, \
							defects_data, \
							main_compound_total_energy, \
							mu_elements, \
							volume, \
							number_species, \
							extrinsic_defect, \
							extrinsic_defect_mu0, \
							extrinsic_defect_deltamu, \
							calculation_type="intrinsic", \
							synthesis_temperature = None	):
	
	k = 8.6173303E-5
	equilibrium_fermi_energy = 0.0
	
	# If not, determine whether it's in the valence or conduction band
	if ( (charge_density_array[0] < charge_density_array[-1]) and (np.sign(charge_density_array[0]) == 1) ) or ( (charge_density_array[0] > charge_density_array[-1]) and (np.sign(charge_density_array[0]) == -1) ):
		
		# Equilibrium Fermi energy is in the VALENCE band
		fermi_energy_increment = EVBM
		while equilibrium_fermi_energy == 0.0:
			
			# Fermi energy 1
			fermi_energy1 = fermi_energy_increment
			
			# Hole concentration
			fE_holes1 = gE_ValenceBand * (1. - 1./( 1. + np.exp( (energies_ValenceBand - fermi_energy1) / (k * temperature) ) ) )
			hole_concentration1 = integrate.simps(fE_holes1, energies_ValenceBand)
			
			# Electron concentration
			fE_electrons1 = gE_ConductionBand / (1. + np.exp( (energies_ConductionBand - fermi_energy1) / (k * temperature) ) )
			electron_concentration1 = integrate.simps(fE_electrons1, energies_ConductionBand)
			
			# Defects concentration (mirroring Calculate_Defect_Carrier_Concentration() function)
			defects_carrier_concentration1 = 0
			
			# Calculate intrinsic defect carrier concentration
			intrinsic_defects_enthalpy_data = Calculate_IntrinsicDefectFormationEnthalpies(	defects_data, \
																							main_compound_total_energy, \
																							fermi_energy1, \
																							mu_elements )
			for intrinsic_defect in intrinsic_defects_enthalpy_data.keys():
				N = Find_SiteMultiplicity(intrinsic_defect, len(number_species), number_species, volume)
				for charge in intrinsic_defects_enthalpy_data[intrinsic_defect].keys():
					if synthesis_temperature is None:
						intrinsic_defect_charge_concentration = float(charge) * N * np.exp( -intrinsic_defects_enthalpy_data[intrinsic_defect][charge] / (k * temperature) )
					elif synthesis_temperature is not None:
						intrinsic_defect_charge_concentration = float(charge) * N * np.exp( -intrinsic_defects_enthalpy_data[intrinsic_defect][charge] / (k * synthesis_temperature) )
					defects_carrier_concentration1 += intrinsic_defect_charge_concentration
			
			if calculation_type == "total":
				if extrinsic_defect != "None":
					
					# Obtain extrinsic defect formation enthalpy
					extrinsic_defects_enthalpy_data = Calculate_ExtrinsicDefectFormationEnthalpies(	defects_data, \
																									main_compound_total_energy, \
																									fermi_energy1, \
																									mu_elements, \
																									extrinsic_defect, \
																									extrinsic_defect_mu0, \
																									extrinsic_defect_deltamu )
					
					# Carrier concentration prefactor
					N = Find_SiteMultiplicity(extrinsic_defect, len(number_species), number_species, volume)
					
					# Defect concentration
					for charge in extrinsic_defects_enthalpy_data[extrinsic_defect].keys():
						if synthesis_temperature is None:
							extrinsic_defect_charge_concentration = float(charge) * N * np.exp( -extrinsic_defects_enthalpy_data[extrinsic_defect][charge] / (k * temperature) )
						elif synthesis_temperature is not None:
							extrinsic_defect_charge_concentration = float(charge) * N * np.exp( -extrinsic_defects_enthalpy_data[extrinsic_defect][charge] / (k * synthesis_temperature) )
						defects_carrier_concentration1 += extrinsic_defect_charge_concentration
			
			# Total charge density at Fermi energy 1
			defect_charge_density1 = hole_concentration1 + electron_concentration1 + defects_carrier_concentration1
			
			# Increment Fermi energy (since we are searching the valence band, decrease the energy to which to search)
			fermi_energy_increment -= (ECBM-EVBM)/100.
			
			# Fermi energy 2 (incremented)
			fermi_energy2 = fermi_energy_increment
			
			# Hole concentration
			fE_holes2 = gE_ValenceBand * (1. - 1./( 1. + np.exp( (energies_ValenceBand - fermi_energy2) / (k * temperature) ) ) )
			hole_concentration2 = integrate.simps(fE_holes2, energies_ValenceBand)
			
			# Electron concentration
			fE_electrons2 = gE_ConductionBand / (1. + np.exp( (energies_ConductionBand - fermi_energy2) / (k * temperature) ) )
			electron_concentration2 = integrate.simps(fE_electrons2, energies_ConductionBand)
			
			# Defects concentration (mirroring Calculate_Defect_Carrier_Concentration() function, with new value of Fermi energy)
			defects_carrier_concentration2 = 0
			intrinsic_defects_enthalpy_data = Calculate_IntrinsicDefectFormationEnthalpies(	defects_data, \
																							main_compound_total_energy, \
																							fermi_energy2, \
																							mu_elements )
			for intrinsic_defect in intrinsic_defects_enthalpy_data.keys():
				N = Find_SiteMultiplicity(intrinsic_defect, len(number_species), number_species, volume)
				for charge in intrinsic_defects_enthalpy_data[intrinsic_defect].keys():
					if synthesis_temperature is None:
						intrinsic_defect_charge_concentration = float(charge) * N * np.exp( -intrinsic_defects_enthalpy_data[intrinsic_defect][charge] / (k * temperature) )
					elif synthesis_temperature is not None:
						intrinsic_defect_charge_concentration = float(charge) * N * np.exp( -intrinsic_defects_enthalpy_data[intrinsic_defect][charge] / (k * synthesis_temperature) )
					defects_carrier_concentration2 += intrinsic_defect_charge_concentration
			if calculation_type == "total":
				if extrinsic_defect != "None":
					extrinsic_defects_enthalpy_data = Calculate_ExtrinsicDefectFormationEnthalpies(	defects_data, \
																									main_compound_total_energy, \
																									fermi_energy2, \
																									mu_elements, \
																									extrinsic_defect, \
																									extrinsic_defect_mu0, \
																									extrinsic_defect_deltamu )
					N = Find_SiteMultiplicity(extrinsic_defect, len(number_species), number_species, volume)
					for charge in extrinsic_defects_enthalpy_data[extrinsic_defect].keys():
						if synthesis_temperature is None:
							extrinsic_defect_charge_concentration = float(charge) * N * np.exp( -extrinsic_defects_enthalpy_data[extrinsic_defect][charge] / (k * temperature) )
						elif synthesis_temperature is not None:
							extrinsic_defect_charge_concentration = float(charge) * N * np.exp( -extrinsic_defects_enthalpy_data[extrinsic_defect][charge] / (k * synthesis_temperature) )
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
		fermi_energy_increment = ECBM
		while equilibrium_fermi_energy == 0.0:
			
			# Fermi energy 1
			fermi_energy1 = fermi_energy_increment
			
			# Hole concentration
			fE_holes1 = gE_ValenceBand * (1. - 1./( 1. + np.exp( (energies_ValenceBand - fermi_energy1) / (k * temperature) ) ) )
			hole_concentration1 = integrate.simps(fE_holes1, energies_ValenceBand)
			
			# Electron concentration
			fE_electrons1 = gE_ConductionBand / (1. + np.exp( (energies_ConductionBand - fermi_energy1) / (k * temperature) ) )
			electron_concentration1 = integrate.simps(fE_electrons1, energies_ConductionBand)
			
			# Defects concentration (mirroring Calculate_Defect_Carrier_Concentration() function)
			defects_carrier_concentration1 = 0
			intrinsic_defects_enthalpy_data = Calculate_IntrinsicDefectFormationEnthalpies(	defects_data, \
																							main_compound_total_energy, \
																							fermi_energy1, \
																							mu_elements )
			for intrinsic_defect in intrinsic_defects_enthalpy_data.keys():
				N = Find_SiteMultiplicity(intrinsic_defect, len(number_species), number_species, volume)
				for charge in intrinsic_defects_enthalpy_data[intrinsic_defect].keys():
					if synthesis_temperature is None:
						intrinsic_defect_charge_concentration = float(charge) * N * np.exp( -intrinsic_defects_enthalpy_data[intrinsic_defect][charge] / (k * temperature) )
					elif synthesis_temperature is not None:
						intrinsic_defect_charge_concentration = float(charge) * N * np.exp( -intrinsic_defects_enthalpy_data[intrinsic_defect][charge] / (k * synthesis_temperature) )
					defects_carrier_concentration1 += intrinsic_defect_charge_concentration
			if calculation_type == "total":
				if extrinsic_defect != "None":
					extrinsic_defects_enthalpy_data = Calculate_ExtrinsicDefectFormationEnthalpies(	defects_data, \
																									main_compound_total_energy, \
																									fermi_energy1, \
																									mu_elements, \
																									extrinsic_defect, \
																									extrinsic_defect_mu0, \
																									extrinsic_defect_deltamu )
					N = Find_SiteMultiplicity(extrinsic_defect, len(number_species), number_species, volume)
					for charge in extrinsic_defects_enthalpy_data[extrinsic_defect].keys():
						if synthesis_temperature is None:
							extrinsic_defect_charge_concentration = float(charge) * N * np.exp( -extrinsic_defects_enthalpy_data[extrinsic_defect][charge] / (k * temperature) )
						elif synthesis_temperature is not None:
							extrinsic_defect_charge_concentration = float(charge) * N * np.exp( -extrinsic_defects_enthalpy_data[extrinsic_defect][charge] / (k * synthesis_temperature) )
						defects_carrier_concentration1 += extrinsic_defect_charge_concentration
			
			# Total charge density at Fermi energy 1
			defect_charge_density1 = hole_concentration1 + electron_concentration1 + defects_carrier_concentration1
			
			# Increment Fermi energy (since we are searching the conduction band, inrease the energy to which to search)
			fermi_energy_increment += (ECBM-EVBM)/100.
			
			# Fermi energy 2 (incremented)
			fermi_energy2 = fermi_energy_increment
			
			# Hole concentration
			fE_holes2 = gE_ValenceBand * (1. - 1./( 1. + np.exp( (energies_ValenceBand - fermi_energy2) / (k * temperature) ) ) )
			hole_concentration2 = integrate.simps(fE_holes2, energies_ValenceBand)
			
			# Electron concentration
			fE_electrons2 = gE_ConductionBand / (1. + np.exp( (energies_ConductionBand - fermi_energy2) / (k * temperature) ) )
			electron_concentration2 = integrate.simps(fE_electrons2, energies_ConductionBand)
			
			# Defects concentration (mirroring Calculate_Defect_Carrier_Concentration() function)
			defects_carrier_concentration2 = 0
			intrinsic_defects_enthalpy_data = Calculate_IntrinsicDefectFormationEnthalpies(	defects_data, \
																							main_compound_total_energy, \
																							fermi_energy2, \
																							mu_elements )
			for intrinsic_defect in intrinsic_defects_enthalpy_data.keys():
				N = Find_SiteMultiplicity(intrinsic_defect, len(number_species), number_species, volume)
				for charge in intrinsic_defects_enthalpy_data[intrinsic_defect].keys():
					if synthesis_temperature is None:
						intrinsic_defect_charge_concentration = float(charge) * N * np.exp( -intrinsic_defects_enthalpy_data[intrinsic_defect][charge] / (k * temperature) )
					elif synthesis_temperature is not None:
						intrinsic_defect_charge_concentration = float(charge) * N * np.exp( -intrinsic_defects_enthalpy_data[intrinsic_defect][charge] / (k * synthesis_temperature) )
					defects_carrier_concentration2 += intrinsic_defect_charge_concentration
			if calculation_type == "total":
				if extrinsic_defect != "None":
					extrinsic_defects_enthalpy_data = Calculate_ExtrinsicDefectFormationEnthalpies(	defects_data, \
																									main_compound_total_energy, \
																									fermi_energy2, \
																									mu_elements, \
																									extrinsic_defect, \
																									extrinsic_defect_mu0, \
																									extrinsic_defect_deltamu )
					N = Find_SiteMultiplicity(extrinsic_defect, len(number_species), number_species, volume)
					for charge in extrinsic_defects_enthalpy_data[extrinsic_defect].keys():
						if synthesis_temperature is None:
							extrinsic_defect_charge_concentration = float(charge) * N * np.exp( -extrinsic_defects_enthalpy_data[extrinsic_defect][charge] / (k * temperature) )
						elif synthesis_temperature is not None:
							extrinsic_defect_charge_concentration = float(charge) * N * np.exp( -extrinsic_defects_enthalpy_data[extrinsic_defect][charge] / (k * synthesis_temperature) )
						defects_carrier_concentration2 += extrinsic_defect_charge_concentration
			
			# Total charge density at Fermi energy 2
			defect_charge_density2 = hole_concentration2 + electron_concentration2 + defects_carrier_concentration2
			
			# Check if Fermi energy admits equilibrium conditions
			if np.sign(defect_charge_density1) != np.sign(defect_charge_density2):
				equilibrium_fermi_energy = fermi_energy1
				equilibrium_hole_concentration = hole_concentration1
				equilibrium_electron_concentration = electron_concentration1
	
	return equilibrium_fermi_energy, equilibrium_hole_concentration, equilibrium_electron_concentration



