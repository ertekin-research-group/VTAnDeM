
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'
__author__ = 'Michael_Lidia_Jiaxing_Elif'

import numpy as np



def Calculate_IntrinsicDefectFormationEnthalpies(	defects_data, \
													main_compound_total_energy, \
													fermi_energy_array, \
													mu_elements	):
	
	# Initialize storage for all charges of all intrinsic defects
	intrinsic_defects_enthalpy_data = {}
	
	# Loop through defects in material
	for defect in defects_data.keys():
		
		# Check that the item is truly a defect
		if ("_" not in defect) and (defect.split("_")[-1] not in mu_elements.keys()):
			continue
		
		# Intrinsic defect formation_enthalpies
		if defects_data[defect]["Extrinsic"] == "No":
			intrinsic_defects_enthalpy_data[defect] = {}
			for charge in defects_data[defect]["charge"].keys():
				defect_formation_enthalpy = defects_data[defect]["charge"][charge]["Energy"] \
											- defects_data["supercellsize"] * main_compound_total_energy \
											+ float(charge) * fermi_energy_array \
											+ defects_data[defect]["charge"][charge]["ECorr"]
				for element in mu_elements.keys():
					defect_formation_enthalpy -= defects_data[defect]["n_"+element] * ( mu_elements[element]["mu0"] + mu_elements[element]["deltamu"] )
				intrinsic_defects_enthalpy_data[defect][charge] = defect_formation_enthalpy
	
	# Return dictionary of formation enthalpies of each defect
	return intrinsic_defects_enthalpy_data



def Calculate_ExtrinsicDefectFormationEnthalpies(	defects_data, \
													main_compound_total_energy, \
													fermi_energy_array, \
													mu_elements, \
													extrinsic_defects, \
													dopant_mu0, \
													dopant_deltamu	):
	
	# Check that the extrinsic defect name truly represents a defect
	for extrinsic_defect in extrinsic_defects:
		if ("_" not in extrinsic_defect) and (extrinsic_defect.split("_")[-1] not in mu_elements.keys()):
			return
	
	# Initialize storage for all charges of the extrinsic defect
	extrinsic_defects_enthalpy_data = {}
	
	for extrinsic_defect in extrinsic_defects:
		
		extrinsic_defects_enthalpy_data[extrinsic_defect] = {}
		
		# Loop through charge states of extrinsic defect
		for charge in defects_data[extrinsic_defect]["charge"].keys():
			defect_formation_enthalpy = defects_data[extrinsic_defect]["charge"][charge]["Energy"] \
										- defects_data["supercellsize"] * main_compound_total_energy \
										- (extrinsic_defect_mu0 + extrinsic_defect_deltamu) \
										+ float(charge) * fermi_energy_array \
										+ defects_data[extrinsic_defect]["charge"][charge]["ECorr"]
			for element in mu_elements.keys():
				# We subtract, since "defects_data[extrinsic_defect]["n_"+element]" is negative
				defect_formation_enthalpy -= defects_data[extrinsic_defect]["n_"+element] * ( mu_elements[element]["mu0"] + mu_elements[element]["deltamu"] )
			extrinsic_defects_enthalpy_data[extrinsic_defect][charge] = defect_formation_enthalpy
	
	return extrinsic_defects_enthalpy_data



def Find_MinimumDefectFormationEnthalpies(defect_formation_enthalpy_data):
	
	# Initialize storage for minimum formation enthalpies of each defect
	minimum_defect_formation_enthalpy_data = {}
	
	# Find minimum formation enthalpies
	for defect in defect_formation_enthalpy_data.keys():
		defect_formation_energy_minimum = np.fromiter(map(min, zip(*defect_formation_enthalpy_data[defect].values())), dtype=np.float)
		minimum_defect_formation_enthalpy_data[defect] = defect_formation_energy_minimum
	
	# Return dictionary of minimum formation enthalpies of each defect
	return minimum_defect_formation_enthalpy_data



def Find_SiteMultiplicity(defect_name, number_species, element_count, volume):
	
	# Check that the number of species matches the number of elements
	if number_species != len(element_count):
		return "Number of species does not match number of elemental species!"
	
	# Find matching site and corresponding stoichiometry of the site
	for species_index in range(number_species):
		if defect_name.split("_")[-1] == "i":	### CHANGE LATER FOR INTERSTITIAL ATOMS
			N = 0
		else:
			N = element_count[defect_name.split("_")[-1]] / volume
	
	return N










