
import numpy as np
import periodictable
import copy
import itertools


all_elements = []
for element in periodictable.elements:
	all_elements.append(str(element))




def Calculate_PhaseDiagram_Projected2D(main_compound, elements_dict: dict, compounds_info, deltamu: dict):
	
	# Plotting the phase diagram from DFT data is one of the main features of this app. This function
	#	single-handedly plots the phase diagram, so it's arguably one of the most important block of
	#	code in this script.
	# This function will be used in two forms:
	#	1) when the user generates the phase diagram by clicking the "Generate Plot" button
	#	2) when the user changes the mu4 value
	# Because there are different forms in which this function will be used, the programmer should take
	#	care when editing this block of code.
	
	
	
	# Number of elements in main_compound
	main_compound_elements_count = {}
	for element_index in sorted(elements_dict.keys()):
		main_compound_elements_count[element_index] = compounds_info[main_compound][elements_dict[element_index]]
	
	
	# Enthalpy of compound
	main_compound_enthalpy = compounds_info[main_compound]["enthalpy"]
	main_compound_enthalpy_adjusted = main_compound_enthalpy	# Enthalpy adjusted for mu4 value (main compound)
	for element_index in sorted(elements_dict.keys()):
		if element_index in [1, 2, 3]:
			continue
		main_compound_enthalpy_adjusted -= main_compound_elements_count[element_index]*deltamu[element_index]
	#main_compound_enthalpy_adjusted = self.main_compound_enthalpy - self.main_compound_number_fourth_specie*self.mu4	# Enthalpy adjusted for mu4 value (main compound)
	
	main_compound_deltamu_first_element = np.linspace(main_compound_enthalpy/main_compound_elements_count[1], 0, 1000)	# Array of possible mu1 values in the plot
	main_compound_deltamu_second_element = (main_compound_enthalpy_adjusted - main_compound_elements_count[1]*main_compound_deltamu_first_element) / main_compound_elements_count[2]	# Array that represents the (diagonal) stability limit of the main compound
	
	
	# Find the bounds of the quaternary phase stability region (for shading the stability region of the quaternary compound)
	stability_minimum_bound = []
	stability_maximum_bound = []
	
	vertical_left_values = []
	vertical_right_values = []
	
	competing_compounds_deltamu_first_element_limit = {}
	competing_compounds_deltamu_second_element_limit = {}
	
	# Loop through all compounds in the database
	for competing_compound in compounds_info.keys():
		
		# Skip if compound is either the main compound or one of the elements
		if (competing_compound in all_elements) or (competing_compound == main_compound):
			continue
		
		competing_compound_elements_count = {}
		for element_index in sorted(elements_dict.keys()):
			try:
				competing_compound_elements_count[element_index] = compounds_info[competing_compound][elements_dict[element_index]]
			except:
				competing_compound_elements_count[element_index] = 0.0
				pass
		
		competing_compound_enthalpy = compounds_info[competing_compound]["enthalpy"]
		#competing_compound_enthalpy_adjusted = competing_compound_enthalpy - competing_compound_number_fourth_specie*self.mu4		# Enthalpy adjusted for mu4 value (competing compound)
		competing_compound_enthalpy_adjusted = competing_compound_enthalpy		# Enthalpy adjusted for mu4 value (competing compound)
		for element_index in sorted(elements_dict.keys()):
			if element_index in [1, 2, 3]:
				continue
			competing_compound_enthalpy_adjusted -= competing_compound_elements_count[element_index]*deltamu[element_index]
		
		
		difference_enthalpy_adjusted = competing_compound_enthalpy_adjusted - (competing_compound_elements_count[3]/main_compound_elements_count[3])*main_compound_enthalpy_adjusted
		
		coefficient_first_specie = competing_compound_elements_count[1] - (main_compound_elements_count[1]*competing_compound_elements_count[3]) / main_compound_elements_count[3]
		coefficient_second_specie = competing_compound_elements_count[2] - (main_compound_elements_count[2]*competing_compound_elements_count[3]) / main_compound_elements_count[3]
		
		if (coefficient_first_specie == 0.0) and (coefficient_second_specie == 0.0):
			continue
			print("OMG SOMETHING'S WRONG!!!")
			print("Compound may not be stoichiometrically balanced???")
			
		elif (coefficient_first_specie != 0.0) and (coefficient_second_specie == 0.0):		# Vertical line
			constant_deltamu_first_element = difference_enthalpy_adjusted / coefficient_first_specie
			if coefficient_first_specie > 0.0:
				vertical_left_values.append(constant_deltamu_first_element)
			elif coefficient_first_specie < 0.0:
				vertical_right_values.append(constant_deltamu_first_element)
			
			competing_compound_deltamu_first_element = np.ones(len(main_compound_deltamu_first_element)) * constant_deltamu_first_element
			competing_compound_deltamu_second_element = np.linspace(main_compound_enthalpy/main_compound_elements_count[2], 0, len(competing_compound_deltamu_first_element))
			
			main_compound_stability_limit_vertical = ( main_compound_enthalpy_adjusted - main_compound_elements_count[1]*constant_deltamu_first_element ) / main_compound_elements_count[2]
			competing_compound_deltamu_first_element_limit = [competing_compound_deltamu_first_element[i] for i in range(len(competing_compound_deltamu_first_element)) if (main_compound_stability_limit_vertical < competing_compound_deltamu_second_element[i])]
			competing_compound_deltamu_second_element_limit = [competing_compound_deltamu_second_element[i] for i in range(len(competing_compound_deltamu_first_element)) if (main_compound_stability_limit_vertical < competing_compound_deltamu_second_element[i])]
			
		elif (coefficient_first_specie == 0.0) and (coefficient_second_specie != 0.0):		# Horizontal line
			competing_compound_deltamu_first_element = copy.deepcopy(main_compound_deltamu_first_element)
			constant_deltamu_second_element = difference_enthalpy_adjusted / coefficient_second_specie
			competing_compound_deltamu_second_element = np.ones(len(competing_compound_deltamu_first_element)) * constant_deltamu_second_element
			
			if coefficient_second_specie > 0.0:
				stability_maximum_bound.append(competing_compound_deltamu_second_element)
			elif coefficient_second_specie < 0.0:
				stability_minimum_bound.append(competing_compound_deltamu_second_element)
			
			competing_compound_deltamu_first_element_limit = [competing_compound_deltamu_first_element[i] for i in range(len(competing_compound_deltamu_first_element)) if (main_compound_deltamu_second_element[i] < competing_compound_deltamu_second_element[i])]
			competing_compound_deltamu_second_element_limit = [competing_compound_deltamu_second_element[i] for i in range(len(competing_compound_deltamu_first_element)) if (main_compound_deltamu_second_element[i] < competing_compound_deltamu_second_element[i])]
		
		elif (coefficient_first_specie != 0.0) and (coefficient_second_specie != 0.0):
			competing_compound_deltamu_first_element = copy.deepcopy(main_compound_deltamu_first_element)
			competing_compound_deltamu_second_element = ( difference_enthalpy_adjusted - coefficient_first_specie*competing_compound_deltamu_first_element ) / coefficient_second_specie
			
			if coefficient_second_specie > 0.0:
				stability_maximum_bound.append(competing_compound_deltamu_second_element)
			elif coefficient_second_specie < 0.0:
				stability_minimum_bound.append(competing_compound_deltamu_second_element)
			
			competing_compound_deltamu_first_element_limit = [competing_compound_deltamu_first_element[i] for i in range(len(competing_compound_deltamu_first_element)) if (main_compound_deltamu_second_element[i] < competing_compound_deltamu_second_element[i])]
			competing_compound_deltamu_second_element_limit = [competing_compound_deltamu_second_element[i] for i in range(len(competing_compound_deltamu_first_element)) if (main_compound_deltamu_second_element[i] < competing_compound_deltamu_second_element[i])]
		
		competing_compounds_deltamu_first_element_limit[competing_compound] = competing_compound_deltamu_first_element_limit
		competing_compounds_deltamu_second_element_limit[competing_compound] = competing_compound_deltamu_second_element_limit
	
	
	stability_minimum_bound.append(main_compound_deltamu_second_element)
	stability_maximum_bound.append(np.zeros(len(main_compound_deltamu_first_element)))
	stability_absolute_minimum = np.fromiter(map(max, zip(*itertools.chain(stability_minimum_bound))), dtype=np.float)
	stability_absolute_maximum = np.fromiter(map(min, zip(*itertools.chain(stability_maximum_bound))), dtype=np.float)
	
	main_compound_deltamu_first_element_cutoff = []
	stability_minimum_cutoff = []
	stability_maximum_cutoff = []
	
	for i in range(len(main_compound_deltamu_first_element)):
		if (vertical_left_values != []) and (vertical_right_values != []):
			if (stability_absolute_minimum[i] < stability_absolute_maximum[i]) and (main_compound_deltamu_first_element[i] < min(vertical_left_values)) and (main_compound_deltamu_first_element[i] > max(vertical_right_values)):
				main_compound_deltamu_first_element_cutoff.append(main_compound_deltamu_first_element[i])
				stability_minimum_cutoff.append(stability_absolute_minimum[i])
				stability_maximum_cutoff.append(stability_absolute_maximum[i])
		if (vertical_left_values != []) and (vertical_right_values == []):
			if (stability_absolute_minimum[i] < stability_absolute_maximum[i]) and (main_compound_deltamu_first_element[i] < min(vertical_left_values)):
				main_compound_deltamu_first_element_cutoff.append(main_compound_deltamu_first_element[i])
				stability_minimum_cutoff.append(stability_absolute_minimum[i])
				stability_maximum_cutoff.append(stability_absolute_maximum[i])
		if (vertical_left_values == []) and (vertical_right_values != []):
			if (stability_absolute_minimum[i] < stability_absolute_maximum[i]) and (main_compound_deltamu_first_element[i] > max(vertical_right_values)):
				main_compound_deltamu_first_element_cutoff.append(main_compound_deltamu_first_element[i])
				stability_minimum_cutoff.append(stability_absolute_minimum[i])
				stability_maximum_cutoff.append(stability_absolute_maximum[i])
		if (vertical_left_values == []) and (vertical_right_values == []):
			if (stability_absolute_minimum[i] < stability_absolute_maximum[i]):
				main_compound_deltamu_first_element_cutoff.append(main_compound_deltamu_first_element[i])
				stability_minimum_cutoff.append(stability_absolute_minimum[i])
				stability_maximum_cutoff.append(stability_absolute_maximum[i])
	
	
	
	return 	main_compound_deltamu_first_element, \
			main_compound_deltamu_second_element, \
			competing_compounds_deltamu_first_element_limit, \
			competing_compounds_deltamu_second_element_limit, \
			main_compound_deltamu_first_element_cutoff, \
			stability_minimum_cutoff, \
			stability_maximum_cutoff













