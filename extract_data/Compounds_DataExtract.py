
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'
__author__ = 'Michael_Jiaxing_Lidia_Benita_Elif'

##-------------------------------------------------------------------##
## This module sets up the initial information about atomic species, ## 
## main compound, competing compounds, formula units.                ##
##-------------------------------------------------------------------## 

import pandas as pd
import numpy as np
import sys
import os
from periodictable import elements


class Phase_Stability_Compounds_Information:
	
	def __init__(self):
		self.species = []			# Keep track of all unique elements involved in the system (e.g. for CuGaTe2, the elements would be Cu, Ga, and Te).
		self.compounds_info = {}	# Keep track of relevant information about all compounds (e.g. for CuGaTe2, compounds may include Cu, CuTe, Ga7Te10, ...)
		
		self.possible_elements = [str(el) for el in elements]	# To check whether element input by user is legitimate
		
		self.project_directory = os.getcwd()
	
	
	###-------------------------------------------------------------------------------###
	###-------------------------------------------------------------------------------###
	###----------- Initialize all elements and compounds in the material -------------###
	###-------------------------------------------------------------------------------###
	###-------------------------------------------------------------------------------###
	
	# Add elements included in the material.
	# The input can either be 1) a string of one element or 2) list of several elements.
	def Add_Species(self, species):
		if isinstance(species, str):					# Check if the input is a string.
			if species in self.possible_elements:		# This is a valid element.
				self.species.append(species)
			else:										# Whatever the user put in is jargon...
				print("WARNING: '"+species+"' is not recognized as a valid element.")
				pass
		elif isinstance(species, list):					# If the input is a list of species, then loop through the list.
			for specie in species:
				self.Add_Species(specie)
		else:											# If the input is neither a string nor list, then ignore it.
			print("WARNING: Species input is not valid. Change this to either a string or a list of species.")
			pass
	
	# Add compounds to the list of (preferably competing) compounds.
	# The input can either be 1) a string of one element or 2) list of several elements.
	# This simply initializes a data storage object (dictionary) for each compound.
	def Add_Compounds(self, new_compound):
		if isinstance(new_compound, str):				# Check if the input is a string.
			if new_compound in self.compounds_info.keys():
				print("WARNING: '"+new_compound+"' is already in the list of compounds! This compound will thus not be recorded.")
			elif new_compound not in os.listdir(self.project_directory):
				print("WARNING: '"+new_compound+"' DFT data cannot be found in this directory. Upload this data into this directory and try again!")
			else:
				self.compounds_info[new_compound] = {}
				number_species_tracker = []
				for specie in self.species:				# Find the number of species in the compound.
					if specie not in new_compound:
						self.compounds_info[new_compound][specie] = 0
						number_species_tracker.append(0)
					else:
						try:
							number_species_in_compound = int(new_compound.split(specie)[-1][:2])
						except:
							try:
								number_species_in_compound = int(new_compound.split(specie)[-1][:1])
							except:
								number_species_in_compound = 1
								pass
							pass
						self.compounds_info[new_compound][specie] = number_species_in_compound
						number_species_tracker.append(number_species_in_compound)
				if all(number_species == 0 for number_species in number_species_tracker):
					self.compounds_info.pop(new_compound, "None")
					print("WARNING: Since none of the species are present in '"+new_compound+"', this was deleted from the list of compounds.")
		elif isinstance(new_compound, list):			# Check if the input is a list.
			for compound in new_compound:
				self.Add_Compounds(compound)
		else:											# If the input is neither a string nor list, then ignore it.
			print("WARNING: Compounds input is not valid. Change this to either a string of a list of compounds.")
			pass
	
	# Add compounds to the list of (preferably competing) compounds from a list.
	# The input CSV file should have the following form:
	#				Element1	Element2	...
	#	Compound1	#			#			...
	#	Compound2	#			#			...
	#	   ...
	def Add_Compounds_From_File(self, compounds_filename):
		compounds_data = pd.read_csv(compounds_filename, header=0, index_col=0)
		header_information = list(compounds_data)										# Header will always be unique elements in the compound.
		for compound, compound_info in compounds_data.iterrows():						# Loop through compounds and their rows.
			if compound in self.compounds_info.keys():									# Check if the compound is already in the list.
				print("WARNING: '"+compound+"' is already in the list of compounds! This compound will thus not be recorded.")
			else:
				self.compounds_info[compound] = {}
				for header_specie, compound_data_info in zip(header_information, list(compound_info)):
					self.compounds_info[compound][header_specie] = compound_data_info	# Record number of each unique element in the compound in the list.
				for specie in self.species:
					if specie not in self.compounds_info[compound].keys():
						self.compounds_info[compound][specie] = 0						# If species is not in compound, record 0 in the list.
	
	
	###-------------------------------------------------------------------------------###
	###-------------------------------------------------------------------------------###
	###------------------------- Obtain DFT information ------------------------------###
	###-------------------------------------------------------------------------------###
	###-------------------------------------------------------------------------------###
	
	# Obtain the elements contained in the DFT POSCAR file of a specific compound.
	def Obtain_POSCAR_Formula(self, compound):
		try:
			open(compound+"/POSCAR", "r")				# Check if DFT data for specified compound exists in the directory.
		except:
			print("WARNING: '"+compound+"' POSCAR cannot be found in your directory.")
			return
		poscar_file = open(compound+"/POSCAR", "r")
		poscar_file_lines = poscar_file.readlines()
		poscar_elements = poscar_file_lines[5].split()	# Obtain the elements in the compound from the POSCAR.
		poscar_formula = poscar_file_lines[6].split()	# Obtain the number of each element in the compound from the POSCAR.
		return poscar_elements, poscar_formula
	
	# Update all POSCAR information in the list of compounds contained in this material object.
	def Update_POSCAR_Formulas(self):
		for compound in self.compounds_info.keys():										# Loop through each compound in the list.
			poscar_elements, poscar_formula = self.Obtain_POSCAR_Formula(compound)		# Invoke "Obtain_POSCAR_Formula" for each compound.
			for element, number_elements in zip(poscar_elements, poscar_formula):
				self.compounds_info[compound]["dft_"+element] = int(number_elements)	# Record under "dft_"(Element Name) how many of each element in involved in the DFT POSCAR.
			for specie in self.species:
				if specie not in poscar_elements:
					self.compounds_info[compound]["dft_"+specie] = 0					# If an element in the species list doesn't exist in the POSCAR, record 0 under "dft_"(Element Name)
	
	# Obtain the total energy of a compound from the finished DFT simulation OUTCAR.
	def Obtain_Total_Energy(self, compound):
		try:
			open(compound+"/OUTCAR", "r")							# Check if DFT data for the specified compound exists in the directory.
		except:
			print("WARNING: '"+compound+"' OUTCAR cannot be found in your directory.")
			return
		outcar_file = open(compound+"/OUTCAR", "r")
		for line in outcar_file.readlines():
			columns = line.split()
			number_columns = len(columns)
			if (number_columns > 4) and (columns[2] == "TOTEN"):	# Find the total energy of the compound.
				total_energy = float(columns[4])
		return total_energy
	
	# Update the total energy of each compound contained in this material.
	def Update_Total_Energies(self):
		for compound in self.compounds_info.keys():							# Loop through each compound in the list.
			total_energy = self.Obtain_Total_Energy(compound)				# Invoke "Obtain_Total_Energy" for each compound.
			self.compounds_info[compound]["total_energy"] = total_energy
	
	
	###-------------------------------------------------------------------------------###
	###-------------------------------------------------------------------------------###
	###-------------- Update compounds list using recorded quantities ----------------###
	###-------------------------------------------------------------------------------###
	###-------------------------------------------------------------------------------###
	
	# Update the number of formula units in the DFT calculation of each compound.
	def Update_Number_FormulaUnits(self):
		for compound in self.compounds_info.keys():
			for specie in self.species:
				if "dft_"+specie not in self.compounds_info[compound].keys():
					print("WARNING: DFT count of '"+specie+"' in '"+compound+"' cannot be found. Make sure you update the DFT data in the list of compounds!")
					print("Suggested command: Update_POSCAR_Formulas()")
				else:
					if self.compounds_info[compound][specie] != 0:
						self.compounds_info[compound]["formula_units"] = self.compounds_info[compound]["dft_"+specie]/self.compounds_info[compound][specie]
	
	# Update the total energy per formula unit of each compound
	def Update_EnergyPerFormulaUnit_Mu0(self):
		for compound in self.compounds_info.keys():
			if "total_energy" not in self.compounds_info[compound].keys():
				print("WARNING: Total energy of '"+compound+"' cannot be found. Make sure you update the total energy!")
				print("Suggested command: Update_Total_Energies()")
				continue
			elif "formula_units" not in self.compounds_info[compound].keys():
				print("WARNING: Number of formula units of '"+compound+"' cannot be found. Make sure you update the formula units!")
				print("Suggested command: Update_Number_FormulaUnits()")
				continue
			else:
				self.compounds_info[compound]["mu0"] = self.compounds_info[compound]["total_energy"] / self.compounds_info[compound]["formula_units"]
	
	# Update the enthalpies of each compound
	def Update_Enthalpies(self):
		for compound in self.compounds_info.keys():
			if "mu0" not in self.compounds_info[compound].keys():
				print("WARNING: Energy per formula unit (mu0) of "+compound+"cannot be found. Make sure you update mu0!")
				print("Suggested command: Update_EnergyPerFormulaUnit_Mu0()")
				continue
			enthalpy_tracker = self.compounds_info[compound]["mu0"]
			for specie in self.species:
				if self.compounds_info[compound][specie] == 0:
					continue
				enthalpy_tracker -= self.compounds_info[compound][specie]*self.compounds_info[specie]["mu0"]
			self.compounds_info[compound]["enthalpy"] = enthalpy_tracker
	
	# If necessary, clear any species that are not necessary
	def Clear_Unnecessary_Species(self):
		for specie in self.species:
			number_specie_tracker = []
			number_dft_specie_tracker = []
			for compound in self.compounds_info.keys():
				number_specie_tracker.append(self.compounds_info[compound][specie])
				number_dft_specie_tracker.append(self.compounds_info[compound]["dft_"+specie])
			if all(number_specie == 0.0 for number_specie in number_specie_tracker) and all(number_dft_specie == 0.0 for number_dft_specie in number_dft_specie_tracker):
				print("WARNING: "+specie+" is not present in any compound. This element will be wiped from this project.")
				for compound in self.compounds_info.keys():
					self.compounds_info[compound].pop(specie, "None")
					self.compounds_info[compound].pop("dft_"+specie, "None")
				self.species = self.species.remove(specie)
	
	def Update_All(self):
		self.Update_POSCAR_Formulas()
		self.Update_Total_Energies()
		self.Update_Number_FormulaUnits()
		self.Update_EnergyPerFormulaUnit_Mu0()
		self.Update_Enthalpies()
		self.Clear_Unnecessary_Species()
	
	# Update all information about the compounds in a CSV file.
	def Upload_Compounds_Info(self):
		compounds_dataframe = pd.DataFrame(self.compounds_info).transpose()
		compounds_dataframe.to_csv("Compounds_Tracker.csv")
	
	def Print_All_Info(self):
		print(pd.DataFrame(self.compounds_info).transpose())

