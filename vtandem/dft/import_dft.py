
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'
__author__ = 'Michael_Jiaxing_Lidia_Benita_Elif'

########################################################################################################################
########################################################################################################################
#################################################### Introduction ######################################################
########################################################################################################################
########################################################################################################################
#
#
#
#



import os
import re
import json
import numpy as np
from shutil import copyfile
from periodictable import elements
from pymatgen.io.vasp.outputs import Vasprun


########################################################################################################################
########################################################################################################################
############################################# Phase Stability Information ##############################################
########################################################################################################################
########################################################################################################################

class Compounds_Import:
	
	def __init__(self):
		
		self.elements = [el.symbol for el in elements]
		
		# Recreate info in Compounds_Tracker.json as dictionary (if JSON exists), otherwise create new JSON
		if "Compounds_Tracker.json" in os.listdir(os.getcwd()):
			with open("Compounds_Tracker.json") as CompoundsTracker:
				self.compounds_info = json.load(CompoundsTracker)
		else:
			self.compounds_info = {"Compounds": {}, "Elements": {}}
	
	
	####################################################################################################################
	####################################################################################################################
	############################################# Add Element to Database ##############################################
	####################################################################################################################
	####################################################################################################################
	
	def Add_Element(self, element_name, directory_name):
		
		# Check if the directory name is legitimate
		if not os.path.isdir(directory_name):
			print("WARNING: Cannot find directory '"+directory_name+"'. Exiting...")
			return
		
		# Check to see that the element is legitimate
		if element_name not in self.elements:
			print(element_name+" is not recognized as a legitimate element. Exiting...")
			return ""
		
		# Check if the element already exists in the database
		if element_name in self.compounds_info["Elements"].keys():
			print(element_name+" is already in the database. The imported data will replace the old data.")
		
		# Initialize
		element_data = {}
		
		# Establish list of elements (should only be one item in list)
		element_data["elements_list"] = [element_name]
		
		# Find stoichiometry of element in POSCAR/CONTCAR
		if ("POSCAR" not in os.listdir(directory_name)) and ("CONTCAR" not in os.listdir(directory_name)):
			print("WARNING: Cannot find structure file (neither POSCAR nor CONTCAR) of "+element_name+". Exiting...")
			return ""
		try:
			structure_file = open(directory_name+"/POSCAR").readlines()
		except:
			structure_file = open(directory_name+"/CONTCAR").readlines()
			pass
		element_data[element_name] = 1.0
		element_data["dft_"+element_name] = float(structure_file[6].strip())
		element_data["formula_units"] = element_data["dft_"+element_name]
		element_data["number_species"] = 1
		
		# Find total energy of element in OUTCAR/OSZICAR
		if ("OUTCAR" not in os.listdir(directory_name)) and ("OSZICAR" not in os.listdir(directory_name)):
			print("WARNING: Cannot find total energy (neither OUTCAR nor OSZICAR) of "+element_name+". Exiting...")
			return ""
		try:
			oszicar_file = open(directory_name+"/OSZICAR").readlines()
			for line in oszicar_file:
				if line.split()[1] == "F=":
					total_energy = float(line.split()[2])
			element_data["total_energy"] = total_energy
		except:
			outcar_file = open(directory_name+"/OUTCAR").readlines()
			for line in outcar_file:
				columns = line.split()
				number_columns = len(columns)
				if (number_columns > 4) and (columns[2] == "TOTEN"):
					total_energy = float(columns[4])
			element_data["total_energy"] = total_energy
			pass
		
		# Calculate total energy per formula unit of compound
		element_data["mu0"] = element_data["total_energy"] / element_data["formula_units"]
		
		# Set enthalpy to zero
		element_data["enthalpy"] = 0.0
		
		# Store data
		self.compounds_info["Elements"][element_name] = element_data
	
	
	####################################################################################################################
	####################################################################################################################
	############################################# Add Compound to Database #############################################
	####################################################################################################################
	####################################################################################################################
	
	def Add_Compound(self, compound_name, directory_name):
		
		# Check if the directory name is legitimate
		if not os.path.isdir(directory_name):
			print("WARNING: Cannot find directory '"+directory_name+"'. Exiting...")
			return
		
		# Check if the compound is actually an element
		if compound_name in self.elements:
			print("WARNING: '"+compound_name+"' is an element, not a compound. Exiting...")
			return
			
		
		# Check if the compound already exists in the database
		if compound_name in self.compounds_info["Compounds"].keys():
			print(compound_name+" is already in the phase stability database (Compounds_Tracker.json). The imported data will replace the old data.")
		
		# Initialize
		compound_data = {}
		
		# Determine elements in compound from compound_name
		elements_list = [ ''.join( [letter for letter in element_segment if not letter.isdigit()] ) for element_segment in re.findall("[A-Z][^A-Z]*", compound_name) ]
		for element in elements_list:
			if element not in self.elements:
				print("'"+element+"' is not a valid element. Exiting...")
				return ""
			else:
				try:
					number_species_in_compound = float(compound_name.split(element)[-1][:2])
				except:
					try:
						number_species_in_compound = float(compound_name.split(element)[-1][:1])
					except:
						number_species_in_compound = 1.0
						pass
					pass
				compound_data[element] = number_species_in_compound
		
		# List of elements in compound
		compound_data["elements_list"] = elements_list
		
		# Find stoichiometry of compound listed in POSCAR/CONTCAR
		if ("POSCAR" not in os.listdir(directory_name)) and ("CONTCAR" not in os.listdir(directory_name)):
			print("WARNING: Cannot find structure file (neither POSCAR nor CONTCAR). Exiting...")
			return ""
		try:
			structure_file = open(directory_name+"/POSCAR").readlines()
		except:
			structure_file = open(directory_name+"/CONTCAR").readlines()
			pass
		species = structure_file[5].split()
		number_species = structure_file[6].split()
		for element, number in zip(species, number_species):
			compound_data["dft_"+element] = float(number)
			compound_data["formula_units"] = compound_data["dft_"+element]/compound_data[element]
		compound_data["number_species"] = len(number_species)
		
		# Find total energy of compound in OUTCAR/OSZICAR
		if ("OUTCAR" not in os.listdir(directory_name)) and ("OSZICAR" not in os.listdir(directory_name)):
			print("WARNING: Cannot find total energy (neither OUTCAR nor OSZICAR). Exiting...")
			return ""
		try:
			oszicar_file = open(directory_name+"/OSZICAR").readlines()
			for line in oszicar_file:
				if line.split()[1] == "F=":
					total_energy = float(line.split()[2])
			compound_data["total_energy"] = total_energy
		except:
			outcar_file = open(directory_name+"/OUTCAR").readlines()
			for line in outcar_file:
				columns = line.split()
				number_columns = len(columns)
				if (number_columns > 4) and (columns[2] == "TOTEN"):
					total_energy = float(columns[4])
			compound_data["total_energy"] = total_energy
			pass
		
		# Calculate total energy per formula unit of compound
		compound_data["mu0"] = compound_data["total_energy"] / compound_data["formula_units"]
		
		# Calculate enthalpy of formation of compound
		enthalpy_tracker = compound_data["mu0"]
		for element in elements_list:
			enthalpy_tracker -= compound_data[element]*self.compounds_info["Elements"][element]["mu0"]
		compound_data["enthalpy"] = enthalpy_tracker
		
		# Find band gap and valence band maximum of compound in vasprun.xml
		if "vasprun.xml" not in os.listdir(directory_name):
			print("WARNING: Cannot find band gap nor valence band maximum (vasprun.xml file). Exiting...")
			return ""
		vasprun = Vasprun(directory_name+"/vasprun.xml")
		(bandgap, cbm, vbm, is_direct) = vasprun.eigenvalue_band_properties
		compound_data["band_gap"] = bandgap
		compound_data["vbm"] = vbm
		
		# Store data
		self.compounds_info["Compounds"][compound_name] = compound_data
	
	
	####################################################################################################################
	####################################################################################################################
	################################################# Update Database ##################################################
	####################################################################################################################
	####################################################################################################################
	
	def Update_Compounds_Database(self):
		
		try:
			copyfile("Compounds_Tracker.json", ".vtandem/Compounds_Tracker_Backup.json")
		except:
			pass
		with open("Compounds_Tracker.json", "w") as jsonfile:
			json.dump(self.compounds_info, jsonfile, indent=4, sort_keys=True)







########################################################################################################################
########################################################################################################################
################################################# Defects Information ##################################################
########################################################################################################################
########################################################################################################################

class Defects_Import:
	
	def __init__(self):
		
		self.elements = [el.symbol for el in elements]
		
		# Recreate defect information as dictionary (if JSON exists), otherwise create new datasheet
		if "Defects_Tracker.json" in os.listdir(os.getcwd()):
			with open("Defects_Tracker.json") as DefectsTracker:
				self.defects_data = json.load(DefectsTracker)
		else:
			self.defects_data = {}
	
	
	####################################################################################################################
	####################################################################################################################
	####################################### Add Defects of Compound to Database ########################################
	####################################################################################################################
	####################################################################################################################
	
	def Add_Defects(self, compound_name, directory_name, supercell_size):
		
		# If compound does not exists in Compounds_Tracker.json, then don't import
		if compound_name not in json.load(open("Compounds_Tracker.json"))["Compounds"].keys():
			print("The compound '"+compound_name+"' does not exist in Compounds_Tracker.json. Skipping...")
			return
		
		# Create list of possible defects in compound
		possible_defect_site_list = [ ''.join( [letter for letter in element_segment if not letter.isdigit()] ) for element_segment in re.findall("[A-Z][^A-Z]*", compound_name) ]
		possible_defect_site_list.append("V")
		
		# Initialize
		defects_data = {}
		
		# Loop through each item in the given directory
		for directory in os.listdir(directory_name):
			
			# Check if the item is indeed a directory
			if not os.path.isdir(directory_name+"/"+directory):
				continue
			
			# If item is not a legitimate defect, skip it
			try:
				directory.split("_")
			except:
				print("The directory name '"+directory+"' is not in the correct format for the defect name. Skipping...")
				continue
			if (directory.split("_")[-1] not in possible_defect_site_list) or (directory.split("_")[0] not in self.elements):
				print("'"+directory+"' is not a legitimate defect name for "+compound_name+". Skipping...")
				continue
			
			defect_name = directory
			defects_data[defect_name] = {}
			
			# Extrinsic dopant?
			if (directory.split("_")[-1]  in possible_defect_site_list) and (directory.split("_")[0] not in possible_defect_site_list):
				defects_data[defect_name]["Extrinsic"] = "Yes"
			else:
				defects_data[defect_name]["Extrinsic"] = "No"
			
			# Obtain "stoichiometry" of defect
			for specie in possible_defect_site_list:
				if specie == "V":
					if defect_name.split("_")[0] == specie:
						defects_data[defect_name]["n_"+defect_name.split("_")[-1]] = -1
				elif specie not in defect_name.split("_"):
					defects_data[defect_name]["n_"+specie] = 0
				else:
					if defect_name.split("_")[0] == specie:
						defects_data[defect_name]["n_"+specie] = 1
					elif defect_name.split("_")[-1] == specie:
						defects_data[defect_name]["n_"+specie] = -1
			
			
			# Get charges and total energy of defect
			defects_data[defect_name]["charge"] = {}
			for subdirectory in os.listdir(directory_name+"/"+defect_name):
				
				# Check to see that the folder is a legitimate charge
				try:
					float(subdirectory.split("q")[-1])
				except:
					print("The directory name '"+subdirectory+"' is not in the correct format for charge. Skipping...")
					continue
				
				charge_state = subdirectory
				
				# Find total energy of defect in OUTCAR/OSZICAR
				if ("OUTCAR" not in os.listdir(directory_name+"/"+defect_name+"/"+charge_state)) and ("OSZICAR" not in os.listdir(directory_name+"/"+defect_name+"/"+charge_state)):
					print("WARNING: Cannot find total energy (neither OUTCAR nor OSZICAR) for '"+defect_name+"' with charge state '"+charge_state+"'. Skipping...")
					continue
				try:
					outcar_file = open(directory_name+"/"+defect_name+"/"+charge_state+"/OUTCAR").readlines()
					for line in outcar_file:
						columns = line.split()
						number_columns = len(columns)
						if (number_columns > 4) and (columns[2] == "TOTEN"):
							total_energy = float(columns[4])
				except:
					oszicar_file = open(directory_name+"/"+defect_name+"/"+charge_state+"/OSZICAR").readlines()
					for line in oszicar_file:
						if "F=" in line:
							total_energy = float(line.split()[2])
					pass
				
				# Store total energy of defect in charge state q
				defects_data[defect_name]["charge"][charge_state.split("q")[-1]] = {}
				defects_data[defect_name]["charge"][charge_state.split("q")[-1]]["Energy"] = total_energy
				
				# UPDATE THIS SECTION LATER
				defects_data[defect_name]["charge"][charge_state.split("q")[-1]]["ECorr"] = 0.0
			
			
			# If no charges were recorded, then delete the defect from the database altogether
			if defects_data[defect_name]["charge"] == {}:
				print("No charges for '"+defect_name+"' were recorded. Erasing data for '"+defect_name+"'.")
				del defects_data[defect_name]
		
		# Store data
		self.defects_data[compound_name] = defects_data
		self.defects_data[compound_name]["supercellsize"] = supercell_size
	
	
	####################################################################################################################
	####################################################################################################################
	############################################# Update Defects Database ##############################################
	####################################################################################################################
	####################################################################################################################
	
	def Update_Defects_Database(self):
		try:
			copyfile("Defects_Tracker.json", ".vtandem/Defects_Tracker_Backup.json")
		except:
			pass
		with open("Defects_Tracker.json", "w") as jsonfile:
			json.dump(self.defects_data, jsonfile, indent=4, sort_keys=True)







########################################################################################################################
########################################################################################################################
################################################### DOS Information ####################################################
########################################################################################################################
########################################################################################################################

class DOS_Import:
	
	def __init__(self):
		
		if "DOS_Tracker.json" in os.listdir(os.getcwd()):
			with open("DOS_Tracker.json") as DOSTracker:
				self.dos_data = json.load(DOSTracker)
		else:
			self.dos_data = {}
	
	
	####################################################################################################################
	####################################################################################################################
	############################################### Add DOS to Database ################################################
	####################################################################################################################
	####################################################################################################################
	
	def Add_DOS(self, compound_name, filename):
		
		# Check if the directory name is legitimate
		if not os.path.isfile(filename):
			print("WARNING: Cannot find file: '"+directory_name+"'. Exiting...")
			return ""
		
		# Check if the compound already exists in the database
		if compound_name in self.dos_data.keys():
			print("The compound '"+compound_name+"' is already in the DOS database (DOS_Tracker.json). The imported data will replace the old data.")
		
		# If compound does not exists in Compounds_Tracker.json, then don't import
		if compound_name not in json.load(open("Compounds_Tracker.json"))["Compounds"].keys():
			print("The compound '"+compound_name+"' does not exist in Compounds_Tracker.json. Skipping...")
			return
		
		# Initialize
		dos_info = {}
		
		# Extract data
		for line in open(filename).readlines():
			
			#Skip lines that describe the partial DOS with respect to orbital projections
			if len(line.split()) != 3:
				continue
			
			try:
				dos_info[float(line.split()[0])] = float(line.split()[1])
			except:
				dos_info[float(line.split()[0])] = 0.0	# Sometimes the DOS can be an extremely small number that VASP outputs e.g. "0.5E-111" as "0.5-111", which is not readable by python
			
		# Store data
		self.dos_data[compound_name] = dos_info
	
	
	####################################################################################################################
	####################################################################################################################
	############################################### Update DOS Database ################################################
	####################################################################################################################
	####################################################################################################################
	
	def Update_DOS_Database(self):
		
		try:
			copyfile("DOS_Tracker.json", ".vtandem/DOS_Tracker_Backup.json")
		except:
			pass
		with open("DOS_Tracker.json", "w") as jsonfile:
			json.dump(self.dos_data, jsonfile, indent=4, sort_keys=True)





