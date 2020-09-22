
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'
__author__ = 'Michael_Jiaxing_Lidia_Elif'

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
	############################################# Add Element to Database ##############################################
	####################################################################################################################
	
	def Add_Element(self, element_name, directory_name):
		
		# Check if the directory name is legitimate
		if not os.path.isdir(directory_name):
			print("WARNING: Cannot find directory '"+directory_name+"'. Exiting...")
			return False
		
		# Check to see that the element is legitimate
		if element_name not in self.elements:
			print(element_name+" is not recognized as a legitimate element. Exiting...")
			return False
		
		# Check if the element already exists in the database
		if element_name in self.compounds_info["Elements"].keys():
			print(element_name+" is already in the database. The imported data will replace the old data.")
			pass
		
		# Initialize
		element_data = {}
		
		# Establish list of elements (should only be one item in list)
		element_data["elements_list"] = [element_name]
		
		# Find stoichiometry of element in POSCAR/CONTCAR
		if ("POSCAR" not in os.listdir(directory_name)) and ("CONTCAR" not in os.listdir(directory_name)):
			print("WARNING: Cannot find structure file (neither POSCAR nor CONTCAR) of "+element_name+". Exiting...")
			return False
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
			return False
		try:
			outcar_file = open(directory_name+"/OUTCAR").readlines()
			for line in outcar_file:
				columns = line.split()
				number_columns = len(columns)
				if (number_columns > 4) and (columns[2] == "TOTEN"):
					total_energy = float(columns[4])
			element_data["total_energy"] = total_energy
		except:
			oszicar_file = open(directory_name+"/OSZICAR").readlines()
			for line in oszicar_file:
				if line.split()[1] == "F=":
					total_energy = float(line.split()[2])
			element_data["total_energy"] = total_energy
			pass
		
		# Calculate total energy per formula unit of compound
		element_data["mu0"] = element_data["total_energy"] / element_data["formula_units"]
		
		# Set enthalpy to zero
		element_data["enthalpy"] = 0.0
		
		# Store data
		self.compounds_info["Elements"][element_name] = element_data
		
		return True
	
	
	####################################################################################################################
	############################################# Add Compound to Database #############################################
	####################################################################################################################
	
	def Add_Compound(self, compound_name, directory_name):
		
		# Check if the directory name is legitimate
		if not os.path.isdir(directory_name):
			print("WARNING: Cannot find directory '"+directory_name+"'. Exiting...")
			return False
		
		# Check if the compound is actually an element
		if compound_name in self.elements:
			print("WARNING: '"+compound_name+"' is an element, not a compound. Exiting...")
			return False
		
		# Check if the compound already exists in the database
		if compound_name in self.compounds_info["Compounds"].keys():
			print(compound_name+" is already in the database. The imported data will replace the old data.")
			pass
		
		# Initialize
		compound_data = {}
		
		# Determine elements in compound from compound_name
		elements_list = [ ''.join( [letter for letter in element_segment if not letter.isdigit()] ) for element_segment in re.findall("[A-Z][^A-Z]*", compound_name) ]
		for element in elements_list:
			if element not in self.elements:
				print("'"+element+"' is not a valid element. Exiting...")
				return False
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
			return False
		try:
			structure_file = open(directory_name+"/CONTCAR").readlines()
		except:
			structure_file = open(directory_name+"/POSCAR").readlines()
			pass
		species = structure_file[5].split()
		number_species = structure_file[6].split()
		for element, number in zip(species, number_species):
			compound_data["dft_"+element] = float(number)
			compound_data["formula_units"] = compound_data["dft_"+element]/compound_data[element]
		compound_data["number_species"] = len(number_species)
		
		# Find volume of compound from lattice vectors in POSCAR/CONTCAR
		lattice_vector_x = np.asarray([float(i) for i in structure_file[2].split()])
		lattice_vector_y = np.asarray([float(i) for i in structure_file[3].split()])
		lattice_vector_z = np.asarray([float(i) for i in structure_file[4].split()])
		volume = np.dot(lattice_vector_x, np.cross(lattice_vector_y, lattice_vector_z))
		compound_data["volume"] = volume*1E-24
		
		# Find total energy of compound in OUTCAR/OSZICAR
		if ("OUTCAR" not in os.listdir(directory_name)) and ("OSZICAR" not in os.listdir(directory_name)):
			print("WARNING: Cannot find total energy (neither OUTCAR nor OSZICAR). Exiting...")
			return False
		try:
			outcar_file = open(directory_name+"/OUTCAR").readlines()
			for line in outcar_file:
				columns = line.split()
				number_columns = len(columns)
				if (number_columns > 4) and (columns[2] == "TOTEN"):
					total_energy = float(columns[4])
			compound_data["total_energy"] = total_energy
		except:
			oszicar_file = open(directory_name+"/OSZICAR").readlines()
			for line in oszicar_file:
				if line.split()[1] == "F=":
					total_energy = float(line.split()[2])
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
			return False
		vasprun = Vasprun(directory_name+"/vasprun.xml")
		(bandgap, cbm, vbm, is_direct) = vasprun.eigenvalue_band_properties
		compound_data["band_gap"] = bandgap
		compound_data["vbm"] = vbm
		
		# Store data
		self.compounds_info["Compounds"][compound_name] = compound_data
		
		return True
	
	
	####################################################################################################################
	################################################# Update Database ##################################################
	####################################################################################################################
	
	def Update_Compounds_Database(self):
		
		try:
			copyfile("Compounds_Tracker.json", ".vtandem/Compounds_Tracker_Backup.json")
		except:
			pass
		
		"""
		with open("Compounds_Tracker.json", "w") as jsonfile:
			json.dump(self.compounds_info, jsonfile, indent=4, sort_keys=True)
		"""
		jsonfile = open("Compounds_Tracker.json", "w")
		json.dump(self.compounds_info, jsonfile, indent=4, sort_keys=True)
		jsonfile.close()







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
	####################################### Add Defects of Compound to Database ########################################
	####################################################################################################################
	
	def Add_Defects(self, compound_name, directory_name, supercell_size):
		
		if not self.Run_Checks(compound_name, directory_name):
			return
		
		# Create list of possible defects in compound
		possible_defect_site_list = [ ''.join( [letter for letter in element_segment if not letter.isdigit()] ) for element_segment in re.findall("[A-Z][^A-Z]*", compound_name) ]
		possible_defect_site_list.append("V")
		possible_defect_site_list.append("i")
		
		# Loop through each item in the given directory
		for directory in os.listdir(directory_name):
			
			# Check that the directory is a legitimate defect name
			if (directory.split("_")[-1] not in possible_defect_site_list) or (directory.split("_")[0] not in self.elements):
				print("'"+directory+"' in '"+directory_name+"' is not a legitimate defect name for '"+compound_name+"'. Skipping...")
				continue
			
			self.Add_Single_Defect(compound_name=compound_name, defect_name=directory, directory_name=directory_name+"/"+directory, supercell_size=supercell_size)
	
	
	####################################################################################################################
	############################################# Add an Individual Defect #############################################
	####################################################################################################################
	
	def Add_Single_Defect(self, compound_name, defect_name, directory_name, supercell_size):
		
		if not self.Run_Checks(compound_name, directory_name):
			return
		
		# Create list of possible defects in compound
		possible_defect_site_list = [ ''.join( [letter for letter in element_segment if not letter.isdigit()] ) for element_segment in re.findall("[A-Z][^A-Z]*", compound_name) ]
		possible_defect_site_list.append("V")
		possible_defect_site_list.append("i")
		
		# Check that the directory is a legitimate defect name
		if (defect_name.split("_")[-1] not in possible_defect_site_list) or (defect_name.split("_")[0] not in self.elements):
			print("'"+defect_name+"' is not a legitimate defect name for '"+compound_name+"'. Skipping...")
			return
		
		for directory in os.listdir(directory_name):
			try:
				float(directory.split("q")[-1])
			except:
				print("The name '"+directory+"' in '"+directory_name+"' is not in the correct format for charge. The correct format for the charge state of a defect is 'q#', where '#' is an integer representing the charge state. Skipping...")
				continue
			self.Add_Defect_Charge(compound_name=compound_name, defect_name=defect_name, charge_state=directory.split("q")[-1], directory_name=directory_name+"/"+directory, supercell_size=supercell_size)
	
	
	
	####################################################################################################################
	###################################### Add Charge State for Individual Defect ######################################
	####################################################################################################################
	
	def Add_Defect_Charge(self, compound_name, defect_name, charge_state, directory_name, supercell_size):
		
		if not self.Run_Checks(compound_name, directory_name):
			return
		
		# Create list of possible defects in compound
		possible_defect_site_list = [ ''.join( [letter for letter in element_segment if not letter.isdigit()] ) for element_segment in re.findall("[A-Z][^A-Z]*", compound_name) ]
		possible_defect_site_list.append("V")
		possible_defect_site_list.append("i")
		
		# Check that the given defect is a legitimate defect name
		if (defect_name.split("_")[-1] not in possible_defect_site_list) or (defect_name.split("_")[0] not in self.elements):
			print("'"+defect_name+"' is not a legitimate defect name for "+compound_name+". Skipping...")
			return
		
		# Check that the data exists
		if ("OUTCAR" not in os.listdir(directory_name)) and ("OSZICAR" not in os.listdir(directory_name)):
			print("WARNING: Cannot find total energy (neither OUTCAR nor OSZICAR files) for defect '"+defect_name+"' with charge state '"+charge_state+"' in '"+directory_name+"'. Skipping")
			return
		
		# Set up data for compound/defects in Defects_Tracker.json
		if compound_name not in self.defects_data.keys():
			self.defects_data[compound_name] = {}
			self.defects_data[compound_name]["supercellsize"] = supercell_size
		if defect_name not in self.defects_data[compound_name].keys():
			self.defects_data[compound_name][defect_name] = {}
			self.Get_Basic_Defect_Info(compound_name, defect_name, possible_defect_site_list)
		if "charge" not in self.defects_data[compound_name][defect_name].keys():
			self.defects_data[compound_name][defect_name]["charge"] = {}
		if charge_state not in self.defects_data[compound_name][defect_name]["charge"].keys():
			self.defects_data[compound_name][defect_name]["charge"][charge_state] = {}
		
		# Find total energy of defect with specified charge state in OUTCAR or OSZICAR
		try:
			outcar_file = open(directory_name+"/OUTCAR").readlines()
			for line in outcar_file:
				columns = line.split()
				number_columns = len(columns)
				if (number_columns > 4) and (columns[2] == "TOTEN"):
					total_energy = float(columns[4])
		except:
			oszicar_file = open(directory_name+"/OSZICAR").readlines()
			for line in oszicar_file:
				if "F=" in line:
					total_energy = float(line.split()[2])
			pass
		
		self.defects_data[compound_name][defect_name]["charge"][charge_state]["Energy"] = total_energy
		self.defects_data[compound_name][defect_name]["charge"][charge_state]["ECorr"] = 0.0
	
	
	####################################################################################################################
	########################################### Find Stoichiometry of Defect ###########################################
	####################################################################################################################
	
	def Get_Basic_Defect_Info(self, compound_name, defect_name, possible_defect_site_list):
		
		# Check if defect is extrinsic
		if (defect_name.split("_")[-1] in possible_defect_site_list) and (defect_name.split("_")[0] not in possible_defect_site_list):
			self.defects_data[compound_name][defect_name]["Extrinsic"] = "Yes"
		else:
			self.defects_data[compound_name][defect_name]["Extrinsic"] = "No"
		
		# Obtain "stoichiometry" of defect
		for specie in possible_defect_site_list:
			if specie == "V":
				if defect_name.split("_")[0] == specie:
					self.defects_data[compound_name][defect_name]["n_"+defect_name.split("_")[-1]] = -1
			elif specie not in defect_name.split("_"):
				self.defects_data[compound_name][defect_name]["n_"+specie] = 0
			else:
				if defect_name.split("_")[0] == specie:
					self.defects_data[compound_name][defect_name]["n_"+specie] = 1
				elif defect_name.split("_")[-1] == specie:
					self.defects_data[compound_name][defect_name]["n_"+specie] = -1
	
	
	####################################################################################################################
	########################################## Check Legitimacy of User Input ##########################################
	####################################################################################################################
	
	def Run_Checks(self, compound_name, directory_name):
		
		# Keep track of legitimacy
		legit = True
		
		# Check if the charge state is indeed a directory
		if not os.path.isdir(directory_name):
			print("The directory '"+directory_name+"' cannot be found. Skipping...")
			legit = False
		
		# Check if the Compounds_Tracker.json file exists
		try:
			open("Compounds_Tracker.json")
		except:
			print("The Compounds_Tracker.json file does not exist. Make sure you import the compounds data first. Exiting...")
			legit = False
			pass
		
		# If compound does not exist in Compounds_Tracker.json, then don't import
		if compound_name not in json.load(open("Compounds_Tracker.json"))["Compounds"].keys():
			print("The compound '"+compound_name+"' does not exist in Compounds_Tracker.json. Skipping...")
			legit = False
		
		return legit
	
	
	####################################################################################################################
	############################################# Update Defects Database ##############################################
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
	############################################### Add DOS to Database ################################################
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
		
		# DOSCAR normalizer
		volume = float(open(filename).readlines()[1].split()[0])
		normalizing_constant = 1E24/volume # VASP outputs DOS values normalized by volume of unit cell in A^3, this needs to be un-normalized and in cm^3
		
		# Extract data
		for line in open(filename).readlines():
			
			#Skip lines that describe the partial DOS with respect to orbital projections
			if len(line.split()) != 3:
				continue
			
			try:
				dos_info[float(line.split()[0])] = float(line.split()[1]) * normalizing_constant
			except:
				dos_info[float(line.split()[0])] = 0.0	# Sometimes the DOS can be an extremely small number that VASP outputs e.g. "0.5E-111" as "0.5-111", which is not readable by python
		
		# Store data
		self.dos_data[compound_name] = dos_info
	
	
	####################################################################################################################
	############################################### Update DOS Database ################################################
	####################################################################################################################
	
	def Update_DOS_Database(self):
		
		try:
			copyfile("DOS_Tracker.json", ".vtandem/DOS_Tracker_Backup.json")
		except:
			pass
		with open("DOS_Tracker.json", "w") as jsonfile:
			json.dump(self.dos_data, jsonfile, indent=4, sort_keys=True)





