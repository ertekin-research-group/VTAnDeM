
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
		
		# Find total energy of element in OUTCAR
		if "OUTCAR" not in os.listdir(directory_name):
			print("WARNING: Cannot find OUTCAR file of "+element_name+". Exiting...")
			return False
		element_data["dft_total_energy"] = self.Get_Total_Energy(directory_name)
		
		# Calculate total energy per formula unit of compound
		element_data["mu0"] = element_data["dft_total_energy"] / element_data["formula_units"]
		
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
		
		# Determine number of each element in compound from compound_name
		elements_list = [ ''.join( [letter for letter in element_segment if not letter.isdigit()] ) for element_segment in re.findall("[A-Z][^A-Z]*", compound_name) ]
		for element in elements_list:
			if element not in self.elements:
				print("'"+element+"' is not a valid element. Exiting...")
				return False
			else:
				try:
					number_species_in_compound = float(compound_name.split(element)[-1][:3])
				except:
					try:
						number_species_in_compound = float(compound_name.split(element)[-1][:2])
					except:
						try:
							number_species_in_compound = float(compound_name.split(element)[-1][:1])
						except:
							number_species_in_compound = 1.0
							pass
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
		
		"""
		# Find volume of compound from lattice vectors in POSCAR/CONTCAR
		lattice_vector_x = np.asarray([float(i) for i in structure_file[2].split()])
		lattice_vector_y = np.asarray([float(i) for i in structure_file[3].split()])
		lattice_vector_z = np.asarray([float(i) for i in structure_file[4].split()])
		volume = np.dot(lattice_vector_x, np.cross(lattice_vector_y, lattice_vector_z))
		compound_data["volume"] = volume*1E-24
		"""
		
		# Find total energy of compound in OUTCAR
		if "OUTCAR" not in os.listdir(directory_name):
			print("WARNING: Cannot find OUTCAR file of '"+compound_name+"' in '"+directory_name+"'. Exiting...")
			return False
		#compound_data["total_energy"] = self.Get_Total_Energy(directory_name)
		compound_data["dft_total_energy"] = self.Get_Total_Energy(directory_name)
		
		"""
		# Calculate total energy per formula unit of compound
		compound_data["mu0"] = compound_data["total_energy"] / compound_data["formula_units"]
		
		# Calculate enthalpy of formation of compound
		enthalpy_tracker = compound_data["mu0"]
		for element in elements_list:
			enthalpy_tracker -= compound_data[element]*self.compounds_info["Elements"][element]["mu0"]
		compound_data["enthalpy"] = enthalpy_tracker
		"""
		
		"""
		# Find band gap and valence band maximum of compound in vasprun.xml
		if "vasprun.xml" not in os.listdir(directory_name):
			print("WARNING: Cannot find band gap nor valence band maximum (vasprun.xml file). Exiting...")
			return False
		vasprun = Vasprun(directory_name+"/vasprun.xml")
		(bandgap, cbm, vbm, is_direct) = vasprun.eigenvalue_band_properties
		compound_data["band_gap"] = bandgap
		compound_data["vbm"] = vbm
		"""
		
		# Store data
		self.compounds_info["Compounds"][compound_name] = compound_data
		
		return True
	
	
	
	def Get_Total_Energy(self, directory_name):
		
		outcar_file = open(directory_name+"/OUTCAR").readlines()
		for line in outcar_file:
			if "entropy" not in line:
				continue
			if line.split()[4] == "energy(sigma->0)":
				total_energy = float(line.split()[-1])
		
		return total_energy
	
	
	
	####################################################################################################################
	################################################# Update Database ##################################################
	####################################################################################################################
	
	def Update_Compounds_Database(self):
		
		try:
			copyfile("Compounds_Tracker.json", ".vtandem/Compounds_Tracker_Backup.json")
		except:
			pass
		
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
		
		# All elements
		self.elements = [el.symbol for el in elements]
		
		# Dictionary to hold all site multiplicities
		self.site_multiplicities = {}
		
		# List of possible defect sites (updated when Add_Defects function is invoked)
		self.possible_defect_site_list = ["i", "V"]
		
		# Recreate defect information as dictionary (if JSON exists), otherwise create new datasheet
		if "Defects_Tracker.json" in os.listdir(os.getcwd()):
			with open("Defects_Tracker.json") as DefectsTracker:
				self.defects_data = json.load(DefectsTracker)
		else:
			self.defects_data = {}
	
	
	####################################################################################################################
	####################################### Add Defects of Compound to Database ########################################
	####################################################################################################################
	
	def Add_Defects(self, compound_name, directory_name):
		
		# Run checks
		self.Run_Checks(compound_name, directory_name)
		
		# Create list of possible defects in compound
		defect_sites = [ ''.join( [letter for letter in element_segment if not letter.isdigit()] ) for element_segment in re.findall("[A-Z][^A-Z]*", compound_name) ]
		for site in defect_sites:
			self.possible_defect_site_list.append(site)
		
		# Get bulk data first (needs to run before importing defects in order to update site multiplicities)
		self.Add_Bulk_Info(compound_name=compound_name, bulk_folder=directory_name+"/Bulk")
		
		# Loop through each item in the given directory
		for directory in os.listdir(directory_name):
			
			# Check that the directory is a legitimate defect name
			if ("_" in directory) and (directory.split("_")[-1] in self.possible_defect_site_list) and (directory.split("_")[0] in self.elements):
				self.Add_Single_Defect(compound_name=compound_name, defect_name=directory, directory_name=directory_name+"/"+directory)
			elif (directory == "Bulk"):
				continue
			else:
				print("Defect '"+directory+"' in directory '"+directory_name+"' is neither 1) a legitimate defect name for compound '"+compound_name+"' nor 2) the 'Bulk' folder. Skipping...")
				continue
	
	
	####################################################################################################################
	############################################# Add an Individual Defect #############################################
	####################################################################################################################
	
	def Add_Single_Defect(self, compound_name, defect_name, directory_name):
		
		for directory in os.listdir(directory_name):
			
			# Check that the folder name is in the correct format
			try:
				float(directory.split("q")[-1])
			except:
				print("The name '"+directory+"' in '"+directory_name+"' is not in the correct format for charge. The correct format for the charge state of a defect is 'q#', where '#' is an integer representing the charge state. Skipping...")
				continue
			
			# Check that the data exists
			if "OUTCAR" not in os.listdir(directory_name+"/"+directory):
				print("WARNING: Cannot find OUTCAR file for defect '"+defect_name+"' with charge state '"+charge_state+"' in '"+directory_name+"/"+directory+"'. Skipping...")
				continue
			
			self.Add_Defect_Charge(compound_name=compound_name, defect_name=defect_name, charge_state=directory.split("q")[-1], directory_name=directory_name+"/"+directory)
	
	
	####################################################################################################################
	###################################### Add Charge State for Individual Defect ######################################
	####################################################################################################################
	
	def Add_Defect_Charge(self, compound_name, defect_name, charge_state, directory_name):
		
		# Get total energy of defect
		total_energy = self.Get_Total_Energy(directory_name)
		
		# Set up data for compound/defects in Defects_Tracker.json
		if compound_name not in self.defects_data.keys():
			self.defects_data[compound_name] = {}
		if defect_name not in self.defects_data[compound_name].keys():
			self.defects_data[compound_name][defect_name] = {}
			self.Get_Basic_Defect_Info(compound_name, defect_name)
		if "charge" not in self.defects_data[compound_name][defect_name].keys():
			self.defects_data[compound_name][defect_name]["charge"] = {}
		if charge_state not in self.defects_data[compound_name][defect_name]["charge"].keys():
			self.defects_data[compound_name][defect_name]["charge"][charge_state] = {}
		
		# Update total energy of defect
		self.defects_data[compound_name][defect_name]["charge"][charge_state]["Energy"] = total_energy
		self.defects_data[compound_name][defect_name]["charge"][charge_state]["ECorr"] = 0.0
	
	
	####################################################################################################################
	########################################### Find Stoichiometry of Defect ###########################################
	####################################################################################################################
	
	def Get_Basic_Defect_Info(self, compound_name, defect_name):
		
		# Check if defect is extrinsic
		if (defect_name.split("_")[-1] in self.possible_defect_site_list) and (defect_name.split("_")[0] not in self.possible_defect_site_list):
			self.defects_data[compound_name][defect_name]["Extrinsic"] = "Yes"
		else:
			self.defects_data[compound_name][defect_name]["Extrinsic"] = "No"
		
		# Obtain "stoichiometry" of defect
		for specie in self.possible_defect_site_list:
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
		
		# Record site multiplicity (from bulk POSCAR/CONTCAR file)
		self.defects_data[compound_name][defect_name]["site_multiplicity"] = self.site_multiplicities[ defect_name.split("_")[-1] ]
	
	
	####################################################################################################################
	########################################### Add Supercell Energy of Bulk ###########################################
	####################################################################################################################
	
	def Add_Bulk_Info(self, compound_name, bulk_folder):
		
		# Check if compound is in Defects_Tracker.json
		if compound_name not in self.defects_data.keys():
			self.defects_data[compound_name] = {}
		
		# Open POSCAR/CONTCAR file
		try:
			structure_file = open(bulk_folder+"/CONTCAR").readlines()
		except:
			structure_file = open(bulk_folder+"/POSCAR").readlines()
			pass
		
		# Initialize
		self.defects_data[compound_name]["Bulk"] = {}
		
		# Update atom counts
		atom_types = structure_file[5].split()
		atom_counts = structure_file[6].split()
		for type, count in zip(atom_types, atom_counts):
			self.defects_data[compound_name]["Bulk"]["dft_"+type] = float(count)
		
		# Update total energy of bulk
		total_energy = self.Get_Total_Energy(bulk_folder)
		self.defects_data[compound_name]["Bulk"]["dft_BulkEnergy"] = total_energy
		
		# Find band gap and valence band maximum of compound in vasprun.xml
		vasprun = Vasprun(bulk_folder+"/vasprun.xml")
		(bandgap, cbm, vbm, is_direct) = vasprun.eigenvalue_band_properties
		self.defects_data[compound_name]["Bulk"]["BandGap"] = bandgap
		self.defects_data[compound_name]["Bulk"]["VBM"] = vbm
		
		# Find volume of compound from lattice vectors in POSCAR/CONTCAR
		lattice_vector_x = np.asarray([float(i) for i in structure_file[2].split()])
		lattice_vector_y = np.asarray([float(i) for i in structure_file[3].split()])
		lattice_vector_z = np.asarray([float(i) for i in structure_file[4].split()])
		volume = np.dot(lattice_vector_x, np.cross(lattice_vector_y, lattice_vector_z))
		self.defects_data[compound_name]["Bulk"]["Volume"] = volume*1E-24
		
		# Update site multiplicities
		self.Update_Site_Multiplicities(compound_name, bulk_folder)
	
	
	
	####################################################################################################################
	################################################# Get Total Energy #################################################
	####################################################################################################################
	
	def Get_Total_Energy(self, directory_name):
		
		outcar_file = open(directory_name+"/OUTCAR").readlines()
		for line in outcar_file:
			if "entropy" not in line:
				continue
			if line.split()[4] == "energy(sigma->0)":
				total_energy = float(line.split()[-1])
		
		return total_energy
	
	
	####################################################################################################################
	######################################### Update Defect Site Multiplicities ########################################
	####################################################################################################################
	
	def Update_Site_Multiplicities(self, compound_name, bulk_folder):
		
		# Open file
		try:
			structure_file = open(bulk_folder+"/POSCAR").readlines()
		except:
			structure_file = open(bulk_folder+"/CONTCAR").readlines()
		
		# Get atom names and number of each atom in structure
		atom_types = structure_file[5].split()
		number_atoms = structure_file[6].split()
		
		# Check that POSCAR/CONTCAR contains the correct atom types
		for atom_type in atom_types:
			if atom_type not in self.possible_defect_site_list:
				raise Exception("Unexpected atom found in POSCAR/CONTCAR file for '"+compound_name+"'. Exiting...")
		
		# Update site multiplicities dictionary
		for atom_type, natoms in zip(atom_types, number_atoms):
			self.site_multiplicities[atom_type] = float(natoms)
		self.site_multiplicities["i"] = 0.0
	
	
	
	####################################################################################################################
	######################################### Record Energy Correction Values ##########################################
	####################################################################################################################
	
	def Add_Energy_Corrections(self, compound_name, csv_filename):
		
		# Check if provided file name is an actual file
		if not os.path.isfile(csv_filename):
			raise Exception("The provided file '"+csv_filename+"' cannot be found/imported. Exiting...")
		
		# Check to see if the compound is actually in the database
		if compound_name not in self.defects_data.keys():
			raise Exception("The compound '"+compound_name+"' has not been imported yet. Exiting...")
		
		# Loop through lines of CSV file
		for line in open(csv_filename).readlines():
			line = line.replace("\n", "")
			data = line.split(",")
			
			# Check to see if the line indeed has 4 entries
			if len(data) != 4:
				print("The line '"+line+"' does not have 4 entries and therefore may not be in the correct format [Compound, Defect, Charge, ECorr]. Skipping...")
				continue
			
			# Check to see if the line is for the compound
			if data[0] != compound_name:
				print("The compound in line '"+line+"' is not '"+compound_name+"'. Skipping...")
				continue
			
			# Check to see if the defect has been imported
			if data[1] not in self.defects_data[compound_name].keys():
				print("The defect in line '"+line+"' cannot be found for compound '"+compound_name+"'. Skipping...")
				continue
			
			# Check to see if the defect charge state has been imported
			if data[2] not in self.defects_data[compound_name][data[1]]["charge"].keys():
				print("The charge state in line '"+line+"' cannot be found for compound '"+compound_name+"'. Skipping...")
				continue
			
			# Update energy correction
			self.defects_data[compound_name][data[1]]["charge"][data[2]]["ECorr"] = float(data[3])
	
	
	
	####################################################################################################################
	########################################## Check Legitimacy of User Input ##########################################
	####################################################################################################################
	
	def Run_Checks(self, compound_name, directory_name):
		
		# Check if the directory is indeed a directory
		if not os.path.isdir(directory_name):
			raise Exception("The directory '"+directory_name+"' cannot be found. Exiting...")
		
		# Check if the Compounds_Tracker.json file exists
		if not os.path.exists("Compounds_Tracker.json"):
			raise Exception("The Compounds_Tracker.json file does not exist. Make sure you import the compounds data first. Exiting...")
		
		"""
		# Check if compound exists in Compounds_Tracker.json
		if compound_name not in json.load(open("Compounds_Tracker.json"))["Compounds"].keys():
			raise Exception("The compound '"+compound_name+"' does not exist in Compounds_Tracker.json. Exiting...")
		"""
		# Check if elements in compound exist in Compounds_Tracker.json
		elements_in_compound = [ ''.join( [letter for letter in element_segment if not letter.isdigit()] ) for element_segment in re.findall("[A-Z][^A-Z]*", compound_name) ]
		elements_in_database = json.load(open("Compounds_Tracker.json"))["Elements"].keys()
		for element in elements_in_compound:
			if element not in elements_in_database:
				raise Exception("The element '"+element+"' of compound '"+compound_name+"' is not in Compounds_Tracker.json. Exiting...")
		
		# Check if directory contains the 'Bulk' folder
		if "Bulk" not in os.listdir(directory_name):
			raise Exception("A folder named 'Bulk' must exist in '"+directory_name+"'. Exiting...")
		
		# Check if POSCAR/CONTCAR, OUTCAR, and vasprun.xml are in the 'Bulk' folder
		if ("POSCAR" not in os.listdir(directory_name+"/Bulk")) and ("CONTCAR" not in os.listdir(directory_name+"/Bulk")):
			raise Exception("POSCAR/CONTCAR missing from 'Bulk' folder. Exiting...")
		if "OUTCAR" not in os.listdir(directory_name+"/Bulk"):
			raise Exception("OUTCAR missing from 'Bulk' folder. Exiting...")
		if "vasprun.xml" not in os.listdir(directory_name+"/Bulk"):
			raise Exception("vasprun.xml missing from 'Bulk' folder. Exiting...")
	
	
	
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
	
	def Add_DOS(self, compound_name, directory_name):
		
		# Check if the directory name is legitimate
		if not os.path.isdir(directory_name):
			raise Exception("The directory '"+directory_name+"' cannot be found. Exiting...")
		
		# Check that both the DOSCAR file and either POSCAR and/or CONTCAR are in the folder
		if ("DOSCAR" not in os.listdir(directory_name)) and ( ("POSCAR" not in os.listdir(directory_name)) and ("CONTCAR" not in os.listdir(directory_name)) ):
			raise Exception("Either the DOSCAR and/or the POSCAR/CONTCAR file is missing from directory '"+directory_name+"'. Exiting...")
		
		# If compound does not exists in Compounds_Tracker.json, then don't import
		if compound_name not in json.load(open("Compounds_Tracker.json"))["Compounds"].keys():
			raise Exception("The compound '"+compound_name+"' does not exist in Compounds_Tracker.json. Skipping...")
		
		# Check if the compound already exists in the database
		if compound_name in self.dos_data.keys():
			print("The compound '"+compound_name+"' is already in the DOS database (DOS_Tracker.json). The imported data will replace the old data.")
		
		# Open files
		if "POSCAR" in os.listdir(directory_name):
			lattice_vectors_str = open(directory_name+"/POSCAR").readlines()[2:5]
		elif "CONTCAR" in os.listdir(directory_name):
			lattice_vectors_str = open(directory_name+"/CONTCAR").readlines()[2:5]
		dos_file = open(directory_name+"/DOSCAR").readlines()
		
		# Initialize
		dos_info = {}
		
		# Volume to normalize DOS
		lat_a = np.array([float(x) for x in lattice_vectors_str[0].split()])
		lat_b = np.array([float(x) for x in lattice_vectors_str[1].split()])
		lat_c = np.array([float(x) for x in lattice_vectors_str[2].split()])
		#volume = float(open(filename).readlines()[1].split()[0])*1E-24	# CHECK LATER (VASP may be reporting the wrong volume)
		volume = np.dot( lat_a, np.cross(lat_b, lat_c) ) * 1E-24
		
		# Fermi energy
		fermi_energy = float( dos_file[5].split()[-2] )
		
		# Extract data
		for line in dos_file:
			
			# Skip lines that don't have the total DOS
			if len(line.split()) != 3:
				continue
			try:
				float(line.split()[0])
				float(line.split()[1])
				float(line.split()[2])
			except:
				continue
			
			try:
				dos_info[float(line.split()[0])-fermi_energy] = float(line.split()[1])
			except:
				dos_info[float(line.split()[0])-fermi_energy] = 0.0	# Sometimes the DOS can be an extremely small number that VASP outputs e.g. "0.5E-111" as "0.5-111", which is not readable by python
		
		# Store data
		self.dos_data[compound_name] = {"Volume": volume, "DOS": dos_info}
	
	
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





