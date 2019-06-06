
__author__ = 'Michael_Lidia_Jiaxing_Benita_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'

import pandas as pd

###############################################################################################
############################# Obtain DFT Data of Compounds ####################################
###############################################################################################

def Obtain_DFT_Data(parent = None):		# For the phase stability diagram
	
	if parent == None:
		print("WARNING: No parent found, output empty set.")
		return {}
	
	# Keep track of DFT data of all compounds in the analysis
	compounds_info = {}
	
	# Import compounds info from CSV created by the script "Compounds_DataExtract.py"
	compounds_data = pd.read_csv("Compounds_Tracker.csv", header=0, index_col=0)
	header_information = list(compounds_data)
	
	for compound, compound_info in compounds_data.iterrows():						# Loop through compounds and their rows.
		compounds_info[compound] = {}
		for header_specie, compound_data_info in zip(header_information, list(compound_info)):
			compounds_info[compound][header_specie] = compound_data_info	# Record number of each unique element in the compound in the list.
		for specie in parent.species_list:
			if specie not in compounds_info[compound].keys():
				compounds_info[compound][specie] = 0
	
	return compounds_info


def Obtain_Defects_Data(parent = None):	# For the defects diagram
	
	if parent == None:
		print("WARNING: No parent found, output empty set.")
		return {}
	
	# Keep track of defects data of the main quaternary compound
	defects_info = pd.read_csv("Defects_Tracker.csv")
	
	return defects_info








