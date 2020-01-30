
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'
__author__ = 'Michael_Lidia_Jiaxing_Elif'

import json
import periodictable

elements = []
for element in periodictable.elements:
	elements.append(element.symbol)

###############################################################################################
############################# Obtain DFT Data of Compounds ####################################
###############################################################################################

def Obtain_Compounds_Data(elements_list = None):		# For the phase stability diagram
	
	# Keep track of DFT data of all compounds in the analysis
	compounds_info = {}
	
	with open("Compounds_Tracker.json") as CompoundsTracker:
		compounds_data = json.load(CompoundsTracker)
	
	# Loop through elements in database
	for element in compounds_data["Elements"].keys():
		# Include information in compounds_info
		compounds_info[element] = compounds_data["Elements"][element]
	
	# Loop through compounds in database
	for compound in compounds_data["Compounds"].keys():
		# Check if compound is competing
		if set(compounds_data["Compounds"][compound]["elements_list"]).issubset(elements_list):
			compounds_info[compound] = compounds_data["Compounds"][compound]
	
	return compounds_info


def Obtain_Defects_Data():	# For the defects diagram
	
	with open("Defects_Tracker.json") as DefectsTracker:
		defects_info = json.load(DefectsTracker)
	
	return defects_info


def Obtain_DOS_Data():	# For the carrier concentration and equilibrium Fermi energy
	
	with open("DOS_Tracker.json") as DOSTracker:
		dos_data = json.load(DOSTracker)
	
	return dos_data







