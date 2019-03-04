
import numpy as np
import os
import sys
import pandas as pd
import Compounds_DataExtract as CDE

current_directory = os.getcwd()

data_directory = sys.argv[1]

# Check if directory exists
try:
	os.chdir(data_directory)
except:
	print("Cannot find directory... Aborting...")
	sys.exit()

# Find compounds in directory
compounds_list = []
for directory_name in os.listdir(data_directory):
	if os.path.isdir(directory_name):
		compounds_list.append(directory_name)

compounds_info = CDE.Phase_Stability_Compounds_Information()	# Access class

# Loop through listed species
try:
	for species in sys.argv[2:]:
		compounds_info.Add_Species(species)
except:
	print("No species listed... Aborting...")
	sys.exit()

# Obtain all DFT data
compounds_info.Add_Compounds(compounds_list)
compounds_info.Update_All()

# Output Compounds_Tracker.csv
os.chdir(current_directory)
compounds_info.Upload_Compounds_Info()

# Print data
compounds_info.Print_All_Info()




