
import os, sys
import json
import click

default_values = {
	"import_phase_stability": 				("None", "./"), \
	"import_defects": 						("None", "./"), \
	"import_defect_energy_corrections":		("None", "./"), \
	"import_dos": 							("None", "./"), \
	"new": 									False, \
	"open":									False, \
	"visualize": 							False
}

import_element_help = 						"Import data for an element for generating the phase stability diagram (see above [1])."
import_compound_help = 						"Import data for a compound for generating the phase stability diagram (see above [2])."
import_defects_help = 						"Import defects data (see above [3])."
import_defect_energy_corrections_help =		"Import defect energy corrections (see above [4])."
import_dos_help = 							"Import density of states data (see above [5])."

@click.command()
@click.option("--import_element", "-e", default=default_values["import_phase_stability"], type=(str, click.Path(exists=True)), help=import_element_help)
@click.option("--import_compound", default=default_values["import_phase_stability"], type=(str, click.Path(exists=True)), help=import_compound_help)
@click.option("--import_defects", default=default_values["import_defects"], type=(str, click.Path(exists=True)), help=import_defects_help)
@click.option("--import_defect_energy_corrections", default=default_values["import_defect_energy_corrections"], type=(str, click.Path(exists=True)), help=import_defect_energy_corrections_help)
@click.option("--import_dos", default=default_values["import_dos"], type=(str, click.Path(exists=True)), help=import_dos_help)
@click.option("--new", "-n", is_flag=True, help="Initializes a new VTAnDeM project.")
@click.option("--open", "-o", is_flag=True, help="Open VTAnDeM import data dialog.")
@click.option("--visualize", "-v", is_flag=True, help="Open material selection dialog.")

def vtandem(import_element, import_compound, import_defects, import_defect_energy_corrections, import_dos, new, open, visualize):
	""" 
	\b
	======================================================================
	\b
	The Visualization Toolkit for Analyzing Defects in Materials (VTAnDeM)
	\b
	======================================================================
	\b
	__        __ ____________   ___             _______           ___      ___
	\ \      / / |____  ____|  /   \           |   __  )         |   \    /   |
	 \ \    / /      |  |     / ( ) \    _____ |  |  )  )   ___  | |\ \  / /| |
	  \ \  / /       |  |    /  ___  \  |  _  ||  |   )  ) / _ \ | | \ \/ / | |
	   \ \/ /        |  |   /  /   \  \ | | | ||  |__)  ) (  __/ | |  \__/  | |
	    \__/         |__|  /__/     \__\|_| |_||_______)   \___) |_|        |_|
	
	\b
	[1] Importing Elements/Compounds
	When importing elements, use the --import_element option. The argument
	<TEXT PATH> should be in the form:
	    'Element_Name /path/to/data/folder'
	'Element_Name' is the case-sensitive name of the element (e.g. Cu).
	/path/to/data/folder should have the file structure:
	
	\b
	  /path/to/data/folder
	    |-- OUTCAR (or OSZICAR)
	    |-- CONTCAR (or POSCAR)
	
	\b
	\b
	[2] Importing Compounds
	When importing compounds, use the --import_compound option. The argument
	<TEXT PATH> should be in the form:
		'Compound_Name /path/to/data/folder'
	'Compound_Name' is the case-sensitive name of the compound (e.g. Cu2HgGeTe4).
	/path/to/data/folder should have the file structure:
	
	\b
	  /path/to/data/folder
	    |-- OUTCAR (or OSZICAR)
	    |-- CONTCAR (or POSCAR)
	
	\b
	\b
	[3] Importing Defects Information
	When importing defects using the --import_defects option, the argument
	<TEXT PATH> should be given in the form:
	    'Compound_Name /path/to/data/folder'
	'Compound_Name' is case-sensitive (e.g. Cu2HgGeTe4).
	/path/to/defects/data/folder should have the file structure:
	
	\b
	  /path/to/defects/data/folder
	    |-- Bulk
	          |-- OUTCAR
	          |-- CONTCAR
	          |-- vasprun.xml
	    |-- Defect1
	          |-- q0
	              |-- OUTCAR
	          |-- q+1
	              |-- OUTCAR
	          |-- ...
	    |-- Defect2
	          |-- q0
	              |-- OUTCAR
	          |-- q-1
	              |-- OUTCAR
	          |-- ...
	    |-- ...
	
	\b
	\b
	[4] Importing Defect Energy Corrections
	When import defect energy corrections, use the --import_defect_energy_corrections
	option. The argument <TEXT PATH> should be given in the form:
	    'Compound_Name /path/to/energy/corrections/csv'
	'Compound_Name' is case-sensitive (e.g. Cu2HgGeTe4).
	/path/to/energy/corrections/csv is a comma-separated values file, where
	    each line has the form:
	
	\b
	  Compound_Name, Defect, Charge, ECorr
	
	\b
	\b
	[5] Importing Density of States
	When importing the DOS using the --import_dos option, the argument  <TEXT
	PATH> should be in the form:
	    'Compound_Name /path/to/DOSCAR'
	'Compound_Name' is case-sensitive (e.g. Cu2HgGeTe4).
	/path/to/DOSCAR is the name of the DOSCAR file containing the DOS info.
	\b
	
	"""
	
	
	print(
"""
__        __ ____________   ___             _______           ___      ___
\ \      / / |____  ____|  /   \           |   __  )         |   \    /   |
 \ \    / /      |  |     / ( ) \    _____ |  |  )  )   ___  | |\ \  / /| |
  \ \  / /       |  |    /  ___  \  |  _  ||  |   )  ) / _ \ | | \ \/ / | |
   \ \/ /        |  |   /  /   \  \ | | | ||  |__)  ) (  __/ | |  \__/  | |
    \__/         |__|  /__/     \__\|_| |_||_______)   \___) |_|        |_|

"""
	)
	
	# Create new VTAnDeM project
	if new:
		if not os.path.isdir(".vtandem"):
			Make_New_VTAnDeM_Project()
		else:
			print("VTAnDeM project already exists!")
		return
	
	# Import element data to Compounds_Tracker.json
	if (import_element[0] != default_values["import_phase_stability"][0]) and (import_element[1] != default_values["import_phase_stability"][1]):
		if not Check_VTAnDeM_Project():
			sys.exit("Cannot find VTAnDeM project. Exiting...")
		from vtandem.dft.import_dft import Compounds_Import
		element_import_object = Compounds_Import()
		"""
		import_element_success = element_import_object.Add_Element(import_element[0], import_element[1])
		if import_element_success:
			element_import_object.Update_Compounds_Database()
			print("Imported element '"+import_element[0]+"' from the folder '"+import_element[1]+"' successfully!")
		"""
		element_import_object.Add_Element(import_element[0], import_element[1])
		element_import_object.Update_Compounds_Database()
		print("Imported element '"+import_element[0]+"' from the folder '"+import_element[1]+"' successfully!")

	# Import compound data to Compounds_Tracker.json
	if (import_compound[0] != default_values["import_phase_stability"][0]) and (import_compound[1] != default_values["import_phase_stability"][1]):
		if not Check_VTAnDeM_Project():
			sys.exit("Cannot find VTAnDeM project. Exiting...")
		from vtandem.dft.import_dft import Compounds_Import
		compound_import_object = Compounds_Import()
		"""
		import_compound_success = compound_import_object.Add_Compound(import_compound[0], import_compound[1])
		if import_compound_success:
			compound_import_object.Update_Compounds_Database()
			print("Imported compound '"+import_compound[0]+"' from the folder '"+import_compound[1]+"' successfully!")
		"""
		compound_import_object.Add_Compound(import_compound[0], import_compound[1])
		compound_import_object.Update_Compounds_Database()
		print("Imported compound '"+import_compound[0]+"' from the folder '"+import_compound[1]+"' successfully!")

	# Import defects data to Defects_Tracker.json
	if (import_defects[0] != default_values["import_defects"][0]) and (import_defects[1] != default_values["import_defects"][1]):
		if not Check_VTAnDeM_Project():
			sys.exit("Cannot find VTAnDeM project. Exiting...")
		from vtandem.dft.import_dft import Defects_Import
		defects_import_object = Defects_Import()
		defects_import_object.Add_Defects(import_defects[0], import_defects[1])
		defects_import_object.Update_Defects_Database()
		print("Imported defects of compound '"+import_defects[0]+"' from the folder '"+import_defects[1]+"' successfully!")
	
	# Import defect energy corrections (ECorr) to Defects_Tracker.json
	if (import_defect_energy_corrections[0] != default_values["import_defect_energy_corrections"][0]) and (import_defect_energy_corrections[1] != default_values["import_defect_energy_corrections"][1]):
		if not Check_VTAnDeM_Project():
			sys.exit("Cannot find VTAnDeM project. Exiting...")
		from vtandem.dft.import_dft import Defects_Import
		DECs_import_object = Defects_Import()
		DECs_import_object.Add_Energy_Corrections(import_defect_energy_corrections[0], import_defect_energy_corrections[1])
		DECs_import_object.Update_Defects_Database()
		print("Imported defect energy corrections of compound '"+import_defect_energy_corrections[0]+"' from the file '"+import_defect_energy_corrections[1]+"' successfully!")
	
	# Import density of states data to DOS_Tracker.json
	if (import_dos[0] != default_values["import_dos"][0]) and (import_dos[1] != default_values["import_dos"][1]):
		if not Check_VTAnDeM_Project():
			sys.exit("Cannot find VTAnDeM project. Exiting...")
		from vtandem.dft.import_dft import DOS_Import
		dos_import_object = DOS_Import()
		dos_import_object.Add_DOS(import_dos[0], import_dos[1])
		dos_import_object.Update_DOS_Database()
		print("Imported density of states of compound '"+import_dos[0]+"' from the folder '"+import_dos[1]+"' successfully!")
	
	# Open VTAnDeM import data dialog
	if open:
		if not Check_VTAnDeM_Project():
			sys.exit("Cannot find VTAnDeM project. Exiting...")
		print("Opening VTAnDeM welcome window...")
		from vtandem.gui_windows import Open_Welcome_VTAnDeM_Window
		Open_Welcome_VTAnDeM_Window()
	
	if visualize:
		if not Check_VTAnDeM_Project():
			sys.exit("Cannot find VTAnDeM project. Exiting...")
		print("Opening VTAnDeM material selection dialog...")
		from vtandem.gui_windows import Open_Material_Selection_Window
		Open_Material_Selection_Window()
	
	if (import_element==default_values["import_phase_stability"]) and (import_compound==default_values["import_phase_stability"]) and (import_defects==default_values["import_defects"]) and (import_dos==default_values["import_dos"]) and (new==default_values["new"]) and (open==default_values["open"]) and (visualize==default_values["visualize"]):
		print("No options declared... Type 'vtandem --help' to show options.")



def Make_New_VTAnDeM_Project():
	
	# Make VTAnDeM project folder
	current_directory = os.getcwd()
	os.mkdir(current_directory+"/.vtandem")
	
	# Create data files in VTAnDeM project folder
	with open(current_directory+"/Compounds_Tracker.json", "w") as compounds_tracker_json:
		json.dump({"Compounds": {}, "Elements": {}}, compounds_tracker_json, indent=4, sort_keys=True)
	with open(current_directory+"/Defects_Tracker.json", "w") as defects_tracker_json:
		json.dump({}, defects_tracker_json, indent=4, sort_keys=True)
	with open(current_directory+"/DOS_Tracker.json", "w") as dos_tracker_json:
		json.dump({}, dos_tracker_json, indent=4, sort_keys=True)
	
	# Print status
	print("Created new VTAnDeM project!")



def Check_VTAnDeM_Project():
	
	# Check to see that the directory is a legitimate VTAnDeM project
	if os.path.isfile("Compounds_Tracker.json") and os.path.isfile("Defects_Tracker.json") and os.path.isfile("DOS_Tracker.json") and os.path.isdir(".vtandem"):
		return True
	else:
		return False


if __name__ == "__main__":
	vtandem()







