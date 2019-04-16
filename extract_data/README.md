Obtaining Necessary DFT Data for VTAnDeM
========================================

Output from VASP, specifically the OUTCAR and POSCAR files, are necessary for VTAnDeM.
The `Generate_Compounds_Data.py` script easily extracts the necessary DFT output (compatible only with VASP outputs for now) and stores them into a VTAnDeM-readable csv file.


Directory Layout
----------------

The DFT data for all compounds should be located in one folder. The main quaternary material and the individual constituents of the quaternary should be in this folder at the very least (e.g. if the quaternary compound is Cu2HgGeTe4, then DFT data for Cu, Hg, Ge, and Te should also be provided).
The data should be divided by compound into different folders. Each folder should be named by the compound name. Each folder should contain the POSCAR and OUTCAR files. For example, in the case of Cu2HgGeTe4, the compounds are Cu, Hg, Ge, Te, CuTe, GeTe, HgTe, and Hg2GeTe4.


Syntax
------

`Generate_Compounds_Data.py` extracts the necessary quantities (e.g. stoichiometries, total energies, etc.) from the provided compounds. The syntax is:

	python3 Generate_Compounds_Data.py FOLDER_NAME SPECIES1 SPECIES2 ...

where `FOLDER_NAME` is the path to the folder containing all the compounds, and `SPECIES#` is a list of the species contained in the quaternary compound.
