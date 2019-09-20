The Visualization Toolkit for Analyzing Defects in Materials (VTAnDeM)
======================================================================

VTAnDeM is a post-processing plotting toolkit for DFT calculations of defects in materials.
The toolkit allows simultaneous visualization of interconnected thermodynamic and electronic properties of materials, including phase stability, defects, and carrier concentrations.


Python Version
--------------
python >= 3.5


Required Packages
-----------------
- Numpy >= 1.16
- Matplotlib >= 3.0
- LabelLines
	- `pip3 install matplotlib-label-lines`
- PeriodicTable
	- `pip3 install periodictable`
- Pymatgen == 2019.5.8
	- `pip3 install pymatgen==2019.5.8`
- PyQt5 == 5.11.3
	- `pip3 install PyQt5==5.11.3`
- PyPolyhedron (Courtesy of Dr. Pearu Peterson and Dr. Sunghyun Kim, https://github.com/frssp/PyPolyhedron)
	- Steps:
		1. `git clone https://github.com/frssp/PyPolyhedron`
		2. `cd /path/to/PyPolyhedron`
		3. `python3 setup.py install`


Installation
---------------
- Download all VTAnDeM files.
	- `git clone https://github.com/mathtoriyama/VTAnDeM`
- Run `python3 setup.py install` in the downloaded VTAnDeM folder.


VTAnDeM Data Input:
-------------------
All necessary data for phase stability is provided in the folder Ertekin_Research_Group/IndividualFolders/Elif-Michael-Jiaxing-Benita-Lidia/VTAnDeM/Materials_Data.tar.gz file on Box.
The data for defects is provided in the same folder under the name Defects_Data.tar.gz.


Usage:
------
VTAnDeM must be called from the command terminal. The steps are as follows:
1. Create a **VTAnDeM project** in a directory of your choice.
	```
	mkdir vtandem_project
	cd vtandem_project
	vtandem --new
	```
2. Import all of your DFT data. This can be done in one of three ways:
	1. Open the VTAnDeM user interface.
		```
		vtandem --open
		```
		(Help can be found in the **help buttons**, which are included in the import dialogs.)
	2. Import data from the command line.
		- Phase Stability Data:
			```
			vtandem --import_element Cu ~/Materials_Data/Cu
			vtandem --import_element Hg ~/Materials_Data/Hg
			vtandem --import_element Ge ~/Materials_Data/Ge
			vtandem --import_element Te ~/Materials_Data/Te
			vtandem --import_compound Cu2HgGeTe4 ~/Materials_Data/Cu2HgGeTe4
			```
		- Defects Data:
			```
			vtandem --import_defects Cu2HgGeTe4 ~/Defects_Data/Cu2HgGeTe4 2 2 2
			```
		- DOS Data:
			```
			vtandem --import_dos Cu2HgGeTe4 ~/Materials_Data/Cu2HgGeTe4/DOSCAR
			```
		- Help can be found with the `vtandem --help` command.
	3. Import data from Python interface.
		- Phase Stability Data:
			```python
			from vtandem.dft import import_dft
			x = import_dft.Compounds_Import()
			x.Add_Element("Cu", "~/Materials_Data/Cu")
			x.Add_Element("Hg", "~/Materials_Data/Hg")
			x.Add_Element("Ge", "~/Materials_Data/Ge")
			x.Add_Element("Te", "~/Materials_Data/Te")
			x.Add_Compound("Cu2HgGeTe4", "~/Materials_Data/Cu2HgGeTe4")
			x.Update_Compounds_Database()
			```
		- Defects Data:
			```python
			from vtandem.dft import import_dft
			x = import_dft.Defects_Import()
			x.Add_Defects("Cu2HgGeTe4", "~/Defects_Data/Cu2HgGeTe4", 2, 2, 2)
			x.Update_Defects_Database()
			```
		- DOS Data:
			```python
			from vtandem.dft import import_dft
			x = import_dft.DOS_Import()
			x.Add_DOS("Cu2HgGeTe4", "~/Materials_Data/Cu2HgGeTe4/DOSCAR")
			x.Update_DOS_Database()
			```
3. Open the VTAnDeM interface.
	```
	vtandem --visualize
	```
	(Note that this step is not necessary if you opened the VTAnDeM UI in step 2.1)

