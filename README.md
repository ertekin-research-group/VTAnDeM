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
- Click
	- `pip3 install click`
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
	- `git clone https://github.com/ertekin-research-group/VTAnDeM`
- Run `python3 setup.py install` in the downloaded VTAnDeM folder.


Examples:
---------
Example VTAnDeM projects can be found in the `Examples` folder for the following materials:
1. Mg<sub>2</sub>Si
2. Hg<sub>2</sub>GeTe<sub>4</sub>
3. Cu<sub>2</sub>HgGeTe<sub>4</sub>


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
			vtandem --import_element Cu PhaseStability/Cu
			vtandem --import_element Hg PhaseStability/Hg
			vtandem --import_element Ge PhaseStability/Ge
			vtandem --import_element Te PhaseStability/Te
			vtandem --import_compound CuTe PhaseStability/CuTe
			vtandem --import_compound GeTe PhaseStability/GeTe
			vtandem --import_compound HgTe PhaseStability/HgTe
			vtandem --import_compound Hg2GeTe4 PhaseStability/Hg2GeTe4
			vtandem --import_compound Cu2GeTe3 PhaseStability/Cu2GeTe3
			```
		- Defects Data:
			```
			vtandem --import_defects Cu2HgGeTe4 Cu2HgGeTe4_Defects
			```
		- Defect Energy Corrections:
			```
			vtandem --import_defect_energy_corrections Cu2HgGeTe4 CHGT_EnergyCorrections.csv
			```
		- DOS Data:
			```
			vtandem --import_dos Cu2HgGeTe4 DOSCAR
			```
	3. Import data from Python interface.
		- Phase Stability Data:
			```python
			from vtandem.dft import import_dft
			x = import_dft.Compounds_Import()
			x.Add_Element("Cu", "PhaseStability/Cu")
			x.Add_Element("Hg", "PhaseStability/Hg")
			x.Add_Element("Ge", "PhaseStability/Ge")
			x.Add_Element("Te", "PhaseStability/Te")
			x.Add_Compound("CuTe", "PhaseStability/CuTe")
			x.Add_Compound("GeTe", "PhaseStability/GeTe")
			x.Add_Compound("HgTe", "PhaseStability/HgTe")
			x.Add_Compound("Hg2GeTe4", "PhaseStability/Hg2GeTe4")
			x.Add_Compound("Cu2GeTe3", "PhaseStability/Cu2GeTe3")
			x.Update_Compounds_Database()
			```
		- Defects Data:
			```python
			from vtandem.dft import import_dft
			x = import_dft.Defects_Import()
			x.Add_Defects("Cu2HgGeTe4", "Cu2HgGeTe4_Defects")
			x.Add_Energy_Corrections("Cu2HgGeTe4", "CHGT_EnergyCorrections.csv")
			x.Update_Defects_Database()
			```
		- DOS Data:
			```python
			from vtandem.dft import import_dft
			x = import_dft.DOS_Import()
			x.Add_DOS("Cu2HgGeTe4", "DOSCAR")
			x.Update_DOS_Database()
			```
3. Open the VTAnDeM interface.
	```
	vtandem --visualize
	```
	(Note that this step is not necessary if you opened the VTAnDeM UI in step 2.i)


Help:
-----
Help with the necessary data structure can be found readily on the command terminal.

> vtandem --help

