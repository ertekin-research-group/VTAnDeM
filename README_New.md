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


VTAnDeM Data Input:
-------------------
All necessary data for phase stability is provided in the Example_Data/Materials_Data.tar.gz file on Box.
The data for defects is provided as Example_Data/Defects_Data.tar.gz.


Examples:
------
Example VTAnDeM projects can be found for the following materials:
1. Mg<sub>2</sub>Si
2. Hg<sub>2</sub>GeTe<sub>4</sub>
3. Cu<sub>2</sub>HgGeTe<sub>4</sub>

