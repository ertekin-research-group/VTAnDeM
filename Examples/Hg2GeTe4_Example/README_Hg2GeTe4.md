# Example: Defects in Hg<sub>2</sub>GeTe<sub>4</sub>

The defect calculations for Hg<sub>2</sub>GeTe<sub>4</sub> in this tutorial have been used in the following publication:

>> J. Qu, et al.,
   [J. Mater. Chem. A **9**, 26189 (2021)](https://doi.org/10.1039/D1TA07410E)


## Steps

1. Create new VTAnDeM project
	```
	vtandem -n
	```

2. Import chemical potentials
	```
	vtandem --import_element Hg PhaseStability/Hg
	vtandem --import_element Ge PhaseStability/Ge
	vtandem --import_element Te PhaseStability/Te
	```

3. Import competing phases
	```
	vtandem --import_compound GeTe PhaseStability/GeTe
	vtandem --import_compound HgTe PhaseStability/HgTe
	```

4. Import defects
	```
	vtandem --import_defects Hg2GeTe4 Hg2GeTe4_Defects
	```

5. Import defect energy corrections
	```
	vtandem --import_defect_energy_corrections Hg2GeTe4 HGT_EnergyCorrections.csv
	```

6. Visualize
	```
	vtandem -v
	```


**NOTE 1**: All defect calculations were performed on 2 x 2 x 2 supercells of Hg<sub>2</sub>GeTe<sub>4</sub>. However, since the bulk structure (in the "Cu2HgGeTe4_Defects/Bulk" folder) is the primitive cell, we need to adjust the volume, bulk energy, number of atoms, and site multiplicity.

**NOTE 2**: Check "VTAnDeM_Project_Reference" folder for reference.


