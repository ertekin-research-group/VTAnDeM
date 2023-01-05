# Example: Defects in Y<sub>2</sub>Te<sub>3</sub>

The defect calculations for Y<sub>2</sub>Te<sub>3</sub> in this tutorial are sourced from the following publication:

>> M.Y. Toriyama, et al.,
   [ACS Appl. Mater. Interfaces **14**, 43517 (2022)](https://doi.org/10.1021/acsami.2c12112)


## Steps

1. Import defects (chemical potentials and competing phases have already been imported to Compounds_Tracker.json)
	```
	vtandem --import_defects Y2Te3 Y2Te3_Defects
	```

2. Import defect energy corrections
	```
	vtandem --import_defect_energy_corrections Y2Te3 Y2Te3_DefectEnergyCorrections.csv
	```

3. Visualize
	```
	vtandem -v
	```


**NOTE 1**: Due to band edge shifting from SOC and GW calculations, in order to reproduce the results in the publication above, please set "BandGap" to 1.779 and "VBM" to 2.817 in the Defects_Tracker.json file.

**NOTE 2**: Check "VTAnDeM_Project_Reference" folder for reference.


