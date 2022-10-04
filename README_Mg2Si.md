# Example: Defects in Mg<sub>2</sub>Si

The defect calculations for Mg<sub>2</sub>Si in this tutorial have been used in the following publications:

>> M.Y. Toriyama, M.K. Brod, and G.J. Snyder,
   [ChemNanoMat **8**, e202200222 (2022)](https://doi.org/10.1002/cnma.202200222)

>> R. Orenstein, et al.,
   [J. Mater. Chem. A **9**, 7208 (2021)](https://doi.org/10.1039/D1TA00115A)


## Steps

1. Create new VTAnDeM project

							vtandem -n

2. Import chemical potentials
						
	"vtandem --import_element Mg PhaseStability/Bulk_Mg"
	"vtandem --import_element Si PhaseStability/Bulk_Si"
	"vtandem --import_element Sn PhaseStability/Bulk_Sn"

3. Import competing phases

	"vtandem --import_compound Mg2Sn PhaseStability/Mg2Sn"

4. Import defects

	"vtandem --import_defects Mg2Si Defects"

5. Import defect energy corrections

	"vtandem --import_defect_energy_corrections Mg2Si EnergyCorrections_Mg2Si.csv"

6. Visualize

	"vtandem -v"


## Optional

7. Import DOS

	"vtandem --import_dos Mg2Si DensityOfStates/DOSCAR"


NOTE 1: Change the "site_multiplicity" field in Defects_Tracker.json for interstitials to 8, otherwise the carrier concentrations will be incorrect.

NOTE 2: Check "VTAnDeM_Project_Reference" folder for reference.

