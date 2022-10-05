
###############################################################
########### Python script to import data to VTAnDeM ###########
###############################################################

### IMPORTANT NOTE:
#	Make sure you create a VTAnDeM project first before
#	running this script.
#		vtandem --new


from vtandem.dft import import_dft

comp_import_obj = import_dft.Compounds_Import()
comp_import_obj.Add_Element("Mg", "PhaseStability/Bulk_Mg")
comp_import_obj.Add_Element("Si", "PhaseStability/Bulk_Si")
comp_import_obj.Add_Element("Sn", "PhaseStability/Bulk_Sn")
comp_import_obj.Add_Compound("Mg2Sn", "PhaseStability/Mg2Sn")
comp_import_obj.Update_Compounds_Database()

def_import_obj = import_dft.Defects_Import()
def_import_obj.Add_Defects("Mg2Si", "Mg2Si_Defects")
def_import_obj.Add_Energy_Corrections("Mg2Si", "EnergyCorrections_Mg2Si.csv")
def_import_obj.Update_Defects_Database()

dos_import_obj = import_dft.DOS_Import()
dos_import_obj.Add_DOS("Mg2Si", "DensityOfStates/DOSCAR")
dos_import_obj.Update_DOS_Database()

