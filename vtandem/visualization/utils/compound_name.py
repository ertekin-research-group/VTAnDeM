
__author__ = 'Michael_Lidia_Jiaxing_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'

import re

###############################################################################################
########################## Rewrite Compound Name Latex-Style ##################################
###############################################################################################

def Compound_Name_Formal(compound_name, type: str):
	
	# Check if 'type' variable is either 'unicode' or 'latex'
	if (type != "unicode") and (type != "latex"):
		raise ValueError("'type' variable must be one of: 'unicode' or 'latex'")
	
	# Obtain the chemistry and stoichiometry of the compound
	compound_stoichiometry_list = re.findall(r'\d+|[A-Z][a-z]*', compound_name)

	# Create fancy compound name
	compound_name_formal = r""
	for i in compound_stoichiometry_list:
		if i.isnumeric():
			if i == "1":
				continue
			if type == "latex":
				compound_name_formal += "$_\mathrm{"+str(int(i))+"}$"
			elif type == "unicode":
				compound_name_formal += "<sub>"+str(int(i))+"</sub>"
		else:
			compound_name_formal += i

	return compound_name_formal

