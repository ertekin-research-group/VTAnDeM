
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'
__author__ = 'Michael_Lidia_Jiaxing_Elif'

import periodictable
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter, PDEntry
from pymatgen.core.composition import Composition
import re

from vtandem.visualization.utils.compound_name import Compound_Name_Formal


# All elements in the periodic table
all_elements = []
for element in periodictable.elements:
	all_elements.append(str(element))
		


#def Create_Compositional_PhaseDiagram(compounds_info, elements_list):
def Create_Compositional_PhaseDiagram(compounds_info, elements_list, main_compound_info, main_compound):
	
	# Copy main compound bulk info to compounds_info (phase stability analysis will be done using Bulk info in Defects_Tracker.json)
	main_compound_elements_list = []  # self.elements_list may contain a dopant atom; we only want elements in main compound
	for element in elements_list:
		if "dft_"+element in main_compound_info.keys():
			main_compound_elements_list.append(element)
	compounds_info[main_compound] = {}
	compounds_info[main_compound]["elements_list"] = main_compound_elements_list
	for element in main_compound_elements_list:
		compounds_info[main_compound][element] = main_compound_info["dft_"+element]
		compounds_info[main_compound]["dft_"+element] = main_compound_info["dft_"+element]
	compounds_info[main_compound]["number_species"] = main_compound_info["number_species"]
	compounds_info[main_compound]["formula_units"] = 1.0
	compounds_info[main_compound]["dft_total_energy"] = main_compound_info["dft_BulkEnergy"]


	# Record all entries for the phase diagram
	phasediagram_entries = []
	
	for compound in compounds_info.keys():

		# Disregard elements not included in main compound
		if (compound in all_elements) and (compound not in elements_list):
			continue
		
		# Get the compound's composition
		compound_composition = {}
		if compound in elements_list:	# Elemental material
			compound_composition[compound] = compounds_info[compound]["dft_"+compound]
		else:	# Compound material
			for element in compounds_info[compound]["elements_list"]:
				compound_composition[element] = compounds_info[compound]["dft_"+element]
		
		# Get the compound's total energy
		compound_total_energy = compounds_info[compound]["dft_total_energy"]

		# Record to list of entries
		phasediagram_entries.append(PDEntry(composition=Composition(compound_composition), energy=compound_total_energy, name=compound))

	# Calculate compositional phase diagram (using pymatgen)
	#	The output data structure is as follows:
	#		lines --> List of arrays, each array is 2x2 for ternary (3x3 for quaternary, etc.), column vector represents point on phase diagram.
	#					ex: array([ [0.3, 0.5], [1.0, 0.0] ]) is a line that goes from point [x=0.3, y=1.0] to point [x=0.5, y=0.0]
	#		labels --> Dictionary with point-PDEntry pairs.
	pmg_phasediagram = PhaseDiagram(phasediagram_entries)
	pmg_phasediagram_plot_object = PDPlotter(pmg_phasediagram)
	(lines, labels, unstable) = pmg_phasediagram_plot_object.pd_plot_data

	return pmg_phasediagram, lines, labels




"""
def Plot_Compositional_PhaseDiagram(axes, type, lines, labels, fontdict, label_stable=True):
	
	# Plot settings
	linewidth = 1
	markersize = 5
	
	# Plot compositional phase diagram
	count = 1
	newlabels = []
	if type == "ternary":
		for x, y in lines:
			axes.plot(x, y, "bo-", linewidth=linewidth, markeredgecolor="b", markerfacecolor="r", markersize=markersize)
		for coords in sorted(labels.keys()):
			entry = labels[coords]
			label = entry.name
			if label_stable:
				if len(entry.composition.elements) == 1:
					axes.text(coords[0], coords[1], label, fontdict=fontdict)
				else:
					axes.text(coords[0], coords[1], str(count), fontdict=fontdict)
					newlabels.append("{} : {}".format(count, label))
					count += 1
	elif type == "quaternary":
		for x, y, z in lines:
			axes.plot(x, y, z, "bo-", linewidth=linewidth, markeredgecolor="b", markerfacecolor="r", markersize=markersize)
		for coords in sorted(labels.keys()):
			entry = labels[coords]
			label = entry.name
			if label_stable:
				if len(entry.composition.elements) == 1:
					axes.text(coords[0], coords[1], coords[2], label, fontdict=fontdict)
				else:
					axes.text(coords[0], coords[1], coords[2], str(count), fontdict=fontdict)
					newlabels.append("{} : {}".format(count, label))
					count += 1
	
	# Remove Cartesian axes
	axes.axis("off")
	
	return newlabels
"""


def Plot_Compositional_PhaseDiagram(axes, type, lines, labels, fontdict):
	
	# Plot settings
	linewidth = 1
	markersize = 5
	
	# Plot compositional phase diagram
	count = 1
	newlabels = []
	if type == "ternary":
		for x, y in lines:
			axes.plot(x, y, "bo-", linewidth=linewidth, markeredgecolor="b", markerfacecolor="r", markersize=markersize)
		for coords in sorted(labels.keys()):
			entry = labels[coords]
			label = entry.name
			if len(entry.composition.elements) == 1:
				axes.text(coords[0], coords[1], label, fontdict=fontdict)
			else:
				compound_label = Compound_Name_Formal(label, type="latex")
				axes.text(coords[0], coords[1], str(count), fontdict=fontdict)
				newlabels.append("{} : {}".format(count, compound_label))
				count += 1
	elif type == "quaternary":
		for x, y, z in lines:
			axes.plot(x, y, z, "bo-", linewidth=linewidth, markeredgecolor="b", markerfacecolor="r", markersize=markersize)
		for coords in sorted(labels.keys()):
			entry = labels[coords]
			label = entry.name
			if len(entry.composition.elements) == 1:
				axes.text(coords[0], coords[1], coords[2], label, fontdict=fontdict)
			else:
				compound_label = Compound_Name_Formal(label, type="latex")
				axes.text(coords[0], coords[1], coords[2], str(count), fontdict=fontdict)
				newlabels.append("{} : {}".format(count, compound_label))
				count += 1
	
	# Remove Cartesian axes
	axes.axis("off")
	
	return newlabels


