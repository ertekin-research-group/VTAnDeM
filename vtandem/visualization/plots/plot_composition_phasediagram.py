
__author__ = 'Michael_Lidia_Jiaxing_Benita_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'


import numpy as np
import itertools
import periodictable
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter, PDEntry
from pymatgen.core.composition import Composition


class Plot_Composition_PhaseDiagram:
	
	def __init__(self, type: str):
		
		# Record type of phase diagram (either ternary or quaternary)
		if type not in ["ternary", "quaternary"]:
			raise ValueError("Argument 'type' can be either 'ternary' or 'quaternary'.")
		self.type = type
		
		# All elements in the periodic table
		self.all_elements = []
		for element in periodictable.elements:
			self.all_elements.append(str(element))
		
		# Font description for phase stability diagram plot
		#self.font = {'family': 'sans-serif', 'color':  'black', 'weight': 'bold', 'size': 14 }
		self.font = {'family': 'sans-serif', 'color':  'black', 'weight': 'bold', 'size': 10 }
		
		# Store all extracted DFT data
		self.compounds_info = {}
		
		# Initialize pymatgen phase diagram objects
		self.pmg_phasediagram = None
		self.pmg_phasediagram_plot_object = None
		
		# Phase diagram (in composition space) object
		self.composition_phasediagram_plot_figure = plt.figure()
		self.composition_phasediagram_plot_canvas = FigureCanvas(self.composition_phasediagram_plot_figure)
		if self.type == "ternary":
			self.composition_phasediagram_plot_drawing = self.composition_phasediagram_plot_figure.add_subplot(111)
		if self.type == "quaternary":
			self.composition_phasediagram_plot_drawing = self.composition_phasediagram_plot_figure.add_subplot(111, projection='3d')
		
		# Store all lines and labels of the compositional phase diagram
		self.lines = []
		self.labels = {}
		
		# Store all plots of phase diagram vertices
		self.vertices_plots = {}
	
	
	
	def Create_Compositional_PhaseDiagram(self):
		
		# Record all entries for the phase diagram
		phasediagram_entries = []
		
		for compound in self.compounds_info.keys():
			
			# Disregard elements not included in main compound
			if (compound in self.all_elements) and (compound not in self.elements_list):
				continue
			
			# Get the compound's composition
			compound_composition = {}
			if compound in self.elements_list:	# Elemental material
				compound_composition[compound] = self.compounds_info[compound]["dft_"+compound]
			else:	# Compound material
				for element in self.compounds_info[compound]["elements_list"]:
					compound_composition[element] = self.compounds_info[compound]["dft_"+element]
			
			# Get the compound's total energy
			compound_total_energy = self.compounds_info[compound]["total_energy"]
			
			# Record to list of entries
			phasediagram_entries.append(PDEntry(compound_composition, compound_total_energy, compound))
		
		# Calculate compositional phase diagram (using pymatgen)
		#	The output data structure is as follows:
		#		lines --> List of arrays, each array is 2x2 for ternary (3x3 for quaternary, etc.), column vector represents point on phase diagram.
		#					ex: array([ [0.3, 0.5], [1.0, 0.0] ]) is a line that goes from point [x=0.3, y=1.0] to point [x=0.5, y=0.0]
		#		labels --> Dictionary with point-PDEntry pairs.
		self.pmg_phasediagram = PhaseDiagram(phasediagram_entries)
		self.pmg_phasediagram_plot_object = PDPlotter(self.pmg_phasediagram)
		(lines, labels, unstable) = self.pmg_phasediagram_plot_object.pd_plot_data
		
		
		y = self.pmg_phasediagram.as_dict()
		x = self.pmg_phasediagram.get_all_chempots(Composition(self.main_compound))
		for key in x.keys():
			print(key, x[key], list(x[key].values()) )
			#print(type(key), key)
		
		elements = []
		key = list(x.keys())[0]
		for element in x[key].keys():
			elements.append(element)
		cu = elements[0]
		bi = elements[1]
		se = elements[2]
		o = elements[3]
		print(bi, cu, se, o)
		phaseregions = []
		chempots = []
		for key in x.keys():
			phaseregions.append(key)
			chempots.append( np.array([ x[key][bi], x[key][cu], x[key][se], x[key][o] ]) )
		chempots = np.asarray(chempots).T
		for x in chempots:
			print(x)
		print(phaseregions)
		
		# Record all lines and points of the compositional phase diagram
		self.lines = lines
		self.labels = labels
	
	
	
	def Plot_Compositional_PhaseDiagram(self, label_stable=True):
		
		# Plot settings
		linewidth = 1
		markersize = 5
		
		# Plot compositional phase diagram
		count = 1
		newlabels = []
		if self.type == "ternary":
			for x, y in self.lines:
				#self.composition_phasediagram_plot_drawing.plot(x, y, "bo-", linewidth=3, markeredgecolor="b", markerfacecolor="r", markersize=10)
				self.composition_phasediagram_plot_drawing.plot(x, y, "bo-", linewidth=linewidth, markeredgecolor="b", markerfacecolor="r", markersize=markersize)
			for coords in sorted(self.labels.keys()):
				entry = self.labels[coords]
				label = entry.name
				if label_stable:
					if len(entry.composition.elements) == 1:
						self.composition_phasediagram_plot_drawing.text(coords[0], coords[1], label, fontdict=self.font)
					else:
						self.composition_phasediagram_plot_drawing.text(coords[0], coords[1], str(count), fontdict=self.font)
						newlabels.append("{} : {}".format(count, label))
						count += 1
		elif self.type == "quaternary":
			for x, y, z in self.lines:
				#self.composition_phasediagram_plot_drawing.plot(x, y, z, "bo-", linewidth=3, markeredgecolor="b", markerfacecolor="r", markersize=10)
				self.composition_phasediagram_plot_drawing.plot(x, y, z, "bo-", linewidth=linewidth, markeredgecolor="b", markerfacecolor="r", markersize=markersize)
			for coords in sorted(self.labels.keys()):
				entry = self.labels[coords]
				label = entry.name
				if label_stable:
					if len(entry.composition.elements) == 1:
						self.composition_phasediagram_plot_drawing.text(coords[0], coords[1], coords[2], label, fontdict=self.font)
					else:
						self.composition_phasediagram_plot_drawing.text(coords[0], coords[1], coords[2], str(count), fontdict=self.font)
						newlabels.append("{} : {}".format(count, label))
						count += 1
		
		# Draw compositional phase diagram
		self.composition_phasediagram_legend = self.composition_phasediagram_plot_figure.text(0.01, 0.01, "\n".join(newlabels))
		self.composition_phasediagram_plot_drawing.axis("off")
		self.composition_phasediagram_plot_canvas.draw()
	
	
	
	
	###############################################################################################
	###################################### Save Figure ############################################
	###############################################################################################
	
	def SaveFigure(self, filename, extension_type):
		
		if filename:
			extension = extension_type.split(".")[-1].split(")")[0]
			if filename.split(".")[-1] == extension:
				self.composition_phasediagram_plot_figure.savefig(filename, bbox_inches='tight')
			else:
				self.composition_phasediagram_plot_figure.savefig(filename+"."+extension, bbox_inches='tight')










