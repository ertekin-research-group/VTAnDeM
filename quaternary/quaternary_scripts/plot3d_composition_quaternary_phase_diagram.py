

import itertools
import periodictable
import matplotlib.pyplot as plt
from scipy.spatial import ConvexHull
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from pymatgen import MPRester
from pymatgen.analysis.phase_diagram import PhaseDiagram, PDPlotter, PDEntry
import numpy as np

import PyQt5
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *


class Composition_Quaternary_PhaseDiagram3D(object):

	def __init__(self, parent = None, main_compound = None, first_element = None, second_element = None, third_element = None, fourth_element = None):
		
		# All elements in the periodic table
		self.all_elements = []
		for element in periodictable.elements:
			self.all_elements.append(str(element))
		
		# Font description for phase stability diagram plot
		self.font = {'family': 'sans-serif', 'color':  'black', 'weight': 'bold', 'size': 14 }
		
		# Establish the first, second, third, and fourth species of the quaternary compound.
		self.main_compound  = main_compound
		self.first_element	= first_element
		self.second_element	= second_element
		self.third_element	= third_element
		self.fourth_element	= fourth_element
		self.elements_list  = [self.first_element, self.second_element, self.third_element, self.fourth_element]
		
		# Store all extracted DFT data
		self.compounds_info = {}
		
		"""
		# Initialize necessary phase diagram objects
		self.phasediagram_entries = []
		"""
		self.phasediagram = None
		self.phasediagram_plot_object = None
		
		
		# Quaternary phase diagram (in composition space) object
		self.composition_phasediagram_3d_plot_figure = plt.figure()
		self.composition_phasediagram_3d_plot_canvas = FigureCanvas(self.composition_phasediagram_3d_plot_figure)
		self.composition_phasediagram_3d_plot_drawing = Axes3D(self.composition_phasediagram_3d_plot_figure)
	
	
	def Create_Compositional_PhaseDiagram3D(self):
		
		phasediagram_entries = []
		
		for compound in self.compounds_info.keys():
			
			compound_composition = {}
			
			# Disregard elements not included in main compound
			if (compound in self.all_elements) and (compound not in self.elements_list):
				continue
			
			# Elements
			if compound in self.elements_list:
				compound_composition[compound] = self.compounds_info[compound]["dft_"+compound]
			
			# Compounds
			else:
				for element in self.compounds_info[compound]["elements_list"]:
					compound_composition[element] = self.compounds_info[compound]["dft_"+element]
			
			compound_total_energy = self.compounds_info[compound]["total_energy"]
			
			phasediagram_entries.append(PDEntry(compound_composition, compound_total_energy))
		
		
		self.phasediagram = PhaseDiagram(phasediagram_entries)
	
	
	
	
	def Plot_Compositional_PhaseDiagram3D(self, label_stable=True):
		
		self.phasediagram_plot_object = PDPlotter(self.phasediagram)
		(lines, labels, unstable) = self.phasediagram_plot_object.pd_plot_data
		
		count = 1
		newlabels = list()
		for x, y, z in lines:
			self.composition_phasediagram_3d_plot_drawing.plot(x, y, z, "bo-", linewidth=3, markeredgecolor="b", markerfacecolor="r", markersize=10)
		for coords in sorted(labels.keys()):
			entry = labels[coords]
			label = entry.name
			if label_stable:
				if len(entry.composition.elements) == 1:
					self.composition_phasediagram_3d_plot_drawing.text(coords[0], coords[1], coords[2], label, fontdict=self.font)
				else:
					self.composition_phasediagram_3d_plot_drawing.text(coords[0], coords[1], coords[2], str(count), fontdict=self.font)
					newlabels.append("{} : {}".format(count, label))
					count += 1
		plt.figtext(0.01, 0.01, "\n".join(newlabels))
		self.composition_phasediagram_3d_plot_drawing.axis("off")
		self.composition_phasediagram_3d_plot_canvas.draw()










