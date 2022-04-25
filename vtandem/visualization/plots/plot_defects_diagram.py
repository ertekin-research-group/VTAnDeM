
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'
__author__ = 'Michael_Lidia_Jiaxing_Elif'

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from labellines import labelLine, labelLines

import PyQt5
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

from vtandem.visualization.utils.defect_formation_energy import *

from vtandem.visualization.plots.save_plot import SaveFigure


class Plot_DefectsDiagram(SaveFigure):
	
	def __init__(self, elements_list):
		
		# Font description for defect formation energy diagram
		self.font = {'family': 'sans-serif', 'color':  'black', 'weight': 'normal', 'size': 12 }
		
		# Store all extracted DFT data
		self.defects_data = {}
		self.main_compound_info = {}
		self.EVBM = 0.0
		self.ECBM = 0.0
		self.fermi_energy_array = None

		# Initialize all mu values
		self.mu_elements = {}
		for element in elements_list:
			self.mu_elements[element] = {"mu0": 0.0, "deltamu": 0.0}

		# Minimum and maximum y-value range
		self.axis_lims = {	"XMin": 0.0,
							"XMax": 1.0,
							"YMin": -2.0,
							"YMax": 2.0
							}
		
		# Store defect formation energy data
		self.intrinsic_defects_enthalpy_data = {}
		self.extrinsic_defects_enthalpy_data = {}
		self.dopant_enthalpy_data = None
		
		# Store defect formation plots and their labels
		self.intrinsic_defect_plots = {}
		self.extrinsic_defect_plots = {}
		self.dopant_plot = None
		self.defect_labels = {}
		
		# Store user-selected dopant
		self.dopant = "None"
		self.dopant_mu0 = 0.0
		self.dopant_deltamu = 0.0
		self.extrinsic_defects = []  # List of extrinsic defects of dopant (e.g. Ge_Bi, Ge_Se, and Ge_O for Ge dopant in Bi2O2Se)
		
		# Defects diagram
		self.defects_diagram_plot_figure = plt.figure()
		self.defects_diagram_plot_figure.subplots_adjust(left=0.225)
		self.defects_diagram_plot_drawing = self.defects_diagram_plot_figure.add_subplot(111)
		self.defects_diagram_plot_canvas = FigureCanvas(self.defects_diagram_plot_figure)
		
		# Save figure feature
		SaveFigure.__init__(self, self.defects_diagram_plot_figure)
		
		# Equilibrium Fermi energy vertical line
		self.equilibrium_fermi_energy_plot = None
		self.equilibrium_fermi_energy_tick = self.defects_diagram_plot_drawing.twiny()
	
	
	
	def Activate_DefectsDiagram_Plot_Axes(self):

		# Set plot axes limits (self.xmin and self.xmax are set in tab_phasediagram_...)
		self.defects_diagram_plot_drawing.set_xlim(self.axis_lims["XMin"], self.axis_lims["XMax"])
		self.defects_diagram_plot_drawing.set_ylim(self.axis_lims["YMin"], self.axis_lims["YMax"])
		
		# Set plot axes labels
		self.defects_diagram_plot_drawing.set_xlabel("Fermi Energy (eV)", fontdict=self.font)
		self.defects_diagram_plot_drawing.set_ylabel("$\Delta E_{D,q}$ (eV)", fontdict=self.font, rotation=90)

		# Set labels for VBM and CBM
		self.defects_diagram_plot_drawing.set_xticks([0.0, self.ECBM-self.EVBM])
		self.defects_diagram_plot_drawing.set_xticklabels(["VBM = 0.0", "CBM = "+str(round(self.ECBM-self.EVBM, 2))])

		# Set placement/direction of ticks and labels
		self.defects_diagram_plot_drawing.xaxis.tick_bottom()
		self.defects_diagram_plot_drawing.yaxis.tick_left()
		self.defects_diagram_plot_drawing.tick_params(axis='both', labelsize=self.font['size']-2)
		self.defects_diagram_plot_drawing.xaxis.set_label_position("bottom")
		self.defects_diagram_plot_drawing.yaxis.set_label_position("left")
		self.defects_diagram_plot_drawing.set_aspect("auto")

		# Color everything outside of band gap and below H=0
		self.defects_diagram_plot_drawing.fill_between(self.fermi_energy_array - self.EVBM, 0, -100, facecolor='#614126', interpolate=True, alpha=.1)
		self.defects_diagram_plot_drawing.fill_between(np.linspace(-1, 0, 100), 100, -100, facecolor='#614126', interpolate=True, alpha=.1)
		self.defects_diagram_plot_drawing.fill_between(np.linspace(self.ECBM-self.EVBM, self.ECBM-self.EVBM+1, 100), 100, -100, facecolor='#614126', interpolate=True, alpha=.1)
		
		# Settings for equilibrium Fermi energy
		#self.equilibrium_fermi_energy_tick.set_xlim(self.xmin, self.xmax)
		self.equilibrium_fermi_energy_tick.set_xlim(self.axis_lims["XMin"], self.axis_lims["XMax"])
		self.equilibrium_fermi_energy_tick.set_xticks([-100])
		self.equilibrium_fermi_energy_tick.tick_params(axis='both', labelsize=self.font['size']-2)
	
	
	
	
	def Update_WindowSize(self, axis_type, axislim_boxes):
		
		# Check if input is legitimate
		try:
			float(axislim_boxes[axis_type].text())
		except:
			axislim_boxes[axis_type].setText(str(self.axis_lims[axis_type]))
			return

		# Check if axes bounds are legimitate
		axis = axis_type[0]
		if float(axislim_boxes[axis+"Min"].text()) > float(axislim_boxes[axis+"Max"].text()):
			axislim_boxes[axis_type].setText(str(self.axis_lims[axis_type]))
			return

		self.axis_lims[axis_type] = float(axislim_boxes[axis_type].text())
		self.defects_diagram_plot_drawing.set_xlim(self.axis_lims["XMin"], self.axis_lims["XMax"])
		self.defects_diagram_plot_drawing.set_ylim(self.axis_lims["YMin"], self.axis_lims["YMax"])
		self.equilibrium_fermi_energy_tick.set_xlim(self.axis_lims["XMin"], self.axis_lims["XMax"])
		self.defects_diagram_plot_canvas.draw()
	

	
	def Update_Deltamus(self, deltamu_values):

		# Args:	
		# 	deltamu_values: Dictionary of deltamu values, in element:value pairs
		for element in self.mu_elements.keys():
			self.mu_elements[element]["deltamu"] = deltamu_values[element]


	def Calculate_DefectFormations(self):

		intrinsic_defects_enthalpy_data = Calculate_IntrinsicDefectFormationEnthalpies(	self.defects_data, \
																						self.main_compound_info, \
																						self.fermi_energy_array, \
																						self.mu_elements	)
		self.intrinsic_defects_enthalpy_data = Find_MinimumDefectFormationEnthalpies(intrinsic_defects_enthalpy_data)
		
		if self.dopant != "None":
			
			extrinsic_defects_enthalpy_data = Calculate_ExtrinsicDefectFormationEnthalpies(	self.defects_data, \
																							self.main_compound_info, \
																							self.fermi_energy_array, \
																							self.mu_elements, \
																							self.extrinsic_defects, \
																							self.dopant, \
																							self.dopant_mu0, \
																							self.dopant_deltamu	)
			self.extrinsic_defects_enthalpy_data = Find_MinimumDefectFormationEnthalpies(extrinsic_defects_enthalpy_data)
	
	
	
	def Initialize_Intrinsic_DefectsDiagram_Plot(self):
		
		# Plot defect formation energy of each intrinsic defect
		for intrinsic_defect in self.intrinsic_defects_enthalpy_data.keys():
			defect_label = r""+intrinsic_defect.split("_")[0]+"$_\mathrm{"+intrinsic_defect.split("_")[-1]+"}$"
			self.intrinsic_defect_plots[intrinsic_defect], = self.defects_diagram_plot_drawing.plot(self.fermi_energy_array - self.EVBM, self.intrinsic_defects_enthalpy_data[intrinsic_defect], label = defect_label)
		
		# Create label for each defect
		try:
			labelLines(list(self.intrinsic_defect_plots.values()), align = False, xvals = np.linspace(0.0, self.ECBM-self.EVBM, len(self.intrinsic_defect_plots.keys())+2)[1:len(self.intrinsic_defect_plots.keys())+1], fontsize = 10, bbox = dict(facecolor = 'white', alpha = 0.8, edgecolor = 'white', pad = 0.5))
		except:
			pass
		
		# Draw defects diagram canvas
		self.defects_diagram_plot_canvas.draw()
	
	
	
	def Update_Intrinsic_DefectsDiagram_Plot(self):
		
		# Update defect formation energy of each intrinsic defect
		for intrinsic_defect in self.intrinsic_defects_enthalpy_data.keys():
			self.intrinsic_defect_plots[intrinsic_defect].set_ydata(self.intrinsic_defects_enthalpy_data[intrinsic_defect])
		
		# Remove labels before redrawing them at new positions
		self.defects_diagram_plot_drawing.texts.clear()
		labelLines(list(self.intrinsic_defect_plots.values()), align = False, xvals = np.linspace(0.0, self.ECBM-self.EVBM, len(self.intrinsic_defect_plots.keys())+2)[1:len(self.intrinsic_defect_plots.keys())+1], fontsize = 10, bbox = dict(facecolor = 'white', alpha = 0.8, edgecolor = 'white', pad = 0.5))
		
		# Draw defects diagram canvas
		self.defects_diagram_plot_canvas.draw()
	
	
	
	def Initialize_Extrinsic_DefectsDiagram_Plot(self):

		for extrinsic_defect in self.extrinsic_defects:
			
			# Check that extrinsic defect involves the dopant atom (e.g. Ge_Bi, Ge_Se, Ge_O if dopant = Ge)
			if extrinsic_defect.split("_")[0] != self.dopant:
				continue

			# Plot defect formation energy of dopant
			defect_label = r""+extrinsic_defect.split("_")[0]+"$_\mathrm{"+extrinsic_defect.split("_")[-1]+"}$"
			self.extrinsic_defect_plots[extrinsic_defect], = self.defects_diagram_plot_drawing.plot(self.fermi_energy_array - self.EVBM, self.extrinsic_defects_enthalpy_data[extrinsic_defect], label = defect_label)
			
			# Create label for each defect
			labelLine(self.extrinsic_defect_plots[extrinsic_defect], x = (self.ECBM-self.EVBM)/2., align = False, fontsize = 10, bbox = dict(facecolor = 'white', alpha = 0.8, edgecolor = 'white', pad = 0.5))
		
		# Draw defects diagram canvas
		self.defects_diagram_plot_canvas.draw()
	

	
	def Update_Extrinsic_DefectsDiagram_Plot(self):
		
		for extrinsic_defect in self.extrinsic_defects:
			
			# Check that extrinsic defect involves the dopant atom (e.g. Ge_Bi, Ge_Se, Ge_O if dopant = Ge)
			if extrinsic_defect.split("_")[0] != self.dopant:
				continue

			# Update defect formation energy of dopant
			self.extrinsic_defect_plots[extrinsic_defect].set_ydata(self.extrinsic_defects_enthalpy_data[extrinsic_defect])
			
			# Remove labels before redrawing them at new positions
			labelLine(self.extrinsic_defect_plots[extrinsic_defect], x = (self.ECBM-self.EVBM)/2., align = False)
		
		# Draw defects diagram canvas
		self.defects_diagram_plot_canvas.draw()
	
	
	
	def Plot_Equilibrium_Fermi_Energy(self, equilibrium_fermi_energy):
		
		# Plot equilibrium Fermi energy
		try:
			self.equilibrium_fermi_energy_plot.remove()
			self.equilibrium_fermi_energy_plot = self.defects_diagram_plot_drawing.axvline(equilibrium_fermi_energy, zorder=1E9, ls='--', color='k')
		except:
			self.equilibrium_fermi_energy_plot = self.defects_diagram_plot_drawing.axvline(equilibrium_fermi_energy, zorder=1E9, ls='--', color='k')
		
		# Place EF^eq text
		try:
			self.equilibrium_fermi_energy_tick.set_xticks([equilibrium_fermi_energy])
			self.equilibrium_fermi_energy_tick.set_xticklabels([r"$E_{f}^{eq}$"])
		except:
			pass

		# Draw defects diagram canvas
		self.defects_diagram_plot_canvas.draw()














