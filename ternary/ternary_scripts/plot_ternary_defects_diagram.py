
__author__ = 'Michael_Lidia_Jiaxing_Benita_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'


import numpy as np
import itertools
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from labellines import labelLine, labelLines


class Ternary_Defects_Diagram(object):
	
	def __init__(self, parent = None, main_compound = None, first_element = None, second_element = None, third_element = None):
		
		# Font description for phase stability diagram plot
		self.font = {'family': 'sans-serif',
				'color':  'black',
				'weight': 'normal',
				'size': 14 }
		
		# Establish the first, second, third, and fourth species of the quaternary compound.
		# Note that this list is subject to change, depending on what the user chooses.
		self.main_compound  = main_compound
		self.first_element  = first_element
		self.second_element = second_element
		self.third_element  = third_element
		self.elements_list = [self.first_element, self.second_element, self.third_element]
		
		
		# Keep track of chemical potential values
		self.mu_elements = {self.first_element: {"mu0": 0.0, "deltamu": 0.0},
							self.second_element: {"mu0": 0.0, "deltamu": 0.0},
							self.third_element: {"mu0": 0.0, "deltamu": 0.0} }
		
		
		"""
		# Establish constants for the species
		self.species_a = first_element
		self.species_b = second_element
		self.species_c = third_element
		"""
		
		# Store all extracted DFT data
		self.ternary_defects_data = None
		self.main_compound_total_energy = 0.0
		self.first_element_mu0 = 0.0
		self.second_element_mu0 = 0.0
		self.third_element_mu0 = 0.0
		self.EVBM = 0.0
		self.ECBM = 0.0
		self.fermi_energy_array = None
		
		"""
		# Store all mu values
		self.mu_values = {}
		self.mu_values[self.first_element]  = 0.0
		self.mu_values[self.second_element] = 0.0
		self.mu_values[self.third_element]  = 0.0
		"""
		
		# Minimum and maximum y-value range
		self.ymin = -2.0
		self.ymax = 2.0
		
		# Store defect formation energy data
		self.intrinsic_defects_enthalpy_data = {}
		self.extrinsic_defects_enthalpy_data = {}
		
		# Store defect formation plots and their labels
		self.intrinsic_defect_plots = {}
		self.extrinsic_defect_plots = {}
		self.defect_labels = {}
		
		# Store user-selected dopant
		self.extrinsic_defect = "None"
		self.extrinsic_defect_mu0 = 0.0
		self.extrinsic_defect_deltamu = 0.0
		
		# Defects diagram plot
		self.ternary_defects_diagram_plot_figure = plt.figure()
		self.ternary_defects_diagram_plot_figure.subplots_adjust(left=0.225)
		self.ternary_defects_diagram_plot_drawing = self.ternary_defects_diagram_plot_figure.add_subplot(111)
		self.ternary_defects_diagram_plot_canvas = FigureCanvas(self.ternary_defects_diagram_plot_figure)
		
		# Equilibrium Fermi energy vertical line
		self.equilibrium_fermi_energy_plot = None
		self.equilibrium_fermi_energy_tick = self.ternary_defects_diagram_plot_drawing.twiny()
	
	
	
	def Activate_DefectsDiagram_Plot_Axes(self):
		
		self.ternary_defects_diagram_plot_drawing.set_xlim(0.0, self.ECBM-self.EVBM)
		self.ternary_defects_diagram_plot_drawing.set_ylim(self.ymin, self.ymax)
		self.ternary_defects_diagram_plot_drawing.set_xlabel("Fermi Energy (eV)", fontdict=self.font)
		self.ternary_defects_diagram_plot_drawing.set_ylabel("$\Delta$H(q,D) (eV)", fontdict=self.font, rotation=90)
		self.ternary_defects_diagram_plot_drawing.set_xticks([0.0, self.ECBM-self.EVBM])
		self.ternary_defects_diagram_plot_drawing.set_xticklabels(["VBM=0.0", "CBM="+str(round(self.ECBM-self.EVBM, 2))])
		self.ternary_defects_diagram_plot_drawing.xaxis.tick_bottom()
		self.ternary_defects_diagram_plot_drawing.yaxis.tick_left()
		self.ternary_defects_diagram_plot_drawing.tick_params(axis='both', labelsize=9)
		self.ternary_defects_diagram_plot_drawing.xaxis.set_label_position("bottom")
		self.ternary_defects_diagram_plot_drawing.yaxis.set_label_position("left")
		self.ternary_defects_diagram_plot_drawing.set_aspect("auto")
		self.ternary_defects_diagram_plot_drawing.fill_between(self.fermi_energy_array - self.EVBM, 0, -100, facecolor='#614126', interpolate=True, alpha=.1)
		self.equilibrium_fermi_energy_tick.set_xlim(0.0, self.ECBM-self.EVBM)
		self.equilibrium_fermi_energy_tick.set_xticks([])
		self.equilibrium_fermi_energy_tick.tick_params(axis='both', labelsize=9)
	
	
	def Calculate_DefectFormations(self):
		
		"""
		# Obtain current user-tunable chemical potentials
		dmu_a = self.mu_values[self.species_a]
		dmu_b = self.mu_values[self.species_b]
		dmu_c = self.mu_values[self.species_c]
		"""
		
		# Initialize storage for all charges of defect
		intrinsic_defects_enthalpy_data = {}
		extrinsic_defects_enthalpy_data = {}
		
		# Loop through defects in ternary
		for defect in self.ternary_defects_data.keys():
			
			# Check that item is truly a defect
			if "_" not in defect:
				continue
			
			# Extrinsic defect
			if self.ternary_defects_data[defect]["Extrinsic"] == "Yes":
				extrinsic_defects_enthalpy_data[defect] = []
				for charge in self.ternary_defects_data[defect]["charge"].keys():
					"""
					defect_formation_enthalpy = self.ternary_defects_data[defect]["charge"][charge]["Energy"] \
												- self.ternary_defects_data["supercellsize"] * self.main_compound_total_energy \
												- self.ternary_defects_data[defect]["n_"+self.species_a] * ( self.first_element_mu0 + float(dmu_a) ) \
												- self.ternary_defects_data[defect]["n_"+self.species_b] * ( self.second_element_mu0 + float(dmu_b) ) \
												- self.ternary_defects_data[defect]["n_"+self.species_c] * ( self.third_element_mu0 + float(dmu_c) ) \
												- (self.extrinsic_defect_mu0 + self.extrinsic_defect_deltamu) \
												+ float(charge) * self.fermi_energy_array \
												+ self.ternary_defects_data[defect]["charge"][charge]["ECorr"] # Extra chemical potential contribution from dopant
					"""
					defect_formation_enthalpy = self.ternary_defects_data[defect]["charge"][charge]["Energy"] \
												- self.ternary_defects_data["supercellsize"] * self.main_compound_total_energy \
												- (self.extrinsic_defect_mu0 + self.extrinsic_defect_deltamu) \
												+ float(charge) * self.fermi_energy_array \
												+ self.ternary_defects_data[defect]["charge"][charge]["ECorr"] # Extra chemical potential contribution from dopant
					for element in self.elements_list:
						defect_formation_enthalpy -= self.ternary_defects_data[defect]["n_"+element] * ( self.mu_elements[element]["mu0"] + self.mu_elements[element]["deltamu"] )
					extrinsic_defects_enthalpy_data[defect].append(defect_formation_enthalpy)
			# Intrinsic defect
			elif self.ternary_defects_data[defect]["Extrinsic"] == "No":
				intrinsic_defects_enthalpy_data[defect] = []
				for charge in self.ternary_defects_data[defect]["charge"].keys():
					"""
					defect_formation_enthalpy = self.ternary_defects_data[defect]["charge"][charge]["Energy"] \
												- self.ternary_defects_data["supercellsize"] * self.main_compound_total_energy \
												- self.ternary_defects_data[defect]["n_"+self.species_a] * ( self.first_element_mu0 + float(dmu_a) ) \
												- self.ternary_defects_data[defect]["n_"+self.species_b] * ( self.second_element_mu0 + float(dmu_b) ) \
												- self.ternary_defects_data[defect]["n_"+self.species_c] * ( self.third_element_mu0 + float(dmu_c) ) \
												+ float(charge) * self.fermi_energy_array \
												+ self.ternary_defects_data[defect]["charge"][charge]["ECorr"]
					"""
					defect_formation_enthalpy = self.ternary_defects_data[defect]["charge"][charge]["Energy"] \
												- self.ternary_defects_data["supercellsize"] * self.main_compound_total_energy \
												+ float(charge) * self.fermi_energy_array \
												+ self.ternary_defects_data[defect]["charge"][charge]["ECorr"]
					for element in self.elements_list:
						defect_formation_enthalpy -= self.ternary_defects_data[defect]["n_"+element] * ( self.mu_elements[element]["mu0"] + self.mu_elements[element]["deltamu"] )
					intrinsic_defects_enthalpy_data[defect].append(defect_formation_enthalpy)
		
		
		# Find minimum formation energy of all charges of each defect
		for intrinsic_defect in intrinsic_defects_enthalpy_data.keys():	# Intrinsic defects
			defect_formation_energy_minimum = np.fromiter(map(min, zip(*itertools.chain(intrinsic_defects_enthalpy_data[intrinsic_defect]))), dtype=np.float)
			self.intrinsic_defects_enthalpy_data[intrinsic_defect] = defect_formation_energy_minimum
		for extrinsic_defect in extrinsic_defects_enthalpy_data.keys():	# Extrinsic defects
			defect_formation_energy_minimum = np.fromiter(map(min, zip(*itertools.chain(extrinsic_defects_enthalpy_data[extrinsic_defect]))), dtype=np.float)
			self.extrinsic_defects_enthalpy_data[extrinsic_defect] = defect_formation_energy_minimum
	
	
	
	
	def Initialize_Intrinsic_DefectsDiagram_Plot(self):
		
		# Plot defect formation energy of each intrinsic defect
		for intrinsic_defect in self.intrinsic_defects_enthalpy_data.keys():
			defect_label = r"$"+intrinsic_defect.split("_")[0]+"_{"+intrinsic_defect.split("_")[-1]+"}$"
			self.intrinsic_defect_plots[intrinsic_defect], = self.ternary_defects_diagram_plot_drawing.plot(self.fermi_energy_array - self.EVBM, self.intrinsic_defects_enthalpy_data[intrinsic_defect], label = defect_label)
		
		# Create label for each defect
		labelLines(list(self.intrinsic_defect_plots.values()), align = False, xvals = np.linspace(0.0, self.ECBM-self.EVBM, len(self.intrinsic_defect_plots.keys())+2)[1:len(self.intrinsic_defect_plots.keys())+1], fontsize = 10, bbox = dict(facecolor = 'white', alpha = 0.8, edgecolor = 'white', pad = 0.5))
		
		# Draw defects diagram canvas
		self.ternary_defects_diagram_plot_canvas.draw()
	
	
	
	def Update_Intrinsic_DefectsDiagram_Plot(self):
		
		# Update defect formation energy of each intrinsic defect
		for intrinsic_defect in self.intrinsic_defects_enthalpy_data.keys():
			self.intrinsic_defect_plots[intrinsic_defect].set_ydata(self.intrinsic_defects_enthalpy_data[intrinsic_defect])
		
		# Remove labels before redrawing them at new positions
		self.ternary_defects_diagram_plot_drawing.texts.clear()
		labelLines(list(self.intrinsic_defect_plots.values()), align = False, xvals = np.linspace(0.0, self.ECBM-self.EVBM, len(self.intrinsic_defect_plots.keys())+2)[1:len(self.intrinsic_defect_plots.keys())+1], fontsize = 10, bbox = dict(facecolor = 'white', alpha = 0.8, edgecolor = 'white', pad = 0.5))
		
		# Draw defects diagram canvas
		self.ternary_defects_diagram_plot_canvas.draw()
	
	
	
	def Initialize_Extrinsic_DefectsDiagram_Plot(self):
		
		# Plot defect formation energy of dopant
		defect_label = r"$"+self.extrinsic_defect.split("_")[0]+"_{"+self.extrinsic_defect.split("_")[-1]+"}$"
		self.extrinsic_defect_plots[self.extrinsic_defect], = self.ternary_defects_diagram_plot_drawing.plot(self.fermi_energy_array - self.EVBM, self.extrinsic_defects_enthalpy_data[self.extrinsic_defect], label = defect_label)
		
		# Create label for each defect	
		labelLine(self.extrinsic_defect_plots[self.extrinsic_defect], x=(self.ECBM-self.EVBM)/2., align=False)
		
		# Draw defects diagram canvas
		self.ternary_defects_diagram_plot_canvas.draw()
	
	
	
	def Update_Extrinsic_DefectsDiagram_Plot(self):
		
		# Update defect formation energy of dopant
		self.extrinsic_defect_plots[self.extrinsic_defect].set_ydata(self.extrinsic_defects_enthalpy_data[self.extrinsic_defect])
		
		# Remove labels before redrawing them at new positions
		#self.ternary_defects_diagram_plot_drawing.texts.clear()
		labelLine(self.extrinsic_defect_plots[self.extrinsic_defect], x=(self.ECBM-self.EVBM)/2., align=False)
		
		# Draw defects diagram canvas
		self.ternary_defects_diagram_plot_canvas.draw()
	
	
	
	
	def Plot_Equilibrium_Fermi_Energy(self, temperature, equilibrium_fermi_energy):
		
		# Plot equilibrium Fermi energy
		try:
			self.equilibrium_fermi_energy_plot.remove()
			self.equilibrium_fermi_energy_plot = self.ternary_defects_diagram_plot_drawing.axvline(equilibrium_fermi_energy, zorder=1E9, ls='--', color='k')
		except:
			self.equilibrium_fermi_energy_plot = self.ternary_defects_diagram_plot_drawing.axvline(equilibrium_fermi_energy, zorder=1E9, ls='--', color='k')
		
		# Place EF^eq text
		"""
		if (equilibrium_fermi_energy > 0.0) and (equilibrium_fermi_energy < self.ECBM-self.EVBM):
			self.equilibrium_fermi_energy_tick.set_xticks([equilibrium_fermi_energy])
			self.equilibrium_fermi_energy_tick.set_xticklabels([r"$E_{f}^{eq}$"])
		else:
			self.equilibrium_fermi_energy_tick.set_xticks([])
			self.equilibrium_fermi_energy_tick.set_xticklabels([])
		"""
		try:
			self.equilibrium_fermi_energy_tick.set_xticks([equilibrium_fermi_energy])
			self.equilibrium_fermi_energy_tick.set_xticklabels([r"$E_{f}^{eq}$"])
		except:
			self.equilibrium_fermi_energy_tick.set_xticks([])
			self.equilibrium_fermi_energy_tick.set_xticklabels([])
		
		# Draw defects diagram canvas
		self.ternary_defects_diagram_plot_canvas.draw()
















