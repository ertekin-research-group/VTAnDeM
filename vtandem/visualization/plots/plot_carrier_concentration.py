
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'
__author__ = 'Michael_Lidia_Jiaxing_Elif'

import numpy as np
import os
import matplotlib.pyplot as plt
from scipy import integrate
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure

# Import functions for carrier concentration
from vtandem.visualization.utils.carrier_concentration import Calculate_CarrierConcentration, Calculate_FreeHole_FreeElectron_Concentrations
from vtandem.visualization.plots.save_plot import SaveFigure

"""
from PyQt5.QtWidgets import QMessageBox, QWidget, QMainWindow
from PyQt5.QtGui import QIcon
"""
from PyQt5.QtWidgets import *
from PyQt5.QtGui import *

script_path = os.path.dirname(__file__)
vtandem_source_path = "/".join(script_path.split("/")[:-1])


class Plot_CarrierConcentration(SaveFigure):
	
	def __init__(self, elements_list):
		
		# Font description for defect formation energy diagram
		self.font = { 'color': 'black', 'weight': 'normal', 'size': 12 }
		
		# Plot settings
		self.max_temperature = 1000  # Maximum temperature in Kelvins
		self.temperature_stepsize = 50

		# Store all extracted DFT data
		self.defects_data = None
		self.main_compound_info = None
		self.dos_data = None
		self.vol = 0.0	# Volume of defect supercell (NOT the DOS cell)
		self.EVBM = 0.0
		self.ECBM = 0.0
		self.fermi_energy_array = None
		self.temperature_array = np.arange(200, self.max_temperature+1, self.temperature_stepsize)
		self.synthesis_temperature = None
		
		# DOS info
		self.energy = None
		self.gE = None
		self.energies_ValenceBand		= None
		self.gE_ValenceBand 			= None
		self.energies_ConductionBand	= None
		self.gE_ConductionBand 			= None

		# Type of DOS (DFT or Effective Mass)
		self.dos_type = "DFT"

		# Effective masses of DOS
		self.effmass_dict = {"holes": 1.0, "electrons": 1.0}
		#self.effmass_dict = {"holes": 0.99, "electrons": 0.99}

		# Distribution statistics for defect concentration
		self.defectconc_stat = "Maxwell-Boltzmann"

		# Equilibrium Fermi energy
		self.intrinsic_equilibrium_fermi_energy = {}
		self.total_equilibrium_fermi_energy = {}
		for temperature in self.temperature_array:
			self.intrinsic_equilibrium_fermi_energy[temperature] = 0.0
			self.total_equilibrium_fermi_energy[temperature] = 0.0
		
		# Initialize all mu values
		self.mu_elements = {}
		for element in elements_list:
			self.mu_elements[element] = {"mu0": 0.0, "deltamu": 0.0}
		
		# Store user-selected dopant
		self.dopant = "None"
		self.dopant_mu0 = 0.0
		self.dopant_deltamu = 0.0
		self.extrinsic_defects = []  # List of extrinsic defects of dopant (e.g. Ge_Bi, Ge_Se, and Ge_O for Ge dopant in Bi2O2Se)
		
		# (WIDGET) Carrier Concentration Plot
		self.carrier_concentration_plot_figure = plt.figure()
		self.carrier_concentration_plot_figure.subplots_adjust(left=0.225)
		self.carrier_concentration_plot_drawing = self.carrier_concentration_plot_figure.add_subplot(111)
		self.carrier_concentration_plot_canvas = FigureCanvas(self.carrier_concentration_plot_figure)
		
		self.carrier_concentration_intrinsic_defect_hole_plot = None
		self.carrier_concentration_intrinsic_defect_electron_plot = None
		self.carrier_concentration_total_hole_plot = None
		self.carrier_concentration_total_electron_plot = None
		
		# Set y-axis minimum and maximum
		self.ymin = 1E16
		self.ymax = 1E23
		
		# Save figure feature
		SaveFigure.__init__(self, self.carrier_concentration_plot_figure)
	
	

	def Activate_CarrierConcentration_Plot_Axes(self):
		
		self.carrier_concentration_plot_drawing.set_xlim(200, self.max_temperature)
		self.carrier_concentration_plot_drawing.set_ylim(self.ymin, self.ymax)
		self.carrier_concentration_plot_drawing.set_xlabel("T(K)", fontdict=self.font)
		self.carrier_concentration_plot_drawing.set_ylabel("Carrier Concentration (cm$^{-3}$)", fontdict=self.font, rotation=90)
		self.carrier_concentration_plot_drawing.set_yscale("log")
		self.carrier_concentration_plot_drawing.xaxis.tick_bottom()
		self.carrier_concentration_plot_drawing.yaxis.tick_left()
		self.carrier_concentration_plot_drawing.tick_params(axis='both', labelsize=self.font['size']-2)
		self.carrier_concentration_plot_drawing.xaxis.set_label_position("bottom")
		self.carrier_concentration_plot_drawing.yaxis.set_label_position("left")
		self.carrier_concentration_plot_drawing.set_aspect("auto")
	
	

	def Update_WindowSize(self, ytype, Ylim_box_object):
		
		# Check to make sure y-limit input is legitimate number
		try:
			float(Ylim_box_object.text())
		except:
			self.Revert_WindowSize_toDefault(ytype, Ylim_box_object)
			return

		# Check to make sure y-limit input is greater than 0
		if float(Ylim_box_object.text()) <= 0.0:
			self.Revert_WindowSize_toDefault(ytype, Ylim_box_object)
			return

		# Modify defects diagram y-axis
		if ytype == "YMin":
			self.ymin = float(Ylim_box_object.text())
		if ytype == "YMax":
			self.ymax = float(Ylim_box_object.text())
		self.carrier_concentration_plot_drawing.set_ylim(self.ymin, self.ymax)
		self.carrier_concentration_plot_canvas.draw()



	def Revert_WindowSize_toDefault(self, ytype, Ylim_box_object):
		ymin_current, ymax_current = self.carrier_concentration_plot_drawing.get_ylim()
		if ytype == "YMin":
			Ylim_box_object.setText(str(ymin_current))
		if ytype == "YMax":
			Ylim_box_object.setText(str(ymax_current))
	


	def Organize_DOS_Data(self):
		
		# Initialize data
		energy = []
		gE = []
		
		# Orgnize data
		for energy_dos in sorted(np.asarray([float(i) for i in self.dos_data["DOS"].keys()])):
			energy.append(energy_dos)
			gE.append(float(self.dos_data["DOS"][str(energy_dos)]))
		
		# Store into global variables
		self.energy = np.asarray(energy)
		self.gE = np.asarray(gE)
	

	
	def Extract_Relevant_Energies_DOSs(self):
		
		# Initialize data
		energies_ValenceBand = []
		gE_ValenceBand = []
		energies_ConductionBand = []
		gE_ConductionBand = []
		
		# The DOS band gap may not be the band gap for the defect formation energy
		#	diagram, especially when band gap corrections are applied. To mitigate
		#	this problem, we use a scissor operator where the VBM and CBM in the
		#	DOSCAR file are repositioned to the band gap of the defect formation
		#	energy diagram.
		past_dos_bandgap = False
		for energy, gE in zip(self.energy, self.gE):
			if energy <= 0.0:
				# Get all energies and corresponding DOSs below VBM
				energies_ValenceBand.append(energy)
				gE_ValenceBand.append(gE)
			else:
				# Get all energies and corresponding DOSs above CBM
				if (not past_dos_bandgap) and (gE <= 1E-4):
					continue
				else:
					past_dos_bandgap = True
					energies_ConductionBand.append(energy)
					gE_ConductionBand.append(gE)
		
		# Make data into numpy arrays
		self.energies_ValenceBand = np.asarray(energies_ValenceBand)
		self.gE_ValenceBand = np.asarray(gE_ValenceBand)
		self.energies_ConductionBand = np.asarray(energies_ConductionBand)
		self.gE_ConductionBand = np.asarray(gE_ConductionBand)

		# Reposition band edges to corrected values (NOT ZERO-ED)
		self.energies_ValenceBand += self.EVBM
		self.energies_ConductionBand += self.ECBM - np.min(self.energies_ConductionBand)

		# Normalize DOS to be per volume
		self.gE_ValenceBand /= self.dos_data["Volume"]
		self.gE_ConductionBand /= self.dos_data["Volume"]



	def Update_DOS_Data(self, dos_option_box, effmass_input_holes, effmass_input_electrons):

		# Check to see if the DOS option changed
		if dos_option_box.currentText() == self.dos_type:
			return

		# Set type of DOS used (DFT or parabolic assumption)
		self.dos_type = dos_option_box.currentText()

		if self.dos_type == "Real DOS":

			effmass_input_holes.setEnabled(False)
			effmass_input_electrons.setEnabled(False)

			self.Extract_Relevant_Energies_DOSs()
			
			self.Recalculate_CarrierConcentrations()

		elif self.dos_type == "Parabolic Approx.":

			# Enable effective mass inputs
			effmass_input_holes.setEnabled(True)
			effmass_input_electrons.setEnabled(True)

			#self.Update_CarrierConcentrations_EffMasses(effmass_input_holes, effmass_input_electrons)
			self.Update_EffMass(effmass_input_holes, effmass_input_electrons, override=True)

		else:

			# Exit the software if neither (for developers)
			sys.exit("DOS (dos_type) must be either 'Real DOS' or 'Parabolix Approx.'. Exiting... ")


	def Check_EffMass_Input(self, carrier_type, effmass_input_object, override):
		
		if override:
			return True

		# Check to make sure effective mass input is a float
		try:
			float(effmass_input_object.text())
		except:
			effmass_input_object.setText("1.00")
			pass

		# Check to see if effective mass input has changed. If not, don't proceed
		if self.effmass_dict[carrier_type] == float(effmass_input_object.text()):
			return False

		# Check to make sure that a plot has been generated
		if self.carrier_concentration_intrinsic_defect_hole_plot is None:
			effmass_input_object.setText("1.00")
			return False

		# Check to make sure that the effective mass is positive
		if float(effmass_input_object.text()) <= 0.0:
			effmass_input_object.setText(str(self.effmass_dict[carrier_type]))
			return False

		return True



	def Update_EffMass(self, effmass_input_holes, effmass_input_electrons, override=False):

		# Check whether we should recalculate carrier concentrations
		effmass_check_holes = self.Check_EffMass_Input("holes", effmass_input_holes, override)
		effmass_check_electrons = self.Check_EffMass_Input("electrons", effmass_input_electrons, override)

		# If any check fails, don't proceed (avoid recalculation)
		if (not effmass_check_holes) and (not effmass_check_electrons):
			return

		self.Update_CarrierConcentrations_EffMasses(effmass_input_holes, effmass_input_electrons)



	def Update_CarrierConcentrations_EffMasses(self, effmass_input_holes, effmass_input_electrons):

		# Update effective masses
		self.effmass_dict["holes"] = float(effmass_input_holes.text())
		self.effmass_dict["electrons"] = float(effmass_input_electrons.text())

		# Update DOSs
		me = 9.109E-31
		e = 1.602E-19
		h = 6.626E-34
		hbar = h / 2 / np.pi
		self.gE_ValenceBand = np.sqrt(2) * (self.effmass_dict["holes"] * me)**(3./2) / np.pi**2 / hbar**3 * np.sqrt((self.EVBM - self.energies_ValenceBand)*e) * 1E-6 * e
		self.gE_ConductionBand = np.sqrt(2) * (self.effmass_dict["electrons"] * me)**(3./2) / np.pi**2 / hbar**3 * np.sqrt((self.energies_ConductionBand - self.ECBM)*e) * 1E-6 * e

		self.Recalculate_CarrierConcentrations()



	def Recalculate_CarrierConcentrations(self):

		message_box = QMainWindow()
		message_box.setWindowTitle("Help")
		message_box.setWindowIcon(QIcon(vtandem_source_path+"/logo/LogoSmall.png"))
		message_box_widget = QWidget()
		message_box_widget_layout = QVBoxLayout(message_box_widget)
		message_box_message = QLabel("Recalculating carrier concentrations...")
		message_box_widget_layout.addWidget(message_box_message)
		message_box.setCentralWidget(message_box_widget)
		message_box.show()
		
		# Recalculate carrier concentrations
		self.Calculate_Hole_Electron_Concentration_Matrices()
		
		message_box.close()
		
		if self.carrier_concentration_intrinsic_defect_hole_plot is not None:
			self.Update_CarrierConcentration_Plot()



	def Update_DefectConcentration_DistributionStatistics(self, stat_box):

		# Update distribution statistics (Maxwell-Boltzmann or Fermi-Dirac)
		self.defectconc_stat = stat_box.currentText()

		# Update carrier concentrations (provided a plot exists)
		if self.carrier_concentration_intrinsic_defect_hole_plot is not None:
			self.Update_CarrierConcentration_Plot()


	
	def Update_Deltamus(self, deltamu_values):

		# Args:	
		# 	deltamu_values: Dictionary of deltamu values, in element:value pairs
		for element in self.mu_elements.keys():
			self.mu_elements[element]["deltamu"] = deltamu_values[element]



	# Free carrier concentrations are calculated separately from defect concentrations. This is to prevent
	#	having to calculate them repeatedly for different thermodynamic conditions (delta mu values) since
	#	they're the same in each condition.
	def Calculate_Hole_Electron_Concentration_Matrices(self):
		self.hole_concentrations_dict, self.electron_concentrations_dict = Calculate_FreeHole_FreeElectron_Concentrations(	self.temperature_array, \
																															self.fermi_energy_array, \
																															self.gE_ValenceBand, \
																															self.energies_ValenceBand, \
																															self.gE_ConductionBand, \
																															self.energies_ConductionBand )
	
	

	def Initialize_CarrierConcentration_Plot(self):
		
		intrinsic_defect_hole_concentration, intrinsic_defect_electron_concentration, total_hole_concentration, total_electron_concentration, intrinsic_equilibrium_fermi_energy_temperature, total_equilibrium_fermi_energy_temperature = Calculate_CarrierConcentration(	EVBM = self.EVBM, \
																																																																			ECBM = self.ECBM, \
																																																																			energies_ValenceBand = self.energies_ValenceBand, \
																																																																			gE_ValenceBand = self.gE_ValenceBand, \
																																																																			energies_ConductionBand = self.energies_ConductionBand, \
																																																																			gE_ConductionBand = self.gE_ConductionBand, \
																																																																			defects_data = self.defects_data, \
																																																																			main_compound_info = self.main_compound_info, \
																																																																			mu_elements = self.mu_elements, \
																																																																			temperature_array = self.temperature_array, \
																																																																			fermi_energy_array = self.fermi_energy_array, \
																																																																			volume = self.vol, \
																																																																			extrinsic_defects = self.extrinsic_defects, \
																																																																			dopant = self.dopant, \
																																																																			dopant_mu0 = self.dopant_mu0, \
																																																																			dopant_deltamu = self.dopant_deltamu, \
																																																																			hole_concentrations_dict = self.hole_concentrations_dict, \
																																																																			electron_concentrations_dict = self.electron_concentrations_dict, \
																																																																			synthesis_temperature = self.synthesis_temperature, \
																																																																			defectconc_stat = self.defectconc_stat )
		
		# Update equilibrium Fermi energy
		self.intrinsic_equilibrium_fermi_energy = intrinsic_equilibrium_fermi_energy_temperature
		self.total_equilibrium_fermi_energy = total_equilibrium_fermi_energy_temperature
		
		try:
			self.carrier_concentration_intrinsic_defect_hole_plot.remove()
			self.carrier_concentration_total_hole_plot.remove()
		except:
			pass

		try:
			self.carrier_concentration_intrinsic_defect_electron_plot.remove()
			self.carrier_concentration_total_electron_plot.remove()
		except:
			pass
		
		self.carrier_concentration_intrinsic_defect_hole_plot, = self.carrier_concentration_plot_drawing.semilogy(self.temperature_array, intrinsic_defect_hole_concentration, 'o-', color='red', label='Hole')
		if self.dopant != "None":
			self.carrier_concentration_total_hole_plot, = self.carrier_concentration_plot_drawing.semilogy(self.temperature_array, total_hole_concentration, 'o-', markerfacecolor='none', markeredgecolor='red', color='red', ls='--', label='Hole (With Dopant)')
		
		self.carrier_concentration_intrinsic_defect_electron_plot, = self.carrier_concentration_plot_drawing.semilogy(self.temperature_array, intrinsic_defect_electron_concentration, 'o-', color='green', label='Electron')
		if self.dopant != "None":
			self.carrier_concentration_total_electron_plot, = self.carrier_concentration_plot_drawing.semilogy(self.temperature_array, total_electron_concentration, 'o-', markerfacecolor='none', markeredgecolor='green', color='green', ls='--', label='Electron (With Dopant)')
		
		self.carrier_concentration_plot_drawing.legend(loc=1, fontsize=self.font['size'])
		self.carrier_concentration_plot_canvas.draw()
	


	def Update_CarrierConcentration_Plot(self):
		
		intrinsic_defect_hole_concentration, intrinsic_defect_electron_concentration, total_hole_concentration, total_electron_concentration, intrinsic_equilibrium_fermi_energy_temperature, total_equilibrium_fermi_energy_temperature = Calculate_CarrierConcentration(	EVBM = self.EVBM, \
																																																																			ECBM = self.ECBM, \
																																																																			energies_ValenceBand = self.energies_ValenceBand, \
																																																																			gE_ValenceBand = self.gE_ValenceBand, \
																																																																			energies_ConductionBand = self.energies_ConductionBand, \
																																																																			gE_ConductionBand = self.gE_ConductionBand, \
																																																																			defects_data = self.defects_data, \
																																																																			main_compound_info = self.main_compound_info, \
																																																																			mu_elements = self.mu_elements, \
																																																																			temperature_array = self.temperature_array, \
																																																																			fermi_energy_array = self.fermi_energy_array, \
																																																																			volume = self.vol, \
																																																																			extrinsic_defects = self.extrinsic_defects, \
																																																																			dopant = self.dopant, \
																																																																			dopant_mu0 = self.dopant_mu0, \
																																																																			dopant_deltamu = self.dopant_deltamu, \
																																																																			hole_concentrations_dict = self.hole_concentrations_dict, \
																																																																			electron_concentrations_dict = self.electron_concentrations_dict, \
																																																																			synthesis_temperature = self.synthesis_temperature, \
																																																																			defectconc_stat = self.defectconc_stat )
		
		# Update equilibrium Fermi energy
		self.intrinsic_equilibrium_fermi_energy = intrinsic_equilibrium_fermi_energy_temperature
		self.total_equilibrium_fermi_energy = total_equilibrium_fermi_energy_temperature
		
		self.carrier_concentration_intrinsic_defect_hole_plot.set_ydata(intrinsic_defect_hole_concentration)
		if self.dopant != "None":
			self.carrier_concentration_total_hole_plot.set_ydata(total_hole_concentration)
		
		self.carrier_concentration_intrinsic_defect_electron_plot.set_ydata(intrinsic_defect_electron_concentration)
		if self.dopant != "None":
			self.carrier_concentration_total_electron_plot.set_ydata(total_electron_concentration)
		
		self.carrier_concentration_plot_canvas.draw()

