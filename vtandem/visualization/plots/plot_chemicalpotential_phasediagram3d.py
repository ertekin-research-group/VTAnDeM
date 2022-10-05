
__author__ = 'Michael_Lidia_Jiaxing_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'

###############################################################################################################################
################################################### Import Libraries ##########################################################
###############################################################################################################################

import numpy as np
import periodictable
import matplotlib.pyplot as plt
from matplotlib import animation
#from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Path3DCollection, Line3DCollection
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from polyhedron import Hrep

from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

from vtandem.visualization.utils.compound_name import Compound_Name_Formal

from vtandem.visualization.plots.save_plot import SaveFigure


class Plot_ChemicalPotential_PhaseDiagram3D(SaveFigure):
	
	def __init__(self, type: str):
		
		# Record type of phase diagram (either ternary or quaternary)
		if type not in ["ternary", "quaternary"]:
			raise ValueError("Argument 'type' can only be either 'ternary' or 'quaternary'.")
		self.type = type
		
		# All elements in the periodic table
		self.all_elements = []
		for element in periodictable.elements:
			self.all_elements.append(str(element))
		
		# Font description for phase stability diagram plot
		self.font = {'color': 'black', 'weight': 'normal', 'size': 14 }
		
		# Store all extracted DFT data
		self.main_compound_info = {}
		self.compounds_info = {}
		self.main_compound_enthalpy = 0.0
		
		# Store colors of all competing compounds
		self.competing_compounds_colorwheel = {}
		
		# 3D phase diagram plot objects
		self.chemicalpotential_phasediagram_plot_figure = plt.figure()
		self.chemicalpotential_phasediagram_plot_canvas = FigureCanvas(self.chemicalpotential_phasediagram_plot_figure)
		#self.chemicalpotential_phasediagram_plot_axes = Axes3D(self.chemicalpotential_phasediagram_plot_figure)		# NOTE: Axes3D MUST be called AFTER calling FigureCanvas
		self.chemicalpotential_phasediagram_plot_axes = self.chemicalpotential_phasediagram_plot_figure.add_subplot(111, projection='3d')
		
		self.path = None
		
		self.chemicalpotential_phasediagram_plot_axes.view_init(elev=15, azim=145)
		
		print("Chemical potential 3d set up!")
		
		# Save figure feature
		SaveFigure.__init__(self, self.chemicalpotential_phasediagram_plot_figure)
	
	
	###############################################################################################
	############################ Create 3D phase diagram inequalities #############################
	###############################################################################################
	
	def Obtain_PhaseDiagram3D_Inequalities(self):
		
		A_matrix = []
		b_vector = []
		compounds_list = []
		
		main_compound_element_count_reduced = {}
		for element in self.elements_list[:-1]:	# Exclude dependent element
			main_compound_element_count_reduced[element] = float(self.main_compound_info["dft_"+element]) / float(self.main_compound_info["dft_"+self.dependent_element])
		main_compound_enthalpy_reduced = float(self.main_compound_enthalpy) / float(self.main_compound_info["dft_"+self.dependent_element])
		
		for competing_compound in self.compounds_info.keys():
			
			# Check whether compound is "competing"
			if (not set(self.compounds_info[competing_compound]["elements_list"]).issubset(self.elements_list_original)) or (competing_compound == self.main_compound):
				continue
			
			competing_compound_element_count = {}
			for element in self.elements_list:
				try:
					competing_compound_element_count[element] = float(self.compounds_info[competing_compound]["dft_"+element])
				except:
					competing_compound_element_count[element] = 0.0
			
			# Get enthalpy of competing compound
			competing_compound_enthalpy_tracker = self.compounds_info[competing_compound]["dft_total_energy"]
			for element in self.elements_list:
				competing_compound_enthalpy_tracker -= competing_compound_element_count[element] * self.compounds_info[element]["mu0"]
			
			coefficients = []
			for element in self.elements_list[:-1]:
				coefficients.append(competing_compound_element_count[element] - competing_compound_element_count[self.dependent_element]*main_compound_element_count_reduced[element])
			b = competing_compound_enthalpy_tracker - competing_compound_element_count[self.dependent_element]*main_compound_enthalpy_reduced
			
			A_matrix.append(np.asarray(coefficients))
			b_vector.append(b)
			compounds_list.append(competing_compound)
		
		return A_matrix, b_vector, compounds_list
	
	
	def Sort_Plane_Vertices(self, ininc_i, adj):
		
		ininc_sorted = []
		ininc_i = list(ininc_i)
		
		while len(ininc_i) > 0:
			
			v = ininc_i.pop()
			ininc_sorted.append(v)
			
			# Find adjacent
			adj_i = adj[v]
			ininc_i = sorted(ininc_i, reverse = True, key = lambda x: np.where(np.concatenate([adj_i, np.arange(1000)]) == x)[0][0])
		
		return ininc_sorted
	
	
	def Draw_PhaseDiagram_Planes(self, verts, ininc, adj, sec_phases=[]):
		
		# find number of polygons with verts 
		if len([half for half in ininc if len(half) > 3]) > 8:
			cmap = plt.get_cmap("tab10_r")
		else:
			cmap = plt.get_cmap("Set1")
		
		# draw plane
		color_counter = 0
		for i, ininc_i in enumerate(ininc):
			
			if len(ininc_i) < len(self.elements_list)-1:
				continue
			
			ininc_i = self.Sort_Plane_Vertices(ininc_i, adj)
			
			x = []
			y = []
			z = []
			if self.type == "ternary":
				for v in ininc_i:
					x.append(verts[v][0])
					y.append(verts[v][1])
					z_point = ( self.main_compound_enthalpy - self.main_compound_info["dft_"+self.element_x]*verts[v][0] - self.main_compound_info["dft_"+self.element_y]*verts[v][1] ) / float(self.main_compound_info["dft_"+self.dependent_element])
					z.append(z_point)
			elif self.type == "quaternary":
				for v in ininc_i:
					x.append(verts[v][0])
					y.append(verts[v][1])
					z.append(verts[v][2])
				x.append(verts[ininc_i[0]][0])
				y.append(verts[ininc_i[0]][1])
				z.append(verts[ininc_i[0]][2])
			
			coord = [list(zip(x, y, z))]
			
			label = sec_phases[i]
			polygon = Poly3DCollection(coord, alpha=0.7, label=Compound_Name_Formal(label, "latex"), closed=True)
			polygon.set_facecolor(cmap(color_counter))
			self.competing_compounds_colorwheel[Compound_Name_Formal(label, "unicode")] = cmap(color_counter)
			color_counter += 1
			
			polygon._facecolors2d = polygon._facecolor3d
			polygon._edgecolors2d = polygon._edgecolor3d
			
			
			self.chemicalpotential_phasediagram_plot_axes.add_collection3d(polygon)
			path = Line3DCollection(coord, lw=1, color=self.competing_compounds_colorwheel[Compound_Name_Formal(label, "unicode")])
			self.chemicalpotential_phasediagram_plot_axes.add_collection3d(path)
		
		# Include competing compounds not included in phase stability bound into colorwheel
		for competing_compound in sec_phases:
			label = Compound_Name_Formal(competing_compound, "unicode")
			if label not in self.competing_compounds_colorwheel.keys():
				self.competing_compounds_colorwheel[label] = cmap(color_counter)
				color_counter += 1
	

	
	def Activate_PhaseDiagram3D_Plot_Axes(self):
		
		# Endpoints of phase diagram
		endpoint_candidates = []
		for element in self.elements_list:
			endpoint_candidates.append( self.main_compound_enthalpy/self.main_compound_info["dft_"+element] )
		phasediagram_endpoints = min(endpoint_candidates)

		print("Hello")
		
		buffer = 0.2 # eV
		self.chemicalpotential_phasediagram_plot_axes.set_xlim([phasediagram_endpoints - buffer, 0 + buffer])
		self.chemicalpotential_phasediagram_plot_axes.set_ylim([phasediagram_endpoints - buffer, 0 + buffer])
		self.chemicalpotential_phasediagram_plot_axes.set_zlim([phasediagram_endpoints - buffer, 0 + buffer])
		
		self.chemicalpotential_phasediagram_plot_axes.set_xlabel(r"$\Delta\mu_{"+self.element_x+"}$ (eV)", fontdict=self.font)
		self.chemicalpotential_phasediagram_plot_axes.set_ylabel(r"$\Delta\mu_{"+self.element_y+"}$ (eV)", fontdict=self.font)
		self.chemicalpotential_phasediagram_plot_axes.set_zlabel(r"$\Delta\mu_{"+self.element_z+"}$ (eV)", fontdict=self.font)
		
		self.chemicalpotential_phasediagram_plot_axes.set_aspect("auto")
		
		self.chemicalpotential_phasediagram_plot_axes.legend(loc='center left', fontsize=self.font['size']-2)
		self.chemicalpotential_phasediagram_plot_axes.legend(loc='upper center', ncol=5, fontsize=self.font['size']-2)
	
	
	
	def Draw_PhaseDiagram3D(self):
		
		A_matrix, b_vector, compounds_list = self.Obtain_PhaseDiagram3D_Inequalities()
		phasediagram_polyhedron = Hrep(A_matrix, b_vector)	# H-representation of Ax <= b
		self.Draw_PhaseDiagram_Planes(phasediagram_polyhedron.generators, phasediagram_polyhedron.ininc, phasediagram_polyhedron.adj, compounds_list)
		
		self.Activate_PhaseDiagram3D_Plot_Axes()
		self.chemicalpotential_phasediagram_plot_canvas.draw()
	
	
	
	def Animate_SpacePotato(self, i):
		
		self.chemicalpotential_phasediagram_plot_axes.clear()
		
		self.Draw_PhaseDiagram3D()
		self.chemicalpotential_phasediagram_plot_axes.view_init(elev=15, azim=i)
		
		self.spacepotato_animation_status.setValue(100./60.*i)
		print("Working on Frame: ", i)
		
		return self.chemicalpotential_phasediagram_plot_figure,
	
	
	
	def Make_SpacePotato_Rotation_Animation(self):
		
		self.spacepotato_animation_generator_window = QMainWindow()
		
		self.spacepotato_animation_generator_widget = QWidget()
		self.spacepotato_animation_generator_widget_layout = QVBoxLayout(self.spacepotato_animation_generator_widget)
		self.spacepotato_animation_generator_window.setCentralWidget(self.spacepotato_animation_generator_widget)
		
		self.spacepotato_animation_status = QProgressBar()
		self.spacepotato_animation_status.setValue(0.0)
		self.spacepotato_animation_generator_widget_layout.addWidget(self.spacepotato_animation_status)
		
		self.spacepotato_animation_generator_window.show()
		
		anim = animation.FuncAnimation(self.chemicalpotential_phasediagram_plot_figure, self.Animate_SpacePotato, frames=360, blit=True)
		anim.save('animation.gif', writer='imagemagick', fps=30)
		print("Done animating!")
		
		self.spacepotato_animation_generator_window.close()


