
__author__ = 'Michael_Lidia_Jiaxing_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'


###############################################################################################################################
###############################################################################################################################
##################################################### Introduction ############################################################
###############################################################################################################################
###############################################################################################################################

# 	This script stores the phase stability diagram of a single quaternary compound as an object within VTAnDeM.
#
# 	The stability of any compound is dictated by thermodynamic principles. All materials are made of atoms, and multicomponent
#		systems in particular (e.g. systems of two or more types of atoms, like GaAs) can exhibit distinct phases/materials
#		depending on how the atoms are arranged. Therefore, for a given set of atoms, there is a myriad of possible materials 
#		that can be formed by e.g. different synthesis methods. Fundamentally, the material with the lowest formation energy in some
#		environmental condition (e.g. GaAs in a 'Ga-rich' environment) is most and will likely form in a lab. Accordingly, if
#		we are studying some quaternary system, there may be other compounds that "compete" with the main compound to become the
#		most stable. For example, if we are interested in studying the quaternary compound Cu2HgGeTe4, in regions that are Cu-rich
#		or Hg-poor (or any combination of the sorts), another compound such as CuTe may be more stable than Cu2HgGeTe4. It is 
#		therefore crucial for material/device processing purposes that we know the exact environmental conditions necessary to create a 
#		stable version of the quaternary compound of interest. The stability of compounds can fortunately be studied using Density 
#		Functional Theory (DFT) calculations.
#	REWRITE ALL OF THIS


###############################################################################################################################
###############################################################################################################################
################################################### Import Libraries ##########################################################
###############################################################################################################################
###############################################################################################################################

import numpy as np
import periodictable
import matplotlib.pyplot as plt
from matplotlib import animation
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Path3DCollection, Line3DCollection
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from polyhedron import Hrep

from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

from vtandem.visualization.compound_name import Compound_Name_Formal



class ChemicalPotential_Quaternary_PhaseDiagram3D(object):
	
	def __init__(self, parent = None, main_compound = None, first_element = None, second_element = None, third_element = None, fourth_element = None):
		
		# All elements in the periodic table
		self.all_elements = []
		for element in periodictable.elements:
			self.all_elements.append(str(element))
		
		# Font description for phase stability diagram plot
		self.font = {'family': 'sans-serif',
				'color':  'black',
				'weight': 'normal',
				'size': 14 }
		
		# Establish the first, second, and third elements along with the dependent element of the quaternary compound.
		self.main_compound  = main_compound
		self.element_x = first_element
		self.element_y = second_element
		self.element_z = third_element
		self.dependent_element = fourth_element
		self.elements_list = [self.element_x, self.element_y, self.element_z, self.dependent_element]	# List of elements that gets updated as user selects element order
		self.elements_list_original = [self.element_x, self.element_y, self.element_z, self.dependent_element]	# Unchanged list of elements for compound naming purposes (e.g. see Compound_Name_Formal)
		
		
		# Store all extracted DFT data
		self.compounds_info = {}
		
		# Track chemical potential of fourth element
		self.mu4 = 0.0
		
		# Store colors of all competing compounds
		self.competing_compounds_colorwheel = {}
		
		# 3D quaternary phase diagram plot objects
		self.chemicalpotential_phasediagram_plot_figure = plt.figure()
		self.chemicalpotential_phasediagram_plot_canvas = FigureCanvas(self.chemicalpotential_phasediagram_plot_figure)
		self.chemicalpotential_phasediagram_plot_axes = Axes3D(self.chemicalpotential_phasediagram_plot_figure)		# NOTE: Axes3D MUST be called AFTER calling FigureCanvas
		
		
		self.path = None
		
		self.chemicalpotential_phasediagram_plot_axes.view_init(elev=15, azim=145)
	
	
	
	###############################################################################################
	################################## Set the dependent element ##################################
	###############################################################################################
	
	def Set_Elements(self, element_x, element_y, element_z, dependent_element):
		
		self.element_x = element_x
		self.element_y = element_y
		self.element_z = element_z
		self.dependent_element = dependent_element
		self.elements_list = [self.element_x, self.element_y, self.element_z, self.dependent_element]
	
	
	
	###############################################################################################
	############################ Create 3D phase diagram inequalities #############################
	###############################################################################################
	
	def Obtain_Quaternary_PhaseDiagram3D_Inequalities(self):
		
		A_matrix = []
		b_vector = []
		compounds_list = []
		
		main_compound_x_reduced = float(self.compounds_info[self.main_compound][self.element_x]) / float(self.compounds_info[self.main_compound][self.dependent_element])
		main_compound_y_reduced = float(self.compounds_info[self.main_compound][self.element_y]) / float(self.compounds_info[self.main_compound][self.dependent_element])
		main_compound_z_reduced = float(self.compounds_info[self.main_compound][self.element_z]) / float(self.compounds_info[self.main_compound][self.dependent_element])
		main_compound_enthalpy_reduced = float(self.compounds_info[self.main_compound]["enthalpy"]) / float(self.compounds_info[self.main_compound][self.dependent_element])
		
		
		for competing_compound in self.compounds_info.keys():
			
			# Check whether compound is "competing"
			if (not set(self.compounds_info[competing_compound]["elements_list"]).issubset(self.elements_list_original)) or (competing_compound == self.main_compound):
				continue
			
			try:
				competing_compound_element_x = float(self.compounds_info[competing_compound][self.element_x])
			except:
				competing_compound_element_x = 0.0
				pass
			
			try:
				competing_compound_element_y = float(self.compounds_info[competing_compound][self.element_y])
			except:
				competing_compound_element_y = 0.0
				pass
			
			try:
				competing_compound_element_z = float(self.compounds_info[competing_compound][self.element_z])
			except:
				competing_compound_element_z = 0.0
				pass
			
			try:
				competing_compound_dependent_element = float(self.compounds_info[competing_compound][self.dependent_element])
			except:
				competing_compound_dependent_element = 0.0
				pass
			
			coefficient_x = competing_compound_element_x - competing_compound_dependent_element*main_compound_x_reduced
			coefficient_y = competing_compound_element_y - competing_compound_dependent_element*main_compound_y_reduced
			coefficient_z = competing_compound_element_z - competing_compound_dependent_element*main_compound_z_reduced
			b = float(self.compounds_info[competing_compound]["enthalpy"]) - competing_compound_dependent_element*main_compound_enthalpy_reduced
			
			A_matrix.append(np.array([coefficient_x, coefficient_y, coefficient_z]))
			b_vector.append(b)
			#compounds_list.append(self.Compound_Name_Formal(competing_compound))
			compounds_list.append(competing_compound)
		
		return A_matrix, b_vector, compounds_list
	
	
	
	def Sort_Plane_Vertices(self, ininc_i, adj):
		ininc_sorted = []
		ininc_i = list(ininc_i)
		while len(ininc_i) > 0:
			v = ininc_i.pop()
			ininc_sorted.append(v)
			# find adj
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
			if len(ininc_i) < 3:
				continue
			ininc_i = self.Sort_Plane_Vertices(ininc_i, adj)
			x = []
			y = []
			z = []
			for v in ininc_i:
				x.append(verts[v][0])
				y.append(verts[v][1])
				z.append(verts[v][2])
			x.append(verts[ininc_i[0]][0])
			y.append(verts[ininc_i[0]][1])
			z.append(verts[ininc_i[0]][2])
			coord = [list(zip(x, y, z))]
			
			label = sec_phases[i]
			polygon = Poly3DCollection(coord, alpha=0.9, label=Compound_Name_Formal(label, self.compounds_info, "latex"), closed=True)
			polygon.set_facecolor(cmap(color_counter))
			self.competing_compounds_colorwheel[Compound_Name_Formal(label, self.compounds_info, "latex")] = cmap(color_counter)
			color_counter += 1
			
			polygon._facecolors2d=polygon._facecolors3d
			polygon._edgecolors2d=polygon._edgecolors3d
			
			self.chemicalpotential_phasediagram_plot_axes.add_collection3d(polygon)
			path = Line3DCollection(coord, lw=2, color='k')
			self.chemicalpotential_phasediagram_plot_axes.add_collection3d(path)
		
		# Include competing compounds not included in phase stability bound into colorwheel
		for competing_compound in sec_phases:
			if competing_compound not in self.competing_compounds_colorwheel.keys():
				self.competing_compounds_colorwheel[Compound_Name_Formal(competing_compound, self.compounds_info, "latex")] = cmap(color_counter)
				color_counter += 1
	
	
	def Activate_PhaseDiagram3D_Plot_Axes(self):
		# calc elemental chemical potential and set them 0
		buffer = 0.2 # eV
		
		# Endpoints of phase diagram
		main_compound_enthalpy = self.compounds_info[self.main_compound]["enthalpy"]
		main_compound_number_first_specie = self.compounds_info[self.main_compound][self.element_x]
		main_compound_number_second_specie = self.compounds_info[self.main_compound][self.element_y]
		main_compound_number_third_specie = self.compounds_info[self.main_compound][self.element_z]
		main_compound_number_fourth_specie = self.compounds_info[self.main_compound][self.dependent_element]
		phasediagram_endpoints = min(main_compound_enthalpy/main_compound_number_first_specie, main_compound_enthalpy/main_compound_number_second_specie, main_compound_enthalpy/main_compound_number_third_specie, main_compound_enthalpy/main_compound_number_fourth_specie)
		
		self.chemicalpotential_phasediagram_plot_axes.set_xlim([phasediagram_endpoints - buffer, 0 + buffer])
		self.chemicalpotential_phasediagram_plot_axes.set_ylim([phasediagram_endpoints - buffer, 0 + buffer])
		self.chemicalpotential_phasediagram_plot_axes.set_zlim([phasediagram_endpoints - buffer, 0 + buffer])
		
		self.chemicalpotential_phasediagram_plot_axes.set_xlabel(r"$\Delta\mu_{"+self.element_x+"}$ (eV)", fontdict=self.font)  
		self.chemicalpotential_phasediagram_plot_axes.set_ylabel(r"$\Delta\mu_{"+self.element_y+"}$ (eV)", fontdict=self.font)  
		self.chemicalpotential_phasediagram_plot_axes.set_zlabel(r"$\Delta\mu_{"+self.element_z+"}$ (eV)", fontdict=self.font)
		
		
		self.chemicalpotential_phasediagram_plot_axes.set_aspect("auto")
		
		self.chemicalpotential_phasediagram_plot_axes.legend(loc='center left')
		self.chemicalpotential_phasediagram_plot_axes.legend(loc='upper center', ncol=5)
	
	
	
	def Draw_PhaseDiagram3D(self):
		
		A_matrix, b_vector, compounds_list = self.Obtain_Quaternary_PhaseDiagram3D_Inequalities()
		
		phasediagram_polyhedron = Hrep(A_matrix, b_vector)	# H-representation of Ax <= b
		
		self.Draw_PhaseDiagram_Planes(phasediagram_polyhedron.generators, phasediagram_polyhedron.ininc, phasediagram_polyhedron.adj, compounds_list)
		
		self.Activate_PhaseDiagram3D_Plot_Axes()
		
		self.chemicalpotential_phasediagram_plot_canvas.draw()
	
	
	
	
	
	
	def Obtain_Mu4_Outline_Inequalities(self):
		
		A_matrix = []
		b_vector = []
		
		main_compound_coefficient1 = self.compounds_info[self.main_compound][self.element_x] / self.compounds_info[self.main_compound][self.element_z]
		main_compound_coefficient2 = self.compounds_info[self.main_compound][self.element_y] / self.compounds_info[self.main_compound][self.element_z]
		
		main_compound_enthalpy_reduced = (self.compounds_info[self.main_compound]["enthalpy"] - self.compounds_info[self.main_compound][self.dependent_element]*self.mu4) / self.compounds_info[self.main_compound][self.element_z]
		
		for competing_compound in self.compounds_info.keys():
			
			try:
				competing_compound_x = self.compounds_info[competing_compound][self.element_x]
			except:
				competing_compound_x = 0.0
				pass
			
			try:
				competing_compound_y = self.compounds_info[competing_compound][self.element_y]
			except:
				competing_compound_y = 0.0
				pass
			
			try:
				competing_compound_z = self.compounds_info[competing_compound][self.element_z]
			except:
				competing_compound_z = 0.0
				pass
			
			try:
				competing_compound_dependent = self.compounds_info[competing_compound][self.dependent_element]
			except:
				competing_compound_dependent = 0.0
				pass
			
			coefficient1 = competing_compound_x - competing_compound_z*main_compound_coefficient1
			coefficient2 = competing_compound_y - competing_compound_z*main_compound_coefficient2
			
			competing_compound_enthalpy_reduced = self.compounds_info[competing_compound]["enthalpy"] - competing_compound_dependent*self.mu4
			b = competing_compound_enthalpy_reduced - competing_compound_z*main_compound_enthalpy_reduced
			
			A_matrix.append(np.array([coefficient1, coefficient2]))
			b_vector.append(b)
		
		return A_matrix, b_vector
	
	
	
	def Draw_Mu4_Outline_Region(self, verts, ininc, adj):
		
		verts_z = (self.compounds_info[self.main_compound]["enthalpy"] - self.compounds_info[self.main_compound][self.dependent_element]*self.mu4 \
					- self.compounds_info[self.main_compound][self.element_x]*verts[:,0] \
					- self.compounds_info[self.main_compound][self.element_y]*verts[:,1] ) \
					/ self.compounds_info[self.main_compound][self.element_z]
		
		verts_z = verts_z.reshape((len(verts),1))
		plane_vertices = np.hstack((verts, verts_z))
		
		for i, ininc_i in enumerate(ininc):
			if len(ininc_i) < 3:
				continue
			ininc_i = self.Sort_Plane_Vertices(ininc_i, adj)
			x = []
			y = []
			z = []
			for v in ininc_i:
				x.append(plane_vertices[v][0])
				y.append(plane_vertices[v][1])
				z.append(plane_vertices[v][2])
			x.append(plane_vertices[ininc_i[0]][0])
			y.append(plane_vertices[ininc_i[0]][1])
			z.append(plane_vertices[ininc_i[0]][2])
			coord = [list(zip(x, y, z))]
			
			try:
				self.path.remove()
			except:
				pass
			
			self.path = Line3DCollection(coord, lw=4, color='k')
			self.chemicalpotential_phasediagram_plot_axes.add_collection3d(self.path)
	
	
	
	
	def Draw_Mu4_Outline(self):
		
		A_matrix, b_vector = self.Obtain_Mu4_Outline_Inequalities()
		
		phasediagram_projection2d = Hrep(A_matrix, b_vector)	# H-representation of Ax <= b
		
		self.Draw_Mu4_Outline_Region(phasediagram_projection2d.generators, phasediagram_projection2d.ininc, phasediagram_projection2d.adj)
		
		self.chemicalpotential_phasediagram_plot_canvas.draw()
	
	
	
	
	
	
	def Animate_SpacePotato(self, i):
		
		self.chemicalpotential_phasediagram_plot_axes.clear()
		
		self.Draw_PhaseDiagram3D()
		self.Draw_Mu4_Outline()
		
		self.chemicalpotential_phasediagram_plot_axes.view_init(elev=15, azim=i)
		
		print(i)
		self.spacepotato_animation_status.setValue(100./60.*i)
		
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
		anim.save('rotating_Cu2HgGeTe4.gif', writer='imagemagick', fps=30)
		
		self.spacepotato_animation_generator_window.close()








