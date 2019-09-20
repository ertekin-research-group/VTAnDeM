
__author__ = 'Michael_Lidia_Jiaxing_Benita_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'


###############################################################################################################################
###############################################################################################################################
##################################################### Introduction ############################################################
###############################################################################################################################
###############################################################################################################################

# 	This script stores the phase stability diagram of a single ternary compound as an object within VTAnDeM.
#
# 	The stability of any compound is dictated by thermodynamic principles. All materials are made of atoms, and multicomponent
#		systems in particular (e.g. systems of two or more types of atoms, like GaAs) can exhibit distinct phases/materials
#		depending on how the atoms are arranged. Therefore, for a given set of atoms, there is a myriad of possible materials 
#		that can be formed by e.g. different synthesis methods. Fundamentally, the material with the lowest formation energy in some
#		environmental condition (e.g. GaAs in a 'Ga-rich' environment) is most and will likely form in a lab. Accordingly, if
#		we are studying some ternary system, there may be other compounds that "compete" with the main compound to become the
#		most stable. For example, if we are interested in studying the ternary compound Hg2GeTe4, in regions that are Cu-rich
#		or Hg-poor (or any combination of the sorts), another compound such as HgTe may be more stable than Hg2GeTe4. It is 
#		therefore crucial for material/device processing purposes that we know the exact environmental conditions necessary to create a 
#		stable version of the ternary compound of interest. The stability of compounds can fortunately be studied using Density 
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



class ChemicalPotential_Ternary_PhaseDiagram3D(object):
	
	def __init__(self, parent = None, main_compound = None, first_element = None, second_element = None, third_element = None):
		
		# All elements in the periodic table
		self.all_elements = []
		for element in periodictable.elements:
			self.all_elements.append(str(element))
		
		# Font description for phase stability diagram plot
		self.font = {'family': 'sans-serif',
				'color':  'black',
				'weight': 'normal',
				'size': 14 }
		
		# Establish the first, second, and third elements along with the dependent element of the ternary compound.
		self.main_compound  = main_compound
		self.element_x = first_element
		self.element_y = second_element
		self.dependent_element = third_element
		self.elements_list = [self.element_x, self.element_y, self.dependent_element]	# List of elements that gets updated as user selects element order
		self.elements_list_original = [self.element_x, self.element_y, self.dependent_element]	# Unchanged list of elements for compound naming purposes (e.g. see Compound_Name_Formal)
		
		
		# Store all extracted DFT data
		self.compounds_info = {}
		
		# Store colors of all competing compounds
		self.competing_compounds_colorwheel = {}
		
		# 3D ternary phase diagram plot objects
		self.ternary_phasediagram_3d_plot_figure = plt.figure()
		self.ternary_phasediagram_3d_plot_canvas = FigureCanvas(self.ternary_phasediagram_3d_plot_figure)
		self.ternary_phasediagram_3d_plot_axes = Axes3D(self.ternary_phasediagram_3d_plot_figure)		# NOTE: Axes3D MUST be called AFTER calling FigureCanvas
		
		self.path = None
		
		self.ternary_phasediagram_3d_plot_axes.view_init(elev=15, azim=145)
	
	
	
	###############################################################################################
	########################## Rewrite Compound Name Latex-Style ##################################
	###############################################################################################
	
	def Compound_Name_Formal(self, compound_name):
		
		if compound_name in self.all_elements:
			return compound_name
		
		# Go into the compounds_info dictionary and obtains the chemistry and stoichiometry of the compound of choice
		compound_species_info = self.compounds_info[compound_name]
		
		compound_name_formal = ""
		for species in compound_species_info["elements_list"]:				# Loop through the list of possible species that can be contained in the compound
			if compound_species_info[species] == 0:
				continue								# Don't add the species to the name if the compound doesn't contain the species
			elif compound_species_info[species] == 1:
				compound_name_formal += species			# Add the species to the name of the compound
			elif compound_species_info[species] > 1:
				compound_name_formal += species+"$_{"+str(int(compound_species_info[species]))+"}$"	# Add the species to the name of the compound
				#compound_name_formal += species+"<sub>"+str(int(compound_species_info[species]))+"</sub>"	# Add the species to the name of the compound
																											#	with a subscript for the stoichiometry
		
		return compound_name_formal
	
	
	
	###############################################################################################
	################################## Set the dependent element ##################################
	###############################################################################################
	
	def Set_Elements(self, element_x, element_y, element_z):
		
		self.element_x = element_x
		self.element_y = element_y
		self.element_z = element_z
		self.elements_list = [self.element_x, self.element_y, self.element_z]
	
	
	
	###############################################################################################
	############################ Create 3D phase diagram inequalities #############################
	###############################################################################################
	
	def Obtain_Ternary_PhaseDiagram3D_Inequalities(self):
		
		A_matrix = []
		b_vector = []
		compounds_list = []
		
		main_compound_x_reduced = float(self.compounds_info[self.main_compound][self.element_x]) / float(self.compounds_info[self.main_compound][self.dependent_element])
		main_compound_y_reduced = float(self.compounds_info[self.main_compound][self.element_y]) / float(self.compounds_info[self.main_compound][self.dependent_element])
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
				competing_compound_dependent_element = float(self.compounds_info[competing_compound][self.dependent_element])
			except:
				competing_compound_dependent_element = 0.0
				pass
			
			coefficient_x = competing_compound_element_x - competing_compound_dependent_element*main_compound_x_reduced
			coefficient_y = competing_compound_element_y - competing_compound_dependent_element*main_compound_y_reduced
			b = float(self.compounds_info[competing_compound]["enthalpy"]) - competing_compound_dependent_element*main_compound_enthalpy_reduced
			
			A_matrix.append(np.array([coefficient_x, coefficient_y]))
			b_vector.append(b)
			compounds_list.append(self.Compound_Name_Formal(competing_compound))
		
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
			
			print(ininc_i)
			
			if len(ininc_i) < 2:
				continue
			ininc_i = self.Sort_Plane_Vertices(ininc_i, adj)
			
			print(ininc_i)
			
			x = []
			y = []
			z = []
			for v in ininc_i:
				x.append(verts[v][0])
				y.append(verts[v][1])
				#z.append(verts[v][2])
				z_point = ( self.compounds_info[self.main_compound]["enthalpy"] - self.compounds_info[self.main_compound][self.element_x]*verts[v][0] - self.compounds_info[self.main_compound][self.element_y]*verts[v][1] ) / float(self.compounds_info[self.main_compound][self.dependent_element])
				z.append(z_point)
			"""
			x.append(verts[ininc_i[0]][0])
			y.append(verts[ininc_i[0]][1])
			#z.append(verts[ininc_i[0]][2])
			z_point = ( self.compounds_info[self.main_compound]["enthalpy"] - self.compounds_info[self.main_compound][self.element_x]*verts[v][0] - self.compounds_info[self.main_compound][self.element_y]*verts[v][1] ) / float(self.compounds_info[self.main_compound][self.dependent_element])
			z.append(z_point)
			"""
			coord = [list(zip(x, y, z))]
			
			print(coord)
			
			label = sec_phases[i]
			polygon = Poly3DCollection(coord, alpha=0.9, label=label, closed=True)
			polygon.set_facecolor(cmap(color_counter))
			self.competing_compounds_colorwheel[label] = cmap(color_counter)
			color_counter += 1
			
			polygon._facecolors2d=polygon._facecolors3d
			polygon._edgecolors2d=polygon._edgecolors3d
			
			self.ternary_phasediagram_3d_plot_axes.add_collection3d(polygon)
			path = Line3DCollection(coord, lw=2, color=self.competing_compounds_colorwheel[label])
			self.ternary_phasediagram_3d_plot_axes.add_collection3d(path)
		
		# Include competing compounds not included in phase stability bound into colorwheel
		for competing_compound in sec_phases:
			if competing_compound not in self.competing_compounds_colorwheel.keys():
				self.competing_compounds_colorwheel[competing_compound] = cmap(color_counter)
				color_counter += 1
	
	
	def Activate_PhaseDiagram3D_Plot_Axes(self):
		# calc elemental chemical potential and set them 0
		buffer = 0.2 # eV
		
		# Endpoints of phase diagram
		main_compound_enthalpy = self.compounds_info[self.main_compound]["enthalpy"]
		main_compound_number_first_specie = self.compounds_info[self.main_compound][self.element_x]
		main_compound_number_second_specie = self.compounds_info[self.main_compound][self.element_y]
		main_compound_number_third_specie = self.compounds_info[self.main_compound][self.dependent_element]
		phasediagram_endpoints = min(main_compound_enthalpy/main_compound_number_first_specie, main_compound_enthalpy/main_compound_number_second_specie, main_compound_enthalpy/main_compound_number_third_specie)
		
		self.ternary_phasediagram_3d_plot_axes.set_xlim([phasediagram_endpoints - buffer, 0 + buffer])
		self.ternary_phasediagram_3d_plot_axes.set_ylim([phasediagram_endpoints - buffer, 0 + buffer])
		self.ternary_phasediagram_3d_plot_axes.set_zlim([phasediagram_endpoints - buffer, 0 + buffer])
		
		self.ternary_phasediagram_3d_plot_axes.set_xlabel(r"$\Delta\mu_{"+self.element_x+"}$ (eV)", fontdict=self.font)  
		self.ternary_phasediagram_3d_plot_axes.set_ylabel(r"$\Delta\mu_{"+self.element_y+"}$ (eV)", fontdict=self.font)  
		self.ternary_phasediagram_3d_plot_axes.set_zlabel(r"$\Delta\mu_{"+self.element_z+"}$ (eV)", fontdict=self.font)
		
		
		self.ternary_phasediagram_3d_plot_axes.set_aspect("auto")
		
		self.ternary_phasediagram_3d_plot_axes.legend(loc='center left')
		self.ternary_phasediagram_3d_plot_axes.legend(loc='upper center', ncol=5)
	
	
	
	
	
	def Draw_PhaseDiagram3D(self):
		
		A_matrix, b_vector, compounds_list = self.Obtain_Ternary_PhaseDiagram3D_Inequalities()
		
		phasediagram_polyhedron = Hrep(A_matrix, b_vector)	# H-representation of Ax <= b
		
		print(phasediagram_polyhedron)
		
		self.Draw_PhaseDiagram_Planes(phasediagram_polyhedron.generators, phasediagram_polyhedron.ininc, phasediagram_polyhedron.adj, compounds_list)
		
		self.Activate_PhaseDiagram3D_Plot_Axes()
		
		self.ternary_phasediagram_3d_plot_canvas.draw()
	
	
	
	
	
	
	
	
	
	
	def Animate_SpacePotato(self, i):
		
		self.ternary_phasediagram_3d_plot_axes.clear()
		
		self.Draw_PhaseDiagram3D()
		
		self.ternary_phasediagram_3d_plot_axes.view_init(elev=15, azim=i)
		
		print(i)
		self.spacepotato_animation_status.setValue(100./60.*i)
		
		return self.ternary_phasediagram_3d_plot_figure,
	
	
	
	def Make_SpacePotato_Rotation_Animation(self):
		
		self.spacepotato_animation_generator_window = QMainWindow()
		
		self.spacepotato_animation_generator_widget = QWidget()
		self.spacepotato_animation_generator_widget_layout = QVBoxLayout(self.spacepotato_animation_generator_widget)
		self.spacepotato_animation_generator_window.setCentralWidget(self.spacepotato_animation_generator_widget)
		
		self.spacepotato_animation_status = QProgressBar()
		self.spacepotato_animation_status.setValue(0.0)
		self.spacepotato_animation_generator_widget_layout.addWidget(self.spacepotato_animation_status)
		
		self.spacepotato_animation_generator_window.show()
		
		anim = animation.FuncAnimation(self.ternary_phasediagram_3d_plot_figure, self.Animate_SpacePotato, frames=360, blit=True)
		anim.save('rotating_Cu2HgGeTe4.gif', writer='imagemagick', fps=30)
		
		self.spacepotato_animation_generator_window.close()








