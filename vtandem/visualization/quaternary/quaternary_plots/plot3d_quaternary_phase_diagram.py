
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

from vtandem.visualization.utils.compound_name import Compound_Name_Formal
from vtandem.visualization.plots.plot_chemicalpotential_phasediagram3d import Plot_ChemicalPotential_PhaseDiagram3D



class ChemicalPotential_Quaternary_PhaseDiagram3D(Plot_ChemicalPotential_PhaseDiagram3D):
	
	def __init__(self, parent = None, main_compound = None, first_element = None, second_element = None, third_element = None, fourth_element = None):
		
		super().__init__(type = "quaternary")
		
		# Establish the first, second, and third elements along with the dependent element of the quaternary compound.
		self.main_compound = main_compound
		self.element_x = first_element
		self.element_y = second_element
		self.element_z = third_element
		self.dependent_element = fourth_element
		self.elements_list = [self.element_x, self.element_y, self.element_z, self.dependent_element]	# List of elements that gets updated as user selects element order
		self.elements_list_original = [first_element, second_element, third_element, fourth_element]	# Unchanged list of elements for compound naming purposes (e.g. see Compound_Name_Formal)
		
		# Track chemical potential of fourth element
		self.mu4 = 0.0
	
	
	
	
	###############################################################################################
	################################## Set the dependent element ##################################
	###############################################################################################
	
	def Set_Elements(self, element_x, element_y, element_z, dependent_element):
		
		self.element_x = element_x
		self.element_y = element_y
		self.element_z = element_z
		self.dependent_element = dependent_element
		self.elements_list = [self.element_x, self.element_y, self.element_z, self.dependent_element]
	
	
	
	def Obtain_Mu4_Outline_Inequalities(self):
		
		A_matrix = []
		b_vector = []
		
		main_compound_coefficient1 = self.main_compound_info["dft_"+self.element_x] / self.main_compound_info["dft_"+self.element_z]
		main_compound_coefficient2 = self.main_compound_info["dft_"+self.element_y] / self.main_compound_info["dft_"+self.element_z]
		
		main_compound_enthalpy_reduced = (self.main_compound_enthalpy - self.main_compound_info["dft_"+self.dependent_element]*self.mu4) / self.main_compound_info["dft_"+self.element_z]
		
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
			
			coefficient1 = competing_compound_element_count[self.element_x] - competing_compound_element_count[self.element_z]*main_compound_coefficient1
			coefficient2 = competing_compound_element_count[self.element_y] - competing_compound_element_count[self.element_z]*main_compound_coefficient2
			competing_compound_enthalpy_reduced = competing_compound_enthalpy_tracker - competing_compound_element_count[self.dependent_element]*self.mu4
			b = competing_compound_enthalpy_reduced - competing_compound_element_count[self.element_z]*main_compound_enthalpy_reduced
			
			A_matrix.append(np.array([coefficient1, coefficient2]))
			b_vector.append(b)
		
		return A_matrix, b_vector
	
	
	
	def Draw_Mu4_Outline_Region(self, verts, ininc, adj):
		
		"""
		verts_z = (self.compounds_info[self.main_compound]["enthalpy"] - self.compounds_info[self.main_compound][self.dependent_element]*self.mu4 \
					- self.compounds_info[self.main_compound][self.element_x]*verts[:,0] \
					- self.compounds_info[self.main_compound][self.element_y]*verts[:,1] ) \
					/ self.compounds_info[self.main_compound][self.element_z]
		"""
		verts_z = (self.main_compound_enthalpy - self.main_compound_info["dft_"+self.dependent_element]*self.mu4 \
					- self.main_compound_info["dft_"+self.element_x]*verts[:,0] \
					- self.main_compound_info["dft_"+self.element_y]*verts[:,1] ) \
					/ self.main_compound_info["dft_"+self.element_z]
		
		
		verts_z = verts_z.reshape((len(verts),1))
		plane_vertices = np.hstack((verts, verts_z))
		
		print(ininc)
		
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








