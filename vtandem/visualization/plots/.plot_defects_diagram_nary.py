
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'
__author__ = 'Michael_Lidia_Jiaxing_Elif'

# Import defect formation energy diagram object
from vtandem.visualization.plots.plot_defects_diagram import Plot_DefectsDiagram

class Plot_DefectsDiagram_n_ary(Plot_DefectsDiagram):

	def __init__(self, elements_list):

		super().__init__()

		self.mu_elements = {}

		for element in elements_list:
			self.mu_elements[element] = {"mu0": 0.0, "deltamu": 0.0}
