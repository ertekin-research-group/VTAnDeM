
__author__ = 'Michael_Lidia_Jiaxing_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'

import PyQt5
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from PyQt5.QtGui import *

###############################################################################################
###################################### Save Figure ############################################
###############################################################################################
	
class SaveFigure(object):
	
	def __init__(self, plot_figure_object):
		
		self.plot_figure_object = plot_figure_object
	
	def SaveFigure(self):
		
		options = QFileDialog.Options()
		options |= QFileDialog.DontUseNativeDialog
		filename, extension_type = QFileDialog.getSaveFileName(caption = "Save Figure", filter = "Portable Network Graphics (*.png);;" \
																								+"Portable Document Format (*.pdf);;" \
																								+"Scalable Vector Graphics (*.svg);;" \
																								+"Encapsulated PostScript (*.eps)", options=options)
		
		if filename:
			extension = extension_type.split(".")[-1].split(")")[0]
			if filename.split(".")[-1] == extension:
				self.plot_figure_object.savefig(filename, bbox_inches='tight')
			else:
				self.plot_figure_object.savefig(filename+"."+extension, bbox_inches='tight')

