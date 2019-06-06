
__author__ = 'Michael_Lidia_Jiaxing_Benita_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from labellines import labelLine, labelLines


class Quaternary_Defects_Diagram(object):
	
	
	def __init__(self, parent = None, main_compound = None, first_species = None, second_species = None, third_species = None, fourth_species = None):
		
		
		# Font description for phase stability diagram plot
		self.font = {'family': 'Arial',
				'color':  'black',
				'weight': 'normal',
				'size': 14 }
		
		# Establish the first, second, third, and fourth species of the quaternary compound.
		# Note that this list is subject to change, depending on what the user chooses.
		self.main_compound  = main_compound
		self.first_species	= first_species
		self.second_species	= second_species
		self.third_species	= third_species
		self.fourth_species	= fourth_species
		
		# Establish constants for the species
		self.species_a = first_species
		self.species_b = second_species
		self.species_c = third_species
		self.species_d = fourth_species
		
		# Store all extracted DFT data
		self.df0 = None
		
		# Store all mu values
		self.mu_values = {}
		self.mu_values[self.first_species]  = 0.0
		self.mu_values[self.second_species] = 0.0
		self.mu_values[self.third_species]  = 0.0
		self.mu_values[self.fourth_species] = 0.0
		
		
		self.ymin = -0.2
		self.ymax = 2.0
		self.EVBM = 0.0
		self.ECBM = 0.0
		self.defect_plots = {}
		self.defect_labels = {}
		
		
		
		### Defects diagram
		self.defects_diagram_plot_figure = plt.figure()
		self.defects_diagram_plot_figure.subplots_adjust(left=0.225)
		self.defects_diagram_plot_drawing = self.defects_diagram_plot_figure.add_subplot(111)
		self.defects_diagram_plot_canvas = FigureCanvas(self.defects_diagram_plot_figure)
	
	
	def Activate_DefectsDiagram_Plot_Axes(self):
		
		self.defects_diagram_plot_drawing.set_xlim(-1.0, 0.0)
		self.defects_diagram_plot_drawing.set_ylim(self.ymin, self.ymax)
		self.defects_diagram_plot_drawing.set_xlabel("Fermi Energy (eV)", fontdict=self.font, labelpad=12)
		self.defects_diagram_plot_drawing.set_ylabel("$\Delta$H(q,D) (eV)", fontdict=self.font, rotation=90, labelpad=20)
		self.defects_diagram_plot_drawing.xaxis.tick_bottom()
		self.defects_diagram_plot_drawing.yaxis.tick_left()
		self.defects_diagram_plot_drawing.tick_params(axis='both', labelsize=6)
		self.defects_diagram_plot_drawing.xaxis.set_label_position("bottom")
		self.defects_diagram_plot_drawing.yaxis.set_label_position("left")
		self.defects_diagram_plot_drawing.set_aspect("auto")
	
	
	
	def Plot_DefectsDiagram(self, xc):
		
		dmu_a = self.mu_values[self.species_a]
		dmu_b = self.mu_values[self.species_b]
		dmu_c = self.mu_values[self.species_c]
		dmu_d = self.mu_values[self.species_d]
		
		df = (self.df0[self.df0['Xc']==xc]).copy()
		df['class'] = ['$'+x+'^{'+str(y)+'}$' for x,y in zip(df['Defect'],df['q'])]
		df['En_prist'] = df['En_prist(1x1x1)']
		df['def_form_en'] = df['En(q,D)']-8*(df['En_prist'])-(df['n_a']*(df['mu_a']+float(dmu_a)) + df['n_b']*(df['mu_b']+float(dmu_b)) + df['n_c']*(df['mu_c']+float(dmu_c)) + df['n_d']*(df['mu_d']+float(dmu_d)))
		
		for idx in df.index:
			defect_name = df.ix[idx]['Defect'] + str(df.ix[idx]['q'])
			Ef_step = 0.02
			self.EVBM = df.ix[idx]['E_VBM']
			self.ECBM = df.ix[idx]['E_VBM'] + df.ix[idx]['Eg']
			self.Ef_pts = np.arange(self.EVBM, self.ECBM+0.1, Ef_step)
			self.E_enthalpy = (df.ix[idx])['def_form_en'] + ((df.ix[idx])['q']*self.Ef_pts)
			
			if 'V' in df.ix[idx]['Defect']:
				if self.first_species in df.ix[idx]['Defect']:
					color = 'blue'
				elif self.second_species in df.ix[idx]['Defect']:
					color = 'purple'
				elif self.third_species in df.ix[idx]['Defect']:
					color = 'black'
				elif self.fourth_species in df.ix[idx]['Defect']:
					color = 'gray'
			else:
				color = 'orange'
			
			if df.ix[idx]['q'] == 0:
				style = '-'
			else:
				style = '--'
			
			# If defect line is there, then update it; otherwise, draw a new defect line
			try:
				self.defect_plots[defect_name].set_data(self.Ef_pts, self.E_enthalpy)
			except:
				self.defect_plots[defect_name], = self.defects_diagram_plot_drawing.plot(self.Ef_pts, self.E_enthalpy, color = color, ls = style, label=str((df.ix[idx])['class']))
		
		# Remove labels before redrawing them at new position
		self.defects_diagram_plot_drawing.texts.clear()
		labelLines(self.defects_diagram_plot_drawing.get_lines(),align=False,xvals=[self.EVBM+0.15, self.EVBM+0.10, self.EVBM+0.15, self.EVBM+0.12, self.EVBM+0.15, self.EVBM+0.10, self.EVBM+0.15, self.EVBM+0.02, self.EVBM+0.17, self.EVBM+0.17, self.EVBM+0.17,self.EVBM+0.17, self.EVBM+0.17, self.EVBM+0.16, self.EVBM+0.14,self.EVBM+0.12, self.EVBM+0.15, self.EVBM+0.14, self.EVBM+0.15, self.EVBM+0.12, self.EVBM+0.10, self.EVBM+0.10],fontsize=10,bbox=dict(facecolor='white',alpha=0.8,edgecolor='white',pad=0.5))
		
		try:
			self.defects_stability_region.remove()
			self.defects_stability_region = self.defects_diagram_plot_drawing.fill_between(self.Ef_pts,0,-0.8,facecolor='#614126', interpolate=True, alpha=.1)
		except:
			self.defects_stability_region = self.defects_diagram_plot_drawing.fill_between(self.Ef_pts,0,-0.8,facecolor='#614126', interpolate=True, alpha=.1)
		
		self.defects_diagram_plot_canvas.draw()



