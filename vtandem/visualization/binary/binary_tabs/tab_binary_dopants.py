
__author__ = 'Michael_Lidia_Jiaxing_Elif'
__name__ = 'VTAnDeM_Visualization-Toolkit-for-Analyzing-Defects-in-Materials'


from vtandem.visualization.tabs.tab_dopants import Tab_Dopants

class Tab_Binary_Dopants(Tab_Dopants):
	
	def __init__(self, 	main_compound = None, \
						first_element = None, \
						second_element = None, \
						compounds_info = None, \
						defects_data = None, \
						main_compound_info = None, \
						dos_data = None, \
						show_defects_diagram = True, \
						show_carrier_concentration = True
						):	# User specifies the main compound and its constituents
		
		# Label the first, second, and third species of the atoms in the ternary compound
		self.first_element = first_element
		self.second_element = second_element
		self.elements_list = [self.first_element, self.second_element]	# Species list (order MAY change)

		super().__init__(		type = "binary", \
								main_compound = main_compound, \
								compounds_info = compounds_info, \
								defects_data = defects_data, \
								main_compound_info = main_compound_info, \
								dos_data = dos_data, \
								show_defects_diagram = show_defects_diagram, \
								show_carrier_concentration = show_carrier_concentration)
