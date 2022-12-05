
from setuptools import setup, find_packages

source_packages = [	"vtandem", \
					"vtandem.dft", \
					"vtandem.visualization", \
					"vtandem.visualization.quaternary", \
					"vtandem.visualization.quaternary.quaternary_plots", \
					"vtandem.visualization.quaternary.quaternary_tabs", \
					"vtandem.visualization.ternary", \
					"vtandem.visualization.ternary.ternary_plots", \
					"vtandem.visualization.ternary.ternary_tabs", \
					"vtandem.visualization.binary", \
					"vtandem.visualization.binary.binary_plots", \
					"vtandem.visualization.binary.binary_tabs", \
					"vtandem.visualization.windows", \
					"vtandem.visualization.plots", \
					"vtandem.visualization.tabs", \
					"vtandem.visualization.utils"
					]
source_image_files = [ 	("logo", ("logo/LogoLong.png", "logo/LogoSmall.png")),
						("icon", ("icon/FolderBrowserIcon.png", "icon/QuestionIcon.png"))
						]

setup(
	name = "vtandem",
	version = "2021.11.18",
	description = "",
	author = "Michael Y. Toriyama, Jiaxing Qu, Lidia C. Gomes, Elif Ertekin",
	author_email = "mathtoriyama@gmail.com",
	url = "",
	packages = source_packages,
	data_files = source_image_files,
	py_modules = ["vtandem"],
	entry_points = {
		"console_scripts": [
			"vtandem = vtandem.vtandem:vtandem"
		]
	},
	install_requires=[
		"numpy>=1.16",
		"matplotlib>=3.0",
		"matplotlib-label-lines",
		"pymatgen",
		"periodictable",
		"PyQt5",
		"click>8.0",
	],
	dependency_links=['https://github.com/frssp/PyPolyhedron/master#egg=polyhedron-0.2.1'],
)




