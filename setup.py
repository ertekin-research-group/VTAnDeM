
from setuptools import setup, find_packages

source_packages = [	"vtandem", \
					"vtandem.dft", \
					"vtandem.visualization", \
					"vtandem.visualization.quaternary", \
					"vtandem.visualization.quaternary.quaternary_scripts", \
					"vtandem.visualization.ternary", \
					"vtandem.visualization.ternary.ternary_scripts", \
					"vtandem.visualization.binary", \
					"vtandem.visualization.binary.binary_scripts", \
					"vtandem.visualization.widgets"
					]
source_image_files = [ 	("logo", ("logo/LogoLong.png", "logo/LogoSmall.png")),
						("icon", ("icon/FolderBrowserIcon.png", "icon/QuestionIcon.png"))
						]

setup(
	name = "vtandem",
	version = "2020.09.08",
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
)




