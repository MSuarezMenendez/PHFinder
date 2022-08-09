from setuptools import setup

setup(
	name='PHFinder',
	version='1.0',
	description='Assisted detection of point heteroplasmy in Sanger sequencing chromatograms',
	author='Marcos Suarez',
	author_email='m.suarez.menendez@rug.nl',
	 install_requires=['biopython==1.73'], #external packages as dependencies
	scripts=['PHFinder/PHFinder.py', 'PHFinder/Bash_h.sh', 'PHFinder/Testing.sh']
)
