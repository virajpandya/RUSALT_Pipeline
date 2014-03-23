#!/usr/bin/python
'''
>>> Rutgers SALT Supernova Spectral Reduction Pipeline <<<

This module is the DRIVER for the pipeline. It imports all
of the necessary Python modules, PyRAF modules, PySALT modules,
and pipeline modules.

The driver assumes that the user is in their working directory
which contains the data to be reduced. It enters a loop
asking the user what she or he wishes to do, and calls another
module (or exits the program) based on the user's choice.

Please read the documentation to gain a thorough understanding
of how the pipeline works synergistically with PyRAF and PySALT
to reduce, extract, and calibrate long-slit spectral data.

*** Modifications ***
Sept. 26, 2013: Created module. -Viraj Pandya

'''

print "Welcome to the Rutgers Supernova Reduction Pipeline for SALT data!"
print "Importing the necessary modules..."

import sys # Standard python module used mainly to exit the pipeline.
import os # Standard python module used mainly for changing directories.
import shutil # Standard shutil module used mainly for copying files.
from glob import glob # For grabbing a list of filenames in working directory.
import numpy as np # NumPy is fundamental to many computations the pipeline does.
from scipy.interpolate import interp1d # For interpolation in 1-D spectra.
from scipy.optimize import leastsq # For fitting functions to 1-D spectra.
import lacosmicx # For removal of cosmic rays from 2-D science images.
import ds9 # To open 2-D images for the user's convenience.
import pyfits # To access and modify FITS data files.

# PyRAF is the main program used to reduce, extract, and calibrate data.
from pyraf import iraf
from iraf import images
from iraf import imutil
from iraf import immatch
from iraf import noao
from iraf import imred
from iraf import ccdred
from iraf import onedspec
from iraf import twodspec
from iraf import apextract
from iraf import longslit
from iraf import pysalt
from iraf import saltspec

# These are the pipeline modules which may be called by this driver or each other.
import dicts # Initializes the pipeline's global dictionaries which will be used across modules.
import tasks # The core of the pipeline: each function represents one step of the entire reduction process.
import params # Customizable parameters for the pipeline.
import pyrafCalls # Contains the functions which set the parameters for and call single PyRAF tasks.
import fullReduce # Automatically calls functions in the tasks module for a full reduction.
import indivReduce # Allows user to run individual reductions tasks.
import sortDicts # Sorts files in working directory into global dictionaries.
import pipeHistory # Informs user about pipeline processes run on each file.

# This section prompts the user for inputs asking them what they want to do.
while True:
	print "Please make sure you are in your working directory."
	print "Here are the options for reducing your data:"
	print "1. Full reduction of the spectral data."
	print "2. Run individual reduction tasks."
	print "3. Quit."
	while True:
		answer = raw_input("Please enter a number from the menu: ")
		if answer == '0' or answer == '1' or answer == '2' or answer == '3':
			break
		else:
			print "Invalid input. You must enter a number from the menu."
	if answer=='1':
		fullReduce.run()
	elif answer=='2':		
		indivReduce.run()
	elif answer=='3':
		sys.exit("Thanks for using this pipeline!")
	else:
		print 'You must pick an option from the menu.'
		continue 
