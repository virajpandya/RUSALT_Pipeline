
'''
>>> Rutgers SALT Supernova Spectral Reduction Pipeline <<<

This module prints the history of pipeline processes run on
each data file in the working directory. 

It also contains functions for modifying specific keywords 
in a FITS header. These keywords, notably 'RUPIPE' and
'RUIMGTYP', let the pipeline sort files into global dictionaries
automatically if the pipeline is run multiple times in a working
directory.

Please refer to the documentation for more information
about the purpose of maintaining a pipeline history.

*** Modifications ***
Sept. 30, 2013: Created module. -Viraj Pandya

'''

import sys # Standard python module used mainly to exit the pipeline.
import os # Standard os module used mainly for changing directories.
import shutil # Standard shutil module used mainly for copying files.
from glob import glob # For grabbing a list of filenames in working directory.
import numpy as np # NumPy is fundamental to many computations the pipeline does.
from scipy.interpolate import interp1d # For interpolation in 1-D spectra.
from scipy.optimize import leastsq # For fitting functions to 1-D spectra.
import lacosmicx # For removal of cosmic rays from 2-D science images.
import ds9 # To open 2-D images for the user's convenience.
import pyfits # To access and modify FITS data files.

# These are the pipeline modules which may be called by this driver or each other.
import dicts # Initializes the pipeline's global dictionaries which will be used across modules.
import params # Customizable parameters for the pipeline.
import pyrafCalls # Contains the functions which set the parameters for and call single PyRAF tasks.

def history():
	# inform user about what 'RUIMGTYP' classes exist.
	print "This task, history(), will print the 'RUIMGTYP' and 'RUPIPE' keywords present in each file's 0-header."
	print "RUIMGTYP can be: "
	print "science, arc, standard, flat, arc_sci, arc_std, bpm_sci, bpm_arc, bpm_std, bwm_sci, bwm_std."
	print "If RUIMGTYP is not found in the header, OBJECT will be printed in its place."
	print "RUPIPE will be shown as characters in a list each of which correspond to some particular pipeline process as defined below."
	print "The combination of RUIMGTYP and RUPIPE will tell you what processes have been run on that image."
	print "Here is what the characters in RUPIPE stand for; the order roughly corresponds to the full data reduction process."
	print "c: combined or count-scaled"
	print "n: normalized or flux-scaled"
	print "f: flat-fielded"
	print "o: original (bad pixel mask)"
	print "l: LAcosmicx-affiliated"
	print "a: specidentified arc (using PySALT)"
	print "r: specrectified (wavelength-calibrated using PySALT)"
	print "b: background-subtracted"
	print "e: extracted"
	print "j: apsum-extracted science for sky and sigma spectra"
	print "i: identified arc spectrum (PyRAF, not PySALT)"
	print "d: dispersion-corrected (using PyRAF's dispcor)"
	print "g: flux-calibrated"
	print "s: flux-normalized (using sarith)"
	print "y: sky spectrum"
	print "z: sigma spectrum"
	print "u: combined dispersion-corrected spectrum"
	print "v: combined flux-calibrated spectrum"
	print "w: combined sky spectrum"
	print "x: combined sigma spectrum"
	print "t: telluric-corrected spectrum"
	print "If RUPIPE is not found in the header, an empty list will be printed in its place."
	while True:
			choice = raw_input("Please enter 0 to see the dictionaries, or 1 to see the pipeline history: ")
			if choice == '0' or choice == '1':
				break
			else:
				print "Invalid input: you must enter 1 to proceed."
	if choice == '0':
		print "original flats:"
		print dicts.flats
		print "original arcs:"
		print dicts.arcs
		print "original standards:"
		print dicts.standards
		print "original sciences:"
		print dicts.sciences   
		print "combined flats:"
		print dicts.combflats
		print "normalized flats:"
		print dicts.normflats
		print "flat-fielded sciences:"
		print dicts.flatsciences
		print "flat-fielded arcs:"
		print dicts.flatarcs
		print "flat-fielded standards:"
		print dicts.flatstandards
		print "bad pixel masks for sciences:"
		print dicts.bpmsciences
		print "bad pixel masks for arcs:"
		print dicts.bpmarcs
		print "bad pixel masks for standards:"
		print dicts.bpmstandards
		print "lacosmicx-corrected sciences:"
		print dicts.laxsciences
		print "cosmic ray pixel masks for sciences:"
		print dicts.bpmlaxsciences
		print "wavelength solutions from PySALT:"
		print dicts.wavesols
		print "wavelength-calibrated arcs:"
		print dicts.wavearcs
		print "wavelength-corrected sciences:"
		print dicts.wavesciences
		print "wavelength-corrected standards:"
		print dicts.wavestandards
		print "background-subtracted sciences:"
		print dicts.backgroundsciences
		print "extracted sciences:"
		print dicts.extractedsciences
		print "extracted standards:"
		print dicts.extractedstandards
		print "extracted arcs for standards:"
		print dicts.extractedarcs_std
		print "extracted arcs for sciences:"
		print dicts.extractedarcs_sci
		print "apsum-extracted sciences:"
		print dicts.apsumsciences
		print "dispersion-corrected sciences:"
		print dicts.dispsciences
		print "dispersion-corrected standards:"
		print dicts.dispstandards
		print "bad wavelength masks for sciences:"
		print dicts.bwmsciences
		print "bad wavelength masks for standards:"
		print dicts.bwmstandards
		print "standard flux calibration files:"
		print dicts.stdfiles
		print "sensfunc flux calibration files:"
		print dicts.sensfiles
		print "flux-calibrated sciences:"
		print dicts.fluxsciences
		print "flux-calibrated standards:"
		print dicts.fluxstandards
		print "flux-scaled sciences:"
		print dicts.scaledfluxsciences
		print "flux-scaled standards:"
		print dicts.scaledstandards
		print "telluric-corrected sciences:"
		print dicts.tellsciences
		print "telluric-corrected standards:"
		print dicts.tellstandards
		print "sky (science) spectra:"
		print dicts.skysciences
		print "sigma (science) spectra:"
		print dicts.sigmasciences
		print "count-scaled sciences:"
		print dicts.scaleddispsciences
		print "combined science spectra:"
		print dicts.combinedspectra
	elif choice == '1':
		images = glob('*.fits')
		for img in images: # print img and 'RUPIPE' (in list form) and 'RUIMGTYP'
			hduimg = pyfits.open(img)
			hdrimg = hduimg[0].header.copy()
			imgtype = hdrimg.get('RUIMGTYP','')
			if imgtype == '':
				imgtype = hdrimg['OBJECT']
			pipestring = hdrimg.get('RUPIPE','')
			pipelist = list(pipestring)
			print img,imgtype,pipelist
			

# This function updates pipeline history keywords in a newly created/modified FITS file.
def updatePipeKeys(inputname,imagetype,procChar):
	hdukey = pyfits.open(inputname,mode='update')
	hdrkey = hdukey[0].header
	hdrkey['RUIMGTYP'] = imagetype
	if procChar != '': # '' is for adding 'science' to mbxgpP*.fits image for future auto-sorting
		try:
			proc = list(hdrkey['RUPIPE']) # add processed character flag if it's not already present
		except:
			proc = [] # in case image has never been processed so 'RUPIPE' keyword is not in the header
		if procChar not in proc:
			proc.append(procChar)
			proc = ''.join(proc)
			hdrkey['RUPIPE'] = proc
	hdukey.flush()
	hdukey.close()

