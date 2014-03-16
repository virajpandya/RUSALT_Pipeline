'''
>>> Rutgers SALT Supernova Spectral Reduction Pipeline <<<

This module is run whenever the data files in the working
directory need to be sorted into the global dictionaries.

The module takes advantage of each FITS data file's header
information, most notably the 'GR-ANGLE' (grating angle)
keyword. It also uses keywords added by the pipeline in
pipeline-processed files such as 'RUPIPE' and 'RUIMGTYP'.

Please refer to the documentation for more information
about how the sorting and accessing via dictionaries works.

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
import pipeHistory # Mainly used to access and modify pipeline processed header keywords in FITS images.

def sort():
	# This returns a list of all the fits filenames in the directory.
	images = glob('*.fits')
	# This creates a list of the "default" standard stars in the relevant pysalt subdirectory:
	# /usr/local/astro64/iraf/extern/pysalt/data/standards/spectroscopic/
	cwd = os.getcwd()
	os.chdir(params.standardsPath)
	possiblestandards = glob('m*.dat')
	os.chdir(cwd)
	# this gets rid of the 'm' in the beginning and '.dat' at the end so comparisons can be made with 'OBJECT'
	for n,s in enumerate(possiblestandards):
		snew = s.split('.')
		possiblestandards[n] = snew[0][1:].lower()
	# This loops through all the images and sorts them into one of the 4 dictionaries from above based on type and gr-angle.
	for img in images:	
		# if necessary, doing header error corrections
		hdulist = pyfits.open(img,mode='update')
		header = hdulist[0].header
		# Skip img if it is a SALTICAM (CCD) acquisition image.
		if header.get('INSTRUME','') == 'SALTICAM' or header.get('OBSMODE','') == 'IMAGING' or img[0:6] == 'mbxgpS':
			continue
		try: # checks if 'MASKTYP' keyword is set to 'LONGSLIT'
			if header['MASKTYP'] != 'LONGSLIT':
				header.update('MASKTYP','LONGSLIT')
				print img+": Changed 'MASKTYP' to 'LONGSLIT'"
		except:
			header.update('MASKTYP','LONGSLIT')
			print img+": Changed 'MASKTYP' to 'LONGSLIT'"
		try: # checks if 'GAIN' keyword exists (SALTMOSAIC may not have transferred it)
			if header.get('GAIN','') == '':
				if header.get('GAINSET','') == 'FAINT' and header.get('ROSPEED','') == 'SLOW':
					header.update('GAIN',2.146) # value from pysalt wiki documentation for this keyword-combo
					print img+": Changed 'GAIN' to 2.146"
		except:
			if header.get('GAINSET','') == 'FAINT' and header.get('ROSPEED','') == 'SLOW':
					header.update('GAIN',2.146) # value from pysalt wiki documentation for this keyword-combo
					print img+": Changed 'GAIN' to 2.146"
		header.update('AIRMASS',1.28) # standard AIRMASS for SALT (to prevent negative or high airmass for some images)
		print img+": Changed 'AIRMASS' to 1.28"
		hdulist.flush()
		hdulist.close()
		# sorting into dictionaries
		hdr = pyfits.getheader(img,0)
		# this truncates the gr-angle so it only has 2 decimal digits (string format)
		angle = hdr.get('GR-ANGLE','')  # since we are maintaining file structure with imcopy, most files should have GR-ANGLE preserved
		if angle != '':
			angle = str('%.3f'%angle)
			a = angle.split('.')
			a[1] = a[1][:2]
			angle = '.'.join(a)
		obj = hdr.get('OBJECT','') # bad pixel masks for example will not have 'OBJECT'
		try: # in case 'RUPIPE' or 'RUIMGTYP' keywords don't exist (image not pipeline-processed)
			processed = list(hdr['RUPIPE']) # returns character list instead of the string sequence
		except:
			processed = []
		try:
			imgclass = hdr['RUIMGTYP']
		except: # possible for mbxgpP raw science or std images with only 'RUIMGTYP' keyword set
			imgclass = ''
		# if non-pipeline-processed image, just check 'OBJECT'
		if processed == [] and imgclass == '':
			firstFour = img[0:4] # first four letters of img
			if firstFour == 'sens': # sensfunc file
				sensAngle = img[4:9]
				dicts.sensfiles[angle] = img
			elif obj=='FLAT':
				temp = dicts.flats.get(angle,[])
				temp.append(img)
				dicts.flats[angle] = temp
				print img+' sorted as flat with angle: '+angle
			elif obj=='ARC':
				dicts.arcs[angle] = img
				print img+' sorted as arc with angle: '+angle
			# check if obj is in list of possible SALT standard stars (lowercase comparisons to ignore case)
			elif obj.lower() in possiblestandards:
				temp = dicts.standards.get(angle,[])
				temp.append(img)
				dicts.standards[angle] = temp
				print img+' sorted as standard star with angle: '+angle
			# if obj is neither flat nor arc nor in possiblestandards, print obj and ask user to classify manually:
			# as standard (short name) or science (long, weird name maybe with 'RU' in it -- obvious)
			else:
				print 'Could not automatically classify '+img+' as a FLAT, ARC, or SALT standard star'
				print 'The object name is: '+obj
				print 'Please enter 0 if this name resembles an abbreviation for a standard star.'
				print 'Please enter 1 if this name is rather long and resembles a science image.'
				while True:
					answer = raw_input("Enter 0 or 1 (see instructions above): ")
					if answer == '0' or answer == '1':
						break
					else:
						print "Invalid input. You must enter either 0 or 1."
				if answer == '0':
					dicts.standards[angle] = img
					print img+' sorted as standard star with angle: '+angle
				elif answer == '1':
					pipeHistory.updatePipeKeys(inputname=img,imagetype='science',procChar='') # for future auto-sorting
					dicts.sciences[angle] = img
					print img+' sorted as science with angle: '+angle
		elif processed == [] and imgclass != '': # mbxgp processed before so 'RUIMGTYP' key present
			if imgclass == 'science':
				dicts.sciences[angle] = img
				print img+' sorted as science with angle: '+angle
			elif imgclass == 'standard':
				dicts.standards[angle] = img
				print img+' sorted as standard star with angle: '+angle
			elif imgclass == 'arc':
				dicts.arcs[angle] = img
				print img+' sorted as arc with angle: '+angle
			elif imgclass == 'flat':
				dicts.fkat[angle] = img
				print img+' sorted as flat with angle: '+angle
		else: # else, sort into global dictionaries based on processed and imgclass
			if 't' in processed and imgclass == 'science': # telluric-corrected science spectrum
				dicts.combinedspectra['tel'] = img
			elif 'k' in processed and imgclass == 'science': # cleaned combined science spectrum
				dicts.combinedspectra['clnsci'] = img
			elif 'k' in processed and imgclass == 'standard': # cleaned combined std star spectrum
				dicts.combinedspectra['clnstd'] = img
			elif 'x' in processed and imgclass == 'science': # combined sigma spectrum
				dicts.combinedspectra['sig'] = img
			elif 'w' in processed and imgclass == 'science': # combined sky spectrum
				dicts.combinedspectra['sky'] = img
			elif 'u' in processed and imgclass == 'science': # combined-disp science spectrum
				dicts.combinedspectra['dsp'] = img
			elif 'v' in processed and imgclass == 'standard': # combined sky spectrum
				dicts.combinedspectra['std'] = img
			elif 'v' in processed and imgclass == 'science': # combined-flux science spectrum
				dicts.combinedspectra['flx'] = img
			elif 'c' in processed and imgclass == 'science': # count-scaled science spectrum
				dicts.scaleddispsciences[angle] = img
			elif 'y' in processed and imgclass == 'science': # sky spectrum
				dicts.skysciences[angle] = img
			elif 'z' in processed and imgclass == 'science': # sigma spectrum
				dicts.sigmasciences[angle] = img
			elif 'n' in processed and imgclass == 'standard': # scaled standard star spectrum
				dicts.scaledstandards[angle] = img
			elif 'n' in processed and imgclass == 'science': # scaled science spectrum
				dicts.scaledfluxsciences[angle] = img
			elif 'g' in processed and imgclass == 'standard': # flux-calibrated standard star spectrum
				dicts.fluxstandards[angle] = img
				dicts.stdfiles[angle] = hdr.get('RUSTD','')
				dicts.sensfiles[angle] = hdr.get('RUSENS','')
			elif 'g' in processed and imgclass == 'science': # flux-calibrated science spectrum
				dicts.fluxsciences[angle] = img
				dicts.stdfiles[angle] = hdr.get('RUSTD','')
				dicts.sensfiles[angle] = hdr.get('RUSENS','')
			elif 'o' in processed and imgclass == 'bwm_sci': # bad wavelength mask for science
				dicts.bwmsciences[angle] = img
			elif 'o' in processed and imgclass == 'bwm_std': # bad wavelength mask for standard
				dicts.bwmstandards[angle] = img
			elif 'd' in processed and imgclass == 'science': # dispersion-corrected science
				dicts.dispsciences[angle] = img
			elif 'd' in processed and imgclass == 'standard': # dispersion-corrected standard
				dicts.dispstandards[angle] = img
			elif 'j' in processed and imgclass == 'science': # apsum-extracted science
				dicts.apsumsciences[angle] = img
			elif 'e' in processed and imgclass == 'science': # extracted science
				dicts.extractedsciences[angle] = img
			elif 'e' in processed and imgclass == 'standard': # extracted standard
				dicts.extractedstandards[angle] = img
			elif 'e' in processed and imgclass == 'arc_std': # extracted arc_std
				dicts.extractedarcs_std[angle] = img
			elif 'e' in processed and imgclass == 'arc_sci': # extracted arc_sci
				dicts.extractedarcs_sci[angle] = img
			elif 'b' in processed and imgclass == 'science': # background-subtracted science
				dicts.backgroundsciences[angle] = img
			elif 'r' in processed and imgclass == 'arc': # rectified arc
				dicts.wavearcs[angle] = img
			elif 'r' in processed and imgclass == 'science': # rectified science
				dicts.wavesciences[angle] = img
			elif 'r' in processed and imgclass == 'standard': # rectified standard
				dicts.wavestandards[angle] = img
			elif 'a' in processed and imgclass == 'arc': # flatarc that had specidentify run on it
				dicts.wavesols[angle] = img
				# since 'a' will always be with 'f', add img to flatarcs[angle] as well
				dicts.flatarcs[angle] = img
			elif 'l' in processed and imgclass == 'science': # lacosmicx-corrected sciences
				dicts.laxsciences[angle] = img
			elif 'l' in processed and imgclass == 'bpm_sci': # cosmic (ray) pixel mask from lacosmicx (for science)
				dicts.bpmlaxsciences[angle] = img
			elif 'o' in processed and imgclass == 'bpm_sci': # bad pixel mask for science
				dicts.bpmsciences[angle] = img
			elif 'o' in processed and imgclass == 'bpm_arc': # bad pixel mask for arc
				dicts.bpmarcs[angle] = img
			elif 'o' in processed and imgclass == 'bpm_std': # bad pixel mask for standard
				dicts.bpmstandards[angle] = img		
			elif 'f' in processed and imgclass == 'science': # flat-fielded science
				dicts.flatsciences[angle] = img
			elif 'f' in processed and imgclass == 'arc': # flat-fielded arc
				# this will actually never run because 'a' will always be in RUPIPE with 'f'
				# flat-fielded arc is added to flatarcs[angle] when the 'a' processing (wavesols) is done
				dicts.flatarcs[angle] = img
			elif 'f' in processed and imgclass == 'standard': # flat-fielded standard
				dicts.flatstandards[angle] = img	
			elif 'n' in processed and imgclass == 'flat': # normalized flat
				dicts.normflats[angle] = img
			elif 'c' in processed and imgclass == 'flat': # combined flat
				dicts.combflats[angle] = img
			
	# printing initial global dictionaries
	print 'The original flat image filenames are:'
	print dicts.flats
	print 'The original science image filenames are:'
	print dicts.sciences
	print 'The original standard star image filenames are:'
	print dicts.standards
	print 'The original arc image filenames are:'
	print dicts.arcs
