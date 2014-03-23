
'''
>>> Rutgers SALT Supernova Spectral Reduction Pipeline <<<

This module allows the user to run individual tasks, i.e., to 
call individual functions from the tasks.py module.

Note that this module will be rewritten to allow the user to
specify custom input images (dictionary arguments).

Please refer to the documentation for more information
about what each individual step in the reduction process does.

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
import params # Customizable parameters for the pipeline.
import pyrafCalls # Contains the functions which set the parameters for and call single PyRAF tasks.
import tasks # Core pipeline module.
import pipeHistory # Pipeline module for tracking processes run on files.
import sortDicts # Pipeline's module for sorting files into global dictionaries.

def run():
	# let the user pick a single task, and once that's done, loop to let them pick another task or return to the main menu
	while True:
		if params.firstIndivReduce == True:
			sortDicts.sort() # sort files in working directory into global dictionaries
			params.firstIndivReduce = False		
		print "Please choose one of the reduction tasks below."
		print "0. Show global dictionaries, or history of pipeline processes for each image."
		print "1. Combine flat images."
		print "2. Normalize flat images."
		print "3. Flat-field your science images."
		print "4. Flat-field your arc images."
		print "5. Flat-field your standard star images."
		print "6. Create a bad pixel mask for each science, arc, and standard star image."
		print "7. Run LAcosmicx on flat-fielded science images."
		print "8. Find a preliminary wavelength solution using PySALT."
		print "9. Rectify your science images using the PySALT wavelength solution."
		print "10. Rectify your standard star images using the PySALT wavelength solution."
		print "11. Background-subtract your science images."
		print "12. Extract science spectra."
		print "13. Extract standard star spectra."
		print "14. Extract arc spectra."
		print "15. Extract 'true' sky and sigma spectra from non-background-subtracted science images."
		print "16. Identify (and Reidentify) lines in your arc spectra for a second wavelength solution." # runs reidentify function too
		print "17. Dispersion-correct your science spectra using the second wavelength solution."
		print "18. Dispersion-correct your standard star spectra using the second wavelength solution."
		print "19. Create bad wavelength masks for dispersion-corrected science (and perhaps standard star) spectra."
		print "20. Flux-calibrate your science and standard star spectra."
		print "21. Sigma-clip science and standard star spectra."
		print "22. Scale your flux-calibrated science spectra to similar flux values at overlapping wavelengths."
		print "23. Scale your standard star spectra to similar flux values at overlapping wavelengths."
		print "24. Copy sky and sigma spectra from flux-calibrated spectra to different spectral files."
		print "25. Scale your dispersion-corrected science spectra to similar count values at overlapping wavelengths."
		print "26. Combine your dispersion-corrected and flux-calibrated science and standard star spectra."
		print "27. Clean combined science and standard star spectra of leftover chip gaps and highly deviant pixels."
		print "28. Telluric-correct combined science spectrum."
		print "29. Return to the main menu."
		options = range(0,31) # list containing integers from 0 to 30
		for i,v in enumerate(options): # convert all integers to string for raw_input comparison
			options[i] = str(v)
		while True:
			choice = raw_input("Please enter a number from the submenu: ")
			if choice in options:
				break
			else:
				print "Invalid input: that choice is not in the submenu."
		if choice == '0': # show 'RUIMGTYP' and 'RUPIPE' for each file that has them
			print "You have chosen to see the history of pipeline processes run on each file."
			pipeHistory.history()
		elif choice=='1': # combine flats
			print "You have chosen to combine flat-field images."
			print "The flat-field images indexed by GR-ANGLE are: "
			print dicts.flats
			if dicts.flats == {}:
				print "There are no flats to combine:"
				continue
			if dicts.combflats != {}:
				print "There are already combined flats: "
				print dicts.combflats
				print "If you continue, the current combined flats will be replaced with new ones."
				print "An alternative is to work in a cloned directory (using Option 0 from the main menu)."
				while True: 
					moveforward = raw_input("Enter 0 to return to the submenu, or 1 to move forward: ")
					if moveforward == '0' or moveforward == '1':
						break
					else:
						print "Invalid input: you must enter either 0 or 1."
				if moveforward == '0':
					print "You chose to return to the submenu instead of replacing the old combined flats."
					continue
			imdelseq = ''
			for combflat in dicts.combflats:
				imdelseq = imdelseq+','+combflat
			pyrafCalls.run_imdel(inputname=imdelseq)
			tasks.combineflats(allflats=dicts.flats,indiv=True)
		elif choice=='2': # normalize flats
			print "You have chosen to normalize combined flat-field images."
			print "The combined flat-field images are:"
			print dicts.combflats
			if dicts.combflats == {}:
				print "There are no combined flats to normalize."
				continue
			if dicts.normflats != {}:
				print "There are already normalized flats: "
				print dicts.normflats
				print "If you continue, the current normalized flats will be replaced with new ones."
				print "An alternative is to work in a cloned directory (using Option 0 from the main menu)."
				while True: 
					moveforward = raw_input("Enter 0 to return to the submenu, or 1 to move forward: ")
					if moveforward == '0' or moveforward == '1':
						break
					else:
						print "Invalid input: you must enter either 0 or 1."
				if moveforward == '0':
					print "You chose to return to the submenu instead of replacing the old normalized flats."
					continue
			imdelseq = ''
			for normflat in dicts.normflats:
				imdelseq = imdelseq+','+normflat
			pyrafCalls.run_imdel(inputname=imdelseq)
			tasks.normalizeflats(combinedflats=dicts.combflats,indiv=True)
		elif choice=='3': # flat-field science images
			print "You have chosen to flat-field science images."
			print "The science images are:"
			print dicts.sciences
			print "The normalized flats are:"
			print dicts.normflats
			if dicts.sciences == {}:
				print "There are no science images to flat-field."
				continue
			if dicts.normflats == {}:
				print "There are no normalized flats to divide by."
				continue
			if dicts.flatsciences != {}:
				print "There are already flat-fielded science images: "
				print dicts.flatsciences
				print "If you continue, the current flat-fielded science images will be replaced with new ones."
				print "An alternative is to work in a cloned directory (using Option 0 from the main menu)."
				while True: 
					moveforward = raw_input("Enter 0 to return to the submenu, or 1 to move forward: ")
					if moveforward == '0' or moveforward == '1':
						break
					else:
						print "Invalid input: you must enter either 0 or 1."
				if moveforward == '0':
					print "You chose to return to the submenu instead of replacing the old flat-fielded science images."
					continue
			imdelseq = ''
			for flatscience in dicts.flatsciences:
				imdelseq = imdelseq+','+dicts.flatsciences[flatscience]
			pyrafCalls.run_imdel(inputname=imdelseq)
			tasks.flattensciences(scienceimages=dicts.sciences,normalizedflats=dicts.normflats,indiv=True)
		elif choice=='4': # flat-field arc images
			print "You have chosen to flat-field arc images."
			print "The original arc images are:"
			print dicts.arcs
			print "The normalized flats are:"
			print dicts.normflats
			if dicts.arcs == {}:
				print "There are no arc images to flat-field."
				continue
			if dicts.normflats == {}:
				print "There are no normalized flats to divide by."
				continue
			if dicts.flatarcs != {}:
				print "There are already flat-fielded arc images: "
				print dicts.flatarcs
				print "If you continue, the current flat-fielded arc images will be replaced with new ones."
				print "An alternative is to work in a cloned directory (using Option 0 from the main menu)."
				while True: 
					moveforward = raw_input("Enter 0 to return to the submenu, or 1 to move forward: ")
					if moveforward == '0' or moveforward == '1':
						break
					else:
						print "Invalid input: you must enter either 0 or 1."
				if moveforward == '0':
					print "You chose to return to the submenu instead of replacing the old flat-fielded arc images."
					continue
			imdelseq = ''
			for flatarc in dicts.flatarcs:
				imdelseq = imdelseq+','+dicts.flatarcs[flatarc]
			pyrafCalls.run_imdel(inputname=imdelseq)
			tasks.flattenarcs(arcimages=dicts.arcs,normalizedflats=dicts.normflats,indiv=True)
		elif choice=='5': # flat-field standard star images
			print "You have chosen to flat-field standard star images."
			print "The original standard star images are: "
			print dicts.standards
			print "The normalized flats are:"
			print dicts.normflats
			if dicts.standards == {}:
				print "There are no standard star images to flat-field."
				continue
			if dicts.normflats == {}:
				print "There are no normalized flats to divide by."
				continue
			if dicts.flatstandards != {}:
				print "There are already flat-fielded standard star images: "
				print dicts.flatstandards
				print "If you continue, the current flat-fielded standard star images will be replaced with new ones."
				print "An alternative is to work in a cloned directory (using Option 0 from the main menu)."
				while True: 
					moveforward = raw_input("Enter 0 to return to the submenu, or 1 to move forward: ")
					if moveforward == '0' or moveforward == '1':
						break
					else:
						print "Invalid input: you must enter either 0 or 1."
				if moveforward == '0':
					print "You chose to return to the submenu instead of replacing the old flat-fielded standard star images."
					continue
			imdelseq = ''
			for flatstandard in dicts.flatstandards:
				imdelseq = imdelseq+','+dicts.flatstandards[flatstandard]
			pyrafCalls.run_imdel(inputname=imdelseq)
			tasks.flattenstandards(standardimages=dicts.standards,normalizedflats=dicts.normflats,indiv=True)
		elif choice=='6': # create bad pixels for all flat-fielded images
			print "You have chosen to create bad pixel masks for all (flat-fielded) science, arc, and standard star images."
			print "The flat-fielded science images are:"
			print dicts.flatsciences
			print "The flat-fielded arc images are:"
			print dicts.flatarcs
			print "The flat-fielded standard star images are:"
			print dicts.flatstandards
			if dicts.flatsciences == {}:
				print "There are no flat-fielded science images to create bad pixel masks for."
				print "Cannot proceed."
				continue
			if dicts.flatarcs == {}:
				print "There are no flat-fielded arc images to create bad pixel masks for."
				print "Will not create bad pixel masks for arc images." 
			if dicts.flatstandards == {}:
				print "There are no flat-fielded standard star images."
				print "Will not create bad pixel masks for standard star images." # no flatstandards => no standards
			if dicts.bpmsciences != {} or dicts.bpmarcs != {} or dicts.bpmstandards != {}:
				print "There are already bad pixel masks for science, arc, or standard star images: "
				print dicts.bpmsciences
				print dicts.bpmarcs
				print dicts.bpmstandards
				print "If you continue, the current bad pixel masks will be replaced with new ones."
				print "An alternative is to work in a cloned directory (using Option 0 from the main menu)."
				while True: 
					moveforward = raw_input("Enter 0 to return to the submenu, or 1 to move forward: ")
					if moveforward == '0' or moveforward == '1':
						break
					else:
						print "Invalid input: you must enter either 0 or 1."
				if moveforward == '0':
					print "You chose to return to the submenu instead of replacing the old bad pixel masks."
					continue
			imdelseq = ''
			for bpmscience in dicts.bpmsciences:
				imdelseq = imdelseq+','+bpmscience
			for bpmarc in dicts.bpmarcs:
				imdelseq = imdelseq+','+bpmarc
			for bpmstandard in dicts.bpmstandards:
				imdelseq = imdelseq+','+bpmstandard
			pyrafCalls.run_imdel(inputname=imdelseq)
			tasks.maskimages(scienceimages=dicts.flatsciences,arcimages=dicts.flatarcs,standardimages=dicts.flatstandards,indiv=True)
		elif choice=='7': # run lacosmicx on flat-fielded science images
			print "You have chosen to run LAcosmicx on (flat-fielded) science images."
			print "The flat-fielded science images are:"
			print dicts.flatsciences
			print "The bad pixel masks for the flat-fielded sciences are:"
			print dicts.bpmsciences
			if dicts.flatsciences == {}:
				print "There are no flat-fielded science images to run LAcosmicx on."
				print "Cannot proceed."
				continue
			if dicts.bpmsciences == {}:
				print "There are no bad pixel masks for the chip gaps in the flat-fielded science images."
				print "Cannot proceed."
				continue
			if dicts.laxsciences != {} or dicts.bpmlaxsciences != {}:
				print "There are already LAcosmicx-corrected science images or LAcosmicx pixel masks: "
				print dicts.flatsciences
				print dicts.bpmlaxsciences
				print "If you continue, the current LAcosmicx-corrected science images and pixel masks will be replaced with new ones."
				print "An alternative is to work in a cloned directory (using Option 0 from the main menu)."
				while True: 
					moveforward = raw_input("Enter 0 to return to the submenu, or 1 to move forward: ")
					if moveforward == '0' or moveforward == '1':
						break
					else:
						print "Invalid input: you must enter either 0 or 1."
				if moveforward == '0':
					print "You chose to return to the submenu instead of replacing the old LAcosmicx-corrected science images and pixel masks."
					continue
			imdelseq = ''
			for laxscience in dicts.laxsciences:
				imdelseq = imdelseq+','+laxscience
			for bpmlaxscience in dicts.bpmlaxsciences:
				imdelseq = imdelseq+','+bpmlaxscience
			pyrafCalls.run_imdel(inputname=imdelseq)
			tasks.LAxsciences(scienceimages=dicts.flatsciences,bpmasks=dicts.bpmsciences,indiv=True)
		elif choice=='8': # use specidentify for a preliminary wavelength solution, only for flat-arcs for now
			print "You have chosen to find a preliminary wavelength solution with flat-fielded arcs using PySALT."
			print "Here are the global dictionaries of flatarcs, solutions, and rectified arcs indexed by GR-ANGLE:"
			print dicts.flatarcs
			print dicts.wavesols
			print dicts.wavearcs
			while True:
				gratingangle = raw_input("Please enter the GR-ANGLE keyword of the arc you wish to run specidentify on: ")
				angles = dicts.flatarcs.keys()
				if gratingangle not in angles: # angles[i] is already a string, not a double/float
					print "Invalid input: that number is not one of the GR-ANGLE keywords for the dictionary of flat arcs."
				else:
					break
			anglewavesol = dicts.wavesols.get(gratingangle,'') # check if wavelength solution (.db file) for this angle already exists
			if anglewavesol != '':
				print "There is already a preliminary wavelength solution for angle "+gratingangle+": "
				print anglewavesol
				print "If you continue, the current preliminary wavelength solution for this angle will be replaced with a new one."
				print "An alternative is to work in a cloned directory (using Option 0 from the main menu)."
				while True: 
					moveforward = raw_input("Enter 0 to return to the submenu, or 1 to move forward: ")
					if moveforward == '0' or moveforward == '1':
						break
					else:
						print "Invalid input: you must enter either 0 or 1."
				if moveforward == '0':
					print "You chose to return to the submenu instead of replacing the old preliminary wavelength solution."
					continue
				os.remove(anglewavesol) # deleting anglewavesol so there are no conflicts in specidentifyarcs()
				wavearc = dicts.wavearcs.get(gratingangle,'')
				if wavearc != '':
					while True: 
						moveforward = raw_input("Rectified arc found; enter 0 to return to submenu, or 1 to delete: ")
						if moveforward == '0' or moveforward == '1':
							break
						else:
							print "Invalid input: you must enter either 0 or 1."
					if moveforward == '0':
						print "You chose to return to the submenu instead of replacing the old preliminary wavelength solution."
						continue
					os.remove(wavearc)
			angleflatarc = dicts.flatarcs[gratingangle] # flat arc with gratingangle
			inputflatarcs = {gratingangle:angleflatarc} # dictionary consisting of only one GR-ANGLE keyword
			tasks.specidentifyarcs(arcimages=inputflatarcs,indiv=True)
		elif choice=='9': # use specrectify on laxsciences with their associated specidentify solution
			print "You have chosen to rectify LAcosmicx-corrected science images with their PySALT wavelength solution."
			if dicts.wavesols == {}:
				print "The dictionary of wavelength solutions obtained via PySALT is completely empty: "
				print dicts.wavesols
				print "You must find preliminary wavelength solutions before rectifying images."
				print "Use Option 6 in the submenu to use PySALT to find wavelength solutions."
				continue
			print "Here are the global dictionaries of LAcosmicx-corrected science images indexed by GR-ANGLE:"
			print dicts.laxsciences
			while True:
				gratingangle = raw_input("Please enter the GR-ANGLE keyword of the science image you wish to run specrectify on: ")
				angles = dicts.laxsciences.keys()
				if gratingangle not in angles: # angles[i] is already a string, not a double/float
					print "Invalid input: that number is not one of the GR-ANGLE keywords for the dictionary of LAcosmicx-corrected sciences."
				else:
					break
			anglewavesol = dicts.wavesols.get(gratingangle,'')
			if anglewavesol == '':  # check if wavelength solution (.db file) for this angle exists
				print "There is no PySALT wavelength solution for this LAcosmicx-corrected science's angle: "
				print wavesols
				print "Cannot rectify image without a PySALT wavelength solution. Run Option 6 first."
			else: # check if the corresponding rectified image already exists
				anglewavescience = dicts.wavesciences.get(gratingangle,'')
				if anglewavescience != '':
					print "There is already a rectified science image for angle "+gratingangle+": "
					print anglewavescience
					print "If you continue, the current rectified science image for this angle will be replaced with a new one."
					print "An alternative is to work in a cloned directory (using Option 0 from the main menu)."
					while True: 
						moveforward = raw_input("Enter 0 to return to the submenu, or 1 to move forward: ")
						if moveforward == '0' or moveforward == '1':
							break
						else:
							print "Invalid input: you must enter either 0 or 1."
					if moveforward == '0':
						print "You chose to return to the submenu instead of replacing the old rectified science image."
						continue
					os.remove(anglewavescience) # deleting anglewavescience so there are no conflicts
				anglelaxscience = dicts.laxsciences[gratingangle] # lacosmicx-corrected science with gratingangle
				inputlaxsciences = {gratingangle:anglelaxscience} # dictionary consisting of only one GR-ANGLE keyword
				inputwavesols = {gratingangle:anglewavesol}
				tasks.wavecalsci(scienceimages=inputlaxsciences,sols=inputwavesols,indiv=True)
		elif choice=='10': # use specrectify on flatstandards with their associated specidentify solution
			print "You have chosen to rectify flat-fielded standard star images with their PySALT wavelength solution."
			if dicts.wavesols == {}:
				print "The dictionary of wavelength solutions obtained via PySALT is completely empty: "
				print dicts.wavesols
				print "You must find preliminary wavelength solutions before rectifying images."
				print "Use Option 6 in the submenu to use PySALT to find wavelength solutions."
				continue
			print "Here are the global dictionaries of flat-fielded science images indexed by GR-ANGLE:"
			print dicts.flatstandards
			while True:
				gratingangle = raw_input("Please enter the GR-ANGLE keyword of the standard star image you wish to run specrectify on: ")
				angles = dicts.flatstandards.keys()
				if gratingangle not in angles: # angles[i] is already a string, not a double/float
					print "Invalid input: that number is not one of the GR-ANGLE keywords for the dictionary of flat standards."
				else:
					break
			anglewavesol = dicts.wavesols.get(gratingangle,'')
			if anglewavesol == '':  # check if wavelength solution (.db file) for this angle exists
				print "There is no PySALT wavelength solution for this flat standard's angle: "
				print dicts.wavesols
				print "Cannot rectify image without a PySALT wavelength solution. Run Option 6 first."
			else: # check if the corresponding rectified image already exists
				anglewavestandard = dicts.wavestandards.get(gratingangle,'')
				if anglewavestandard != '':
					print "There is already a rectified standard star image for angle "+gratingangle+": "
					print anglewavestandard
					print "If you continue, the current rectified standard star image for this angle will be replaced with a new one."
					print "An alternative is to work in a cloned directory (using Option 0 from the main menu)."
					while True: 
						moveforward = raw_input("Enter 0 to return to the submenu, or 1 to move forward: ")
						if moveforward == '0' or moveforward == '1':
							break
						else:
							print "Invalid input: you must enter either 0 or 1."
					if moveforward == '0':
						print "You chose to return to the submenu instead of replacing the old rectified standard star image."
						continue
					os.remove(anglewavestandard) # deleting anglewavestandard so there are no conflicts
				angleflatstandard = dicts.flatstandards[gratingangle] # flat standard with gratingangle
				inputflatstandards = {gratingangle:angleflatstandard} # dictionary consisting of only one GR-ANGLE keyword
				inputwavesols = {gratingangle:anglewavesol}
				tasks.wavecalstd(standardimages=inputflatstandards,sols=inputwavesols,indiv=True)
		elif choice=='11': # background-subtract sciences
			print "You have chosen to background-subtract rectified science images."
			print "The rectified science images are:"
			print dicts.wavesciences
			if dicts.wavesciences == {}:
				print "There are no rectified science images to background-subtract."
				continue
			if dicts.backgroundsciences != {}: # check if the corresponding background-subtracted image already exists
				print "There are already background-subtracted science images: "
				print dicts.backgroundsciences
				print "If you continue, the current background-subtracted science images will be replaced with new ones."
				print "An alternative is to work in a cloned directory (using Option 0 from the main menu)."
				while True: 
					moveforward = raw_input("Enter 0 to return to the submenu, or 1 to move forward: ")
					if moveforward == '0' or moveforward == '1':
						break
					else:
						print "Invalid input: you must enter either 0 or 1."
				if moveforward == '0':
					print "You chose to return to the submenu instead of replacing the old background-subtracted science images."
					continue
			imdelseq = ''
			for wavescience in dicts.wavesciences:
				imdelseq = imdelseq+','+wavescience
			pyrafCalls.run_imdel(inputname=imdelseq)
			tasks.subtractbackground(scienceimages=dicts.wavesciences,indiv=True)
		elif choice=='12': # extract sciences (only from background-subtracted for now)
			print "You have chosen to extract background-subtracted science spectra."
			print "Here is the global dictionary of background-subtracted science images indexed by GR-ANGLE:"
			print dicts.backgroundsciences
			while True:
				gratingangle = raw_input("Please enter the GR-ANGLE keyword of the background-subtracted science image you wish to extract: ")
				angles = dicts.backgroundsciences.keys()
				if gratingangle not in angles: # angles[i] is already a string, not a double/float
					print "Invalid input: that number is not one of the GR-ANGLE keywords for the dictionary of background-subtracted sciences."
				else:
					break
			angleextscience = dicts.extractedsciences.get(gratingangle,'')
			if angleextscience != '': # check if the corresponding extracted image already exists
				print "There is already an extracted science spectrum for angle "+gratingangle+": "
				print angleextscience
				print "If you continue, the current extracted science spectrum for this angle will be replaced with a new one."
				print "An alternative is to work in a cloned directory (using Option 0 from the main menu)."
				while True: 
					moveforward = raw_input("Enter 0 to return to the submenu, or 1 to move forward: ")
					if moveforward == '0' or moveforward == '1':
						break
					else:
						print "Invalid input: you must enter either 0 or 1."
				if moveforward == '0':
					print "You chose to return to the submenu instead of replacing the old extracted science spectrum."
					continue
				os.remove(angleextscience) # deleting angleextscience so there are no conflicts
			anglebkgscience = dicts.backgroundsciences[gratingangle] 
			inputbackgroundsciences = {gratingangle:anglebkgscience}
			tasks.extractsciences(scienceimages=inputbackgroundsciences,indiv=True)
		elif choice=='13': # extract standard stars (only from rectified for now)
			print "You have chosen to extract rectified standard star spectra."
			print "Here is the global dictionary of rectified standard star images indexed by GR-ANGLE:"
			print dicts.wavestandards
			while True:
				gratingangle = raw_input("Please enter the GR-ANGLE keyword of the rectified standard star image you wish to extract: ")
				angles = dicts.wavestandards.keys()
				if gratingangle not in angles: # angles[i] is already a string, not a double/float
					print "Invalid input: that number is not one of the GR-ANGLE keywords for the dictionary of rectified standards."
				else:
					break
			angleextstandard = dicts.extractedstandards.get(gratingangle,'')
			if angleextstandard != '': # check if the corresponding extracted image already exists
				print "There is already an extracted standard star spectrum for angle "+gratingangle+": "
				print angleextstandard
				print "If you continue, the current extracted standard star spectrum for this angle will be replaced with a new one."
				print "An alternative is to work in a cloned directory (using Option 0 from the main menu)."
				while True: 
					moveforward = raw_input("Enter 0 to return to the submenu, or 1 to move forward: ")
					if moveforward == '0' or moveforward == '1':
						break
					else:
						print "Invalid input: you must enter either 0 or 1."
				if moveforward == '0':
					print "You chose to return to the submenu instead of replacing the old extracted standard star spectrum."
					continue
				os.remove(angleextstandard) # deleting angleextstandard so there are no conflicts
			anglewavestandard = dicts.wavestandards[gratingangle] 
			inputwavestandards = {gratingangle:anglewavestandard}
			tasks.extractstandards(standardimages=inputwavestandards,indiv=True)
		elif choice=='14': # extract arcs (only from rectified for now) -- extracts twice (for sci and for std)
			print "You have chosen to extract rectified arc spectra."
			print "Here is the global dictionary of rectified arc images indexed by GR-ANGLE:"
			print dicts.wavearcs
			while True:
				gratingangle = raw_input("Please enter the GR-ANGLE keyword of the rectified arc image you wish to extract: ")
				angles = dicts.wavearcs.keys()
				if gratingangle not in angles: # angles[i] is already a string, not a double/float
					print "Invalid input: that number is not one of the GR-ANGLE keywords for the dictionary of rectified arcs."
				else:
					break
			anglebkgscience = dicts.backgroundsciences.get(gratingangle,'') # existence of reference science 2D image
			anglewavestandard = dicts.wavestandards.get(gratingangle,'') # existence of reference standard 2D image
			if anglebkgscience == '' or anglewavestandard == '': 
				print "A reference (science or standard star) image does not exist for apsum: "
				print anglebkgscience
				print anglewavestandard
				print "Cannot proceed."
				continue
			angleextarc_sci = dicts.extractedarcs_sci.get(gratingangle,'') #sci-extracted arc
			angleextarc_std = dicts.extractedarcs_std.get(gratingangle,'') #std-extracted arc
			if angleextarc_sci != '' or angleextarc_std != '': # check if one or the other (sci or std) corresponding extracted arcs already exist
				print "There is already an extracted arc spectrum for angle "+gratingangle+": "
				print angleextarc_sci
				print angleextarc_std
				print "If you continue, both of the current extracted arc spectra for this angle will be replaced with new ones."
				print "An alternative is to work in a cloned directory (using Option 0 from the main menu)."
				while True: 
					moveforward = raw_input("Enter 0 to return to the submenu, or 1 to move forward: ")
					if moveforward == '0' or moveforward == '1':
						break
					else:
						print "Invalid input: you must enter either 0 or 1."
				if moveforward == '0':
					print "You chose to return to the submenu instead of replacing the old extracted arc spectra."
					continue
				os.remove(angleextarc_sci) # deleting both sci- and std-extracted arc spectra since they're both done by extractarcs() anyway
				os.remove(angleextarc_std)
			anglewavearc = dicts.wavearcs[gratingangle]
			inputwavearcs = {gratingangle:anglewavearc}
			inputrefsciences = {gratingangle:anglebkgscience}
			inputrefstandards = {gratingangle:anglewavestandard}
			tasks.extractarcs(arcimages=inputwavearcs,refsciences=inputrefsciences,refstandards=inputrefstandards,indiv=True)
		elif choice=='15': # extract non-background subtracted science spectra
			print "You have chosen to extract the true science sky and sigma spectra."
			print "The non-background-subtracted science images will be used for this task."
			print "Here is the global dictionary of rectified science images indexed by GR-ANGLE:"
			print dicts.wavesciences
			print "Here is the global dictionary of background-subtracted science 2D IMAGES indexed by GR-ANGLE:"
			print dicts.backgroundsciences
			print "Here is the global dictionary of background-subtracted science SPECTRA indexed by GR-ANGLE:"
			print dicts.extractedsciences
			if dicts.wavesciences == {} or dicts.backgroundsciences == {} or dicts.extractedsciences == {}:
				print "Cannot proceed: one of the required dictionaries printed above is empty."
				continue
			if dicts.apsumsciences != {}: # check if apsum-extracted sciences already exist
				print "There are already apsum-extracted science spectra."
				print dicts.apsumsciences
				print "If you continue, any apsum-extracted science spectra will be replaced with new ones."
				print "An alternative is to work in a cloned directory (using Option 0 from the main menu)."
				while True: 
					moveforward = raw_input("Enter 0 to return to the submenu, or 1 to move forward: ")
					if moveforward == '0' or moveforward == '1':
						break
					else:
						print "Invalid input: you must enter either 0 or 1."
				if moveforward == '0':
					print "You chose to return to the submenu instead of replacing the old apsum-extracted science spectra."
					continue
				for k in dicts.apsumsciences.keys(): # delete each old apsum-extracted science spectrum
					os.remove(dicts.apsumsciences[k])
			tasks.extractSciSigma(scienceimages=dicts.wavesciences,refsciences=dicts.backgroundsciences,mainsciences=dicts.extractedsciences,indiv=True)
		elif choice=='16': # identify (and reidentify) lines in extracted arcs (both sci and std)
			print "You have chosen to identify lines in the arc spectra for a second wavelength solution."
			print "Here are the global dictionaries of arc spectra indexed by GR-ANGLE:"
			print dicts.extractedarcs_sci
			print dicts.extractedarcs_std
			print "This task runs identify on the science-arcs interactively, and reidentify on the standard-star-arcs non-interactively."
			while True:
				gratingangle = raw_input("Please enter the GR-ANGLE keyword of the arcs (science) you wish to run identify on: ")
				angles = dicts.extractedarcs_sci.keys()
				if gratingangle not in angles: # angles[i] is already a string, not a double/float
					print "Invalid input: that number is not one of the GR-ANGLE keywords for the dictionary of science-arcs."
				else:
					break
			angleextscience = dicts.extractedsciences.get(gratingangle,'') # existence of corresponding science spectrum
			angleextstandard = dicts.extractedstandards.get(gratingangle,'') # existence of corresponding standard star spectrum
			if angleextscience == '': 
				print "The corresponding science spectrum does not exist for dispcor (used to apply the new wavelength solution): "
				print angleextscience
				print "Cannot proceed."
				continue
			if angleextstandard == '':
				print "The corresponding standard star spectrum does not exist for dispcor (used to apply the new wavelength solution): "
				print angleextstandard
				print "Skipping identification of arc lines for standard star spectrum."
				angleextarc_std = ''
				inputextarcs_std = {}
				inputstandards = {}
			elif angleextstandard != '':
				angleextarc_std = dicts.extractedarcs_std[gratingangle]
				inputextarcs_std = {gratingangle:angleextarc_std}
				inputstandards = {gratingangle:angleextstandard}
			# identify and reidentify don't produce any files; just header keywords or update the database directory
			angleextarc_sci = dicts.extractedarcs_sci[gratingangle]
			inputextarcs_sci = {gratingangle:angleextarc_sci}
			inputsciences = {gratingangle:angleextscience}
			tasks.identifyarcs(arcsci=inputextarcs_sci,scienceimages=inputsciences,indiv=True) # runs identify on science-arcs
			tasks.reidentifyarcs(arcstd=inputextarcs_std,arcsci=inputextarcs_sci,standardimages=inputstandards,indiv=True) # reidentify for standard-star-arcs
		elif choice=='17': # dispersion-correct science spectra with second wavelength solution
			print "You have chosen to dispersion-correct science spectra using the second wavelength solution."
			print "Here is the global dictionary of science spectra indexed by GR-ANGLE:"
			print dicts.extractedsciences
			while True:
				gratingangle = raw_input("Please enter the GR-ANGLE keyword of the science spectrum you wish to dispersion-correct: ")
				angles = dicts.extractedsciences.keys()
				if gratingangle not in angles: # angles[i] is already a string, not a double/float
					print "Invalid input: that number is not one of the GR-ANGLE keywords for the dictionary of science spectra."
				else:
					break
			angleextarc_sci = dicts.extractedarcs_sci.get(gratingangle,'') # existence of identified science-arc spectrum
			if angleextarc_sci == '': 
				print "The corresponding science-extracted arc spectrum does not exist: "
				print angleextarc_sci
				print "Cannot proceed."
				continue
			hdulist = pyfits.open(angleextarc_sci) # need to check if lines of science-extracted arc spectrum have been identified
			hdr = hdulist[0].header
			try:
				proc = list(hdr['RUPIPE']) 
			except:
				proc = [] # in case image has never been processed so 'RUPIPE' keyword is not in the header
			if 'i' not in proc:
				print "Lines in the corresponding science-extracted arc spectrum have not been identified."
				print "Cannot proceed without a second wavelength solution (via Option 14)."
				continue
			hdulist.close() # else, continue onwards
			angledispscience = dicts.dispsciences.get(gratingangle,'') # check for existing dispersion-corrected science
			if angledispscience != '':
				print "There is already a dispersion-corrected science spectrum for angle "+gratingangle+": "
				print angledispscience
				print "If you continue, the dispersion-corrected science spectrum for this angle will be replaced with a new one."
				print "An alternative is to work in a cloned directory (using Option 0 from the main menu)."
				while True: 
					moveforward = raw_input("Enter 0 to return to the submenu, or 1 to move forward: ")
					if moveforward == '0' or moveforward == '1':
						break
					else:
						print "Invalid input: you must enter either 0 or 1."
				if moveforward == '0':
					print "You chose to return to the submenu instead of replacing the old extracted dispersion-corrected science spectrum."
					continue
				os.remove(angledispscience)
			angleextscience = dicts.extractedsciences[gratingangle]
			inputsciences = {gratingangle:angleextscience} # already checked for REFSPEC1 and identified arcs above (via 'i' in 'RUPIPE' header keyword)
			tasks.dispcorsci(scienceimages=inputsciences,indiv=True)
		elif choice=='18': # dispersion-correct standard star spectra with second wavelength solution
			print "You have chosen to dispersion-correct standard star spectra using the second wavelength solution."
			print "Here is the global dictionary of standard star spectra indexed by GR-ANGLE:"
			print dicts.extractedstandards
			while True:
				gratingangle = raw_input("Please enter the GR-ANGLE keyword of the standard star spectrum you wish to dispersion-correct: ")
				angles = dicts.extractedstandards.keys()
				if gratingangle not in angles: # angles[i] is already a string, not a double/float
					print "Invalid input: that number is not one of the GR-ANGLE keywords for the dictionary of standard star spectra."
				else:
					break
			angleextarc_std = dicts.extractedarcs_std.get(gratingangle,'') # existence of identified standard-star-arc spectrum
			if angleextarc_std == '': 
				print "The corresponding standard-star-extracted arc spectrum does not exist: "
				print angleextarc_std
				print "Cannot proceed."
				continue
			hdulist = pyfits.open(angleextarc_std) # need to check if lines of std-star-extracted arc spectrum have been identified
			hdr = hdulist[0].header
			try:
				proc = list(hdr['RUPIPE']) 
			except:
				proc = [] # in case image has never been processed so 'RUPIPE' keyword is not in the header
			if 'i' not in proc:
				print "Lines in the corresponding standard-star-extracted arc spectrum have not been identified."
				print "Cannot proceed without a second wavelength solution (via Option 14)."
				continue
			hdulist.close() # else, continue onwards
			angledispstandard = dicts.dispstandards.get(gratingangle,'') # check for existing dispersion-corrected std star
			if angledispstandard != '':
				print "There is already a dispersion-corrected standard star spectrum for angle "+gratingangle+": "
				print angledispstandard
				print "If you continue, the dispersion-corrected standard star spectrum for this angle will be replaced with a new one."
				print "An alternative is to work in a cloned directory (using Option 0 from the main menu)."
				while True: 
					moveforward = raw_input("Enter 0 to return to the submenu, or 1 to move forward: ")
					if moveforward == '0' or moveforward == '1':
						break
					else:
						print "Invalid input: you must enter either 0 or 1."
				if moveforward == '0':
					print "You chose to return to the submenu instead of replacing the old extracted dispersion-corrected standard star spectrum."
					continue
				os.remove(angledispstandard)
			angleextstandard = dicts.extractedstandards[gratingangle]
			inputstandards = {gratingangle:angleextstandard} # already checked for REFSPEC1 and identified arcs above (via 'i' in 'RUPIPE' header keyword)
			tasks.dispcorstd(standardimages=inputstandards,indiv=True)
		elif choice=='19': # create bad wavelength masks for all dispersion-corrected science and standard star spectra
			print "You have chosen to create bad wavelength masks for dispersion-corrected science spectra."
			print "Bad wavelength masks will also automatically be created for dispersion-corrected standard star spectra should they exist."
			print "The dispersion-corrected science spectra are:"
			print dicts.dispsciences
			print "The dispersion-corrected standard star spectra are:"
			print dicts.dispstandards
			if dicts.dispsciences == {}:
				print "There are no dispersion-corrected science spectra to create bad wavelength masks for."
				print "Cannot proceed."
				continue
			if dicts.dispstandards == {}:
				print "There are no dispersion-corrected standard star spectra."
				print "Will not create bad wavelength masks for standard star images." # no dispstandards => no standards
			if dicts.bwmsciences != {} or dicts.bwmstandards != {}:
				print "There are already bad wavelength masks for science or standard star spectra: "
				print dicts.bwmsciences
				print dicts.bwmstandards
				print "If you continue, the current bad wavelength masks will be replaced with new ones."
				print "An alternative is to work in a cloned directory (using Option 0 from the main menu)."
				while True: 
					moveforward = raw_input("Enter 0 to return to the submenu, or 1 to move forward: ")
					if moveforward == '0' or moveforward == '1':
						break
					else:
						print "Invalid input: you must enter either 0 or 1."
				if moveforward == '0':
					print "You chose to return to the submenu instead of replacing the old bad wavelength masks."
					continue
			imdelseq = ''
			for bwmscience in dicts.bwmsciences:
				imdelseq = imdelseq+','+bwmscience
			for bwmstandard in dicts.bwmstandards:
				imdelseq = imdelseq+','+bwmstandard
			pyrafCalls.run_imdel(inputname=imdelseq)
			tasks.maskspectra(dispscispectra=dicts.dispsciences,dispstdspectra=dicts.dispstandards,indiv=True)
		elif choice=='20': # flux-calibrate dispersion-corrected science and standard star spectra
			print "You have chosen to flux-calibrate science and standard star spectra."
			print "Here are the global dictionaries of dispersion-corrected science and standard star spectra indexed by GR-ANGLE:"
			print dicts.fluxsciences
			print dicts.fluxstandards
			print dicts.stdfiles
			print dicts.sensfiles
			if dicts.fluxsciences != {} or dicts.fluxstandards != {} or dicts.sensfiles != {}: # check for existence of flux-calibrated science and standard star spectra
				print "There are already flux-calibrated science or standard star spectra: "
				print dicts.fluxsciences
				print dicts.fluxstandards
				print dicts.sensfiles
				print "If you continue, the flux-calibrated science or standard star spectra will be replaced with new ones."
				print "An alternative is to work in a cloned directory (using Option 0 from the main menu)."
				while True: 
					moveforward = raw_input("Enter 0 to return to the submenu, or 1 to move forward: ")
					if moveforward == '0' or moveforward == '1':
						break
					else:
						print "Invalid input: you must enter either 0 or 1."
				if moveforward == '0':
					print "You chose to return to the submenu instead of replacing the old flux-calibrated science or standard star spectra."
					continue
				for k in dicts.fluxsciences.keys():
					os.remove(dicts.fluxsciences[k])
				for k in dicts.fluxstandards.keys():
					os.remove(dicts.fluxstandards[k])
				for k in dicts.stdfiles.keys():
					try:
						os.remove(dicts.stdfiles[k])
					except:
						print "no std bandpass file for angle "+k
				print "DELETE old standard star SENSFUNC files?"
				while True: 
					moveforward = raw_input("Enter 0 to KEEP, or 1 to DELETE: ")
					if moveforward == '0' or moveforward == '1':
						break
					else:
						print "Invalid input: you must enter either 0 or 1."
				if moveforward == 1:
					for k in dicts.sensfiles.keys():
						os.remove(dicts.sensfiles[k])
			inputsciences = dicts.dispsciences
			inputstandards = dicts.dispstandards
			tasks.fluxcal(dspsci=inputsciences,dspstd=inputstandards,indiv=True)
		elif choice=='21': # sigma-clip science and standard star spectra
			print "You have chosen to sigma-clip science and standard star spectra."
			print "Here are the global dictionaries of relevant science and standard star spectra:"
			print dicts.dispsciences
			print dicts.fluxsciences
			print dicts.fluxstandards
			if dicts.dispsciences == {} and dicts.fluxsciences == {} and dicts.fluxstandards == {}:
				print "There are no science or standard star spectra to sigma-clip."
				print "Cannot proceed."
				continue
			print "WARNING: If these spectra have already been sigma-clipped, and"
			print "you try to sigma-clip them again, you could overclean the spectra."
			print "You would have to re-run flux calibration to undo the overcleaning."
			print "Note: sigma-clipped pixels are added to the bad wavelength mask, not deleted in this step."
			while True: 
				moveforward = raw_input("Enter 0 to return to the submenu, or 1 to move forward: ")
				if moveforward == '0' or moveforward == '1':
					break
				else:
					print "Invalid input: you must enter either 0 or 1."
			if moveforward == '0':
				print "You chose to return to the submenu instead of cleaning the spectra."
				continue
			tasks.sigmaClipSpectra(dispspectra=dicts.dispsciences,fluxspectra=dicts.fluxsciences,standardspectra=dicts.fluxstandards,masksciences=dicts.bwmsciences,maskstandards=dicts.bwmstandards,indiv=True)
		elif choice=='22': # scale flux science spectra across all gr-angles
			print "You have chosen to scale flux-calibrated science spectra."
			print "Here is the global dictionary of flux-calibrated science spectra indexed by GR-ANGLE:"
			print dicts.fluxsciences
			if dicts.fluxsciences == {}:
				print "There are no flux-calibrated science spectra to scale."
				print "Cannot proceed."
				continue
			if len(dicts.fluxsciences) == 1:
				print "There is only one flux-calibrated science image."
				print "Cannot proceed."
				continue
			if dicts.scaledfluxsciences != {}:
				print "There are already scaled flux-calibrated science spectra: "
				print dicts.scaledfluxsciences
				print "If you continue, the scaled flux-calibrated science spectra will be replaced with new ones."
				print "An alternative is to work in a cloned directory (using Option 0 from the main menu)."
				while True: 
					moveforward = raw_input("Enter 0 to return to the submenu, or 1 to move forward: ")
					if moveforward == '0' or moveforward == '1':
						break
					else:
						print "Invalid input: you must enter either 0 or 1."
				if moveforward == '0':
					print "You chose to return to the submenu instead of replacing the old scaled science spectra."
					continue
				for k in dicts.scaledfluxsciences.keys(): # delete each old scaled science spectrum
					os.remove(dicts.scaledfluxsciences[k])
			inputsciences = dicts.fluxsciences
			tasks.scaleFluxScienceSpectra(sciencespectra=inputsciences,indiv=True)		
		elif choice=='23': # scale standard star spectra across all gr-angles
			print "You have chosen to scale standard star spectra."
			print "Here is the global dictionary of flux-calibrated standard star spectra indexed by GR-ANGLE:"
			print dicts.fluxstandards
			if dicts.fluxstandards == {}:
				print "There are no flux-calibrated standard star spectra to scale."
				print "Cannot proceed."
				continue
			if len(dicts.fluxstandards) == 1:
				print "There is only one flux-calibrated standard star image."
				print "Cannot proceed."
				continue
			if dicts.scaledstandards != {}:
				print "There are already scaled standard star spectra: "
				print dicts.scaledstandards
				print "If you continue, the scaled standard star spectra will be replaced with new ones."
				print "An alternative is to work in a cloned directory (using Option 0 from the main menu)."
				while True: 
					moveforward = raw_input("Enter 0 to return to the submenu, or 1 to move forward: ")
					if moveforward == '0' or moveforward == '1':
						break
					else:
						print "Invalid input: you must enter either 0 or 1."
				if moveforward == '0':
					print "You chose to return to the submenu instead of replacing the old scaled standard star spectra."
					continue
				for k in dicts.scaledstandards.keys(): # delete each old scaled standard star spectrum
					os.remove(dicts.scaledstandards[k])
			inputsciences = dicts.fluxstandards
			tasks.scaleStdStarSpectra(stdspectra=dicts.fluxstandards,indiv=True)
		elif choice=='24': # copy the sky and sigma spectral bands (3 and 4) of flux-cal science spectra into separate fits files
			print "You have chosen to copy the sky and sigma spectral bands into separate spectral files."
			print "Here is the global dictionary of flux-calibrated science spectra indexed by GR-ANGLE:"
			print dicts.fluxsciences
			if dicts.fluxsciences == {}: # check for emptiness of fluxsciences 
				print "Cannot proceed: the flux-calibrated science spectra dictionary is empty."
				continue
			# check for existing skysciences and sigmasciences
			if dicts.skysciences != {} or dicts.sigmasciences != {}:
				print "There are already sky and sigma spectra files: "
				print dicts.skysciences
				print dicts.sigmasciences
				print "If you continue, any sky and sigma spectra will be replaced with new ones."
				print "An alternative is to work in a cloned directory (using Option 0 from the main menu)."
				while True: 
					moveforward = raw_input("Enter 0 to return to the submenu, or 1 to move forward: ")
					if moveforward == '0' or moveforward == '1':
						break
					else:
						print "Invalid input: you must enter either 0 or 1."
				if moveforward == '0':
					print "You chose to return to the submenu instead of replacing the old sky and sigma spectra."
					continue
				for k in dicts.skysciences.keys():
					os.remove(dicts.skysciences[k])
				for j in dicts.sigmasciences.keys():
					os.remove(dicts.sigmasciences[j])
			tasks.copySciSkySig(fluxspectra=dicts.fluxsciences,indiv=True)
		elif choice=='25': # scale disp science spectra across all gr-angles
			print "You have chosen to scale dispersion-corrected science spectra."
			print "Here is the global dictionary of dispersion-corrected science spectra indexed by GR-ANGLE:"
			print dicts.dispsciences
			if dicts.dispsciences == {}:
				print "There are no dispersion-corrected science spectra to scale."
				print "Cannot proceed."
				continue
			if len(dicts.dispsciences) == 1:
				print "There is only one dispersion-corrected science image."
				print "Cannot proceed."
				continue
			if dicts.scaleddispsciences != {}:
				print "There are already scaled dispersion-corrected science spectra: "
				print dicts.scaleddispsciences
				print "If you continue, the scaled dispersion-corrected  science spectra will be replaced with new ones."
				print "An alternative is to work in a cloned directory (using Option 0 from the main menu)."
				while True: 
					moveforward = raw_input("Enter 0 to return to the submenu, or 1 to move forward: ")
					if moveforward == '0' or moveforward == '1':
						break
					else:
						print "Invalid input: you must enter either 0 or 1."
				if moveforward == '0':
					print "You chose to return to the submenu instead of replacing the old scaled science spectra."
					continue
				for k in dicts.scaleddispsciences.keys(): # delete each old scaled science spectrum
					os.remove(dicts.scaleddispsciences[k])
			inputsciences = dicts.dispsciences
			tasks.scaleDispScienceSpectra(sciencespectra=inputsciences,indiv=True)
		elif choice=='26': # odcombine dispsciences and fluxsciences (separately)
			print "You have chosen to combine dispersion-corrected and flux-calibrated spectra."
			print "Here are the global dictionaries of relevant spectra indexed by GR-ANGLE:"
			print dicts.dispsciences
			print dicts.fluxsciences
			print dicts.scaleddispsciences
			print dicts.scaledfluxsciences
			print dicts.fluxstandards
			print dicts.scaledstandards
			if dicts.dispsciences == {} and dicts.fluxsciences == {}: # check for emptiness of at least dispsciences and fluxsciences dictionaries 
				print "The dispersion-corrected AND flux-calibrated science spectra dictionaries are both empty: "
				print dicts.dispsciences
				print dicts.fluxsciences
				print "Cannot proceed."
				continue
			combdispscience = dicts.combinedspectra.get('dsp','') # check for existing combined flux or disp spectra
			combfluxscience = dicts.combinedspectra.get('flx','')
			combsigscience = dicts.combinedspectra.get('sig','')
			combskyscience = dicts.combinedspectra.get('sky','')
			combfluxstandard = dicts.combinedspectra.get('std','')
			if combdispscience != '' or combfluxscience != '' or combsigscience != '' or combskyscience != '' or combfluxstandard != '':
				print "There are already combined spectra: "
				print combdispscience, combfluxscience, combsigscience, combskyscience, combfluxstandard
				print "If you continue, any combined spectra will be replaced with new ones."
				print "An alternative is to work in a cloned directory (using Option 0 from the main menu)."
				while True: 
					moveforward = raw_input("Enter 0 to return to the submenu, or 1 to move forward: ")
					if moveforward == '0' or moveforward == '1':
						break
					else:
						print "Invalid input: you must enter either 0 or 1."
				if moveforward == '0':
					print "You chose to return to the submenu instead of replacing the old combined spectra."
					continue
				try:
					os.remove(combdispscience)
				except:
					print "Did not delete combined dispersion-corrected science spectrum because it does not exist."
				try:
					os.remove(combfluxscience)
				except:	
					print "Did not delete combined flux-calibrated science spectrum because it does not exist."
				try:
					os.remove(combsigscience)
				except:	
					print "Did not delete combined sigma (science-based) spectrum because it does not exist."
				try:
					os.remove(combskyscience)
				except:	
					print "Did not delete combined sky (science-based) spectrum because it does not exist."
				try:
					os.remove(combfluxstandard)
				except:	
					print "Did not delete combined flux-calibrated standard star spectrum because it does not exist."
			tasks.combinespectra(dispspectra=dicts.dispsciences,fluxspectra=dicts.fluxsciences,fluxstd=dicts.fluxstandards,sclstd=dicts.scaledstandards,indiv=True) # input spectra are just the original, complete dictionaries
		elif choice=='27': # clean combined sci and standard star spectra
			print "You have chosen to clean up combined science and standard star spectra."
			print "Here are the global dictionaries of relevant spectra indexed by GR-ANGLE:"
			print dicts.combinedspectra
			if dicts.combinedspectra == {}: # check for lack of combined sci and std star spectra 
				print "The combined spectra dictionary is empty: "
				print dicts.combinedspectra
				print "Cannot proceed."
				continue
			clnsci = dicts.combinedspectra.get('clnsci','') # check for existing cleaned flux or disp spectra
			clnstd = dicts.combinedspectra.get('clnstd','')
			if clnsci != '' or clnstd != '':
				print "There are already cleaned spectra: "
				print clnsci,clnstd
				print "If you continue, any cleaned spectra will be replaced with new ones."
				print "An alternative is to work in a cloned directory (using Option 0 from the main menu)."
				while True: 
					moveforward = raw_input("Enter 0 to return to the submenu, or 1 to move forward: ")
					if moveforward == '0' or moveforward == '1':
						break
					else:
						print "Invalid input: you must enter either 0 or 1."
				if moveforward == '0':
					print "You chose to return to the submenu instead of replacing the old cleaned spectra."
					continue
				try:
					os.remove(clnsci)
				except:
					print "Did not delete cleaned science spectrum because it does not exist."
				try:
					os.remove(clnstd)
				except:	
					print "Did not delete cleaned standard star spectrum because it does not exist."
			tasks.cleanCombinedSpectra(combspectra=dicts.combinedspectra,indiv=True)
		elif choice=='28': # telluric-correct combined sci and standard star spectra
			print "You have chosen to telluric-correct a science spectrum."
			print "Here are the global dictionaries of relevant spectra indexed by GR-ANGLE:"
			print dicts.combinedspectra
			clnsci = dicts.combinedspectra.get('clnsci','') # check for existing cleaned flux or disp spectra
			clnstd = dicts.combinedspectra.get('clnstd','')
			if clnsci == '' or clnstd == '': # check for lack of combined sci and std star spectra 
				print "There is no cleaned and combined science or standard star spectrum: "
				print clnsci,clnstd
				print "Cannot proceed."
				continue
			telsci = dicts.combinedspectra.get('tel','') # check for existing telluric-corrected spectrum
			if telsci != '':
				print "There is already a telluric-corrected spectrum: "
				print telsci
				print "If you continue, this spectrum will be replaced with a new one."
				print "An alternative is to work in a cloned directory (using Option 0 from the main menu)."
				while True: 
					moveforward = raw_input("Enter 0 to return to the submenu, or 1 to move forward: ")
					if moveforward == '0' or moveforward == '1':
						break
					else:
						print "Invalid input: you must enter either 0 or 1."
				if moveforward == '0':
					print "You chose to return to the submenu instead of replacing the old spectrum."
					continue
				try:
					os.remove(telsci)
				except:
					print "Did not delete telluric-corrected science spectrum because it does not exist."
			tasks.telluricCorrect(combspectra=dicts.combinedspectra,indiv=True)
		elif choice=='29':
			print "You chose to return to the main menu."
			break
		else:
			print 'You must pick an option from the menu.'
			continue 

