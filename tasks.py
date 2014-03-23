
'''
>>> Rutgers SALT Supernova Spectral Reduction Pipeline <<<

This is the CORE MODULE of the pipeline. It contains the pipeline 
functions which are used to prepare the spectral data which will 
then be fed to PyRAF for reduction, extraction, and calibration.

Each function is relatively self-contained in the sense that
it represents only one step of the entire reduction process. 
Naturally, then, each function has within it one main call
corresponding to the relevant PyRAF task.

Please refer to the documentation for more information
about how exactly each function prepares the data for PyRAF.

*** Modifications ***
Sept. 26, 2013: Created module. -Viraj Pandya

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
import pipeHistory # Responsible for updating FITS headers with pipeline-processed indicator keywords.

# This function runs flatcombine on flats with the same gr-angle (and thus same exptime)
# default input is the dictionary flats created in the function dictionaries()
# argument dictionary should have gr-angle as keywords, and lists of image file names as the values of those keywords
def combineflats(allflats=dicts.flats,indiv=False):
	# This calls dictionaries() if there does not exist the dictionary flats.
	try:
		len(allflats)
	except:
		dictionaries()
		allflats=dicts.flats
	# This will run flatcombine on each list of flats in the dictionary flats indexed by gr-angle.
	# Making an assumption about the flats in each list: same gr-angle => same exptime.
	angles = allflats.keys()
	for angle in angles:
		flatlist = allflats[angle]
		print "the flatlist for "+angle+" is: "
		print flatlist
		# set output names
		name = 'flt'+str(angle)+'cmb.fits'
		auxname = 'flt'+str(angle)+'cmbAUX.fits'
		# If flatlist doesn't have multiple flats in it, can't "combine", just imcopy and use new one as combined flat
		if len(flatlist)==0:
			print "WARNING: no flat for angle"+angle
			continue
		elif len(flatlist)==1:
			print "WARNING: only one flat for angle"+angle
			pyrafCalls.run_imcopy(inputname=flatlist[0],outputname=name) # will make copy, preserve 2-extension file structure
			pipeHistory.updatePipeKeys(inputname=name,imagetype='flat',procChar='c')
			dicts.combflats[angle] = name
			continue
		# Loop through the headers of each flat in the list to find the average EXPTIME and AIRMASS
		avgEXPTIME = 0
		avgAIRMASS = 0
		for flat in flatlist:
			print flat
			hdr = pyfits.getheader(flat,0)
			avgEXPTIME = avgEXPTIME + abs(hdr['EXPTIME'])
			avgAIRMASS = avgAIRMASS + abs(hdr['AIRMASS'])
			hdulist = pyfits.open(flat,mode='update')
			hdrold = hdulist[0].header.copy()
			hdulist.close()
		avgEXPTIME = avgEXPTIME/len(flatlist)
		avgAIRMASS = avgAIRMASS/len(flatlist)
		# this calls the run_flatcombine(flatstocombine,combflatname) function
		pyrafCalls.run_flatcombine(flatstocombine=flatlist,combflatname=auxname,customRun=indiv)
		# This copies one original file to the final file (to preserve file structure), and replaces data and header keys
		pyrafCalls.run_imcopy(inputname=flatlist[0],outputname=name)
		hdulist = pyfits.open(auxname)
		datnew = hdulist[0].data.copy()
		hdrnew = hdulist[0].header.copy()
		hdulist.close()
		hdunew = pyfits.open(name,mode='update')
		hdunew[1].data = datnew
		hdunew[1].header = hdrnew
		hdunew.flush()
		hdr = hdunew[0].header
		hdr['EXPTIME'] = avgEXPTIME
		hdr['AIRMASS'] = avgAIRMASS
		hdr1 = hdunew[1].header
		hdr1['EXTVER'] = 1
		hdr1['EXTNAME'] = 'SCI'
		hdunew.flush()
		hdunew.close()
		pipeHistory.updatePipeKeys(inputname=name,imagetype='flat',procChar='c')
		# This adds the newly created combined flat name to the combflats dictionary indexed by gr-angle
		dicts.combflats[angle] = name
		# this deletes the aux image file
		pyrafCalls.run_imdel(auxname)
		# This informs the user of the newly created combined flat-field image
		print 'output: '+name
	# printing new dictionaries
	print 'The combined flat image filenames are:'
	print dicts.combflats

# This function normalizes the combined flats in the combflats dictionary.
def normalizeflats(combinedflats=dicts.combflats,indiv=False):
	# This calls combineflats() if there does not exist the dictionary combflats.
	try:
		len(combinedflats)
	except:
		combineflats()
		combinedflats=dicts.combflats
	# This runs response on each combined flat and adds the normalized flat to the normflats dictionary indexed by gr-angle.
	angles = combinedflats.keys()
	for angle in angles:
		combflat = combinedflats[angle]
		hdulist = pyfits.open(combflat,mode='update')
		hdr = hdulist[0].header
		hdr.update('DISPAXIS',1)
		hdulist.flush()
		hdrold = hdulist[0].header.copy()
		hdulist.close()
		# output names
		name = 'flt'+str(angle)+'nrm.fits'
		auxname = 'flt'+str(angle)+'nrmAUX.fits'
		# calls the run_response function below
		pyrafCalls.run_response(combinedflat=combflat+'[1]',normflatname=auxname,customRun=indiv)
		# This copies one original file to the final file (to preserve file structure), and replaces data and header keys
		pyrafCalls.run_imcopy(inputname=combflat,outputname=name)
		hdulist = pyfits.open(auxname)
		datnew = hdulist[0].data.copy()
		hdrnew = hdulist[0].header.copy()
		hdulist.close()
		hdunew = pyfits.open(name,mode='update')
		hdunew[1].data = datnew
		hdunew[1].header = hdrnew
		hdunew.flush()
		hdr = hdunew[0].header
		hdr1 = hdunew[1].header
		hdr1['EXTVER'] = 1
		hdr1['EXTNAME'] = 'SCI'
		hdunew.flush()
		hdunew.close()
		pipeHistory.updatePipeKeys(inputname=name,imagetype='flat',procChar='n')
		# Adds the name of the newly created normalized flat to the normflats dictionary indexed by gr-angle
		dicts.normflats[angle] = name
		# this deletes the aux image file
		pyrafCalls.run_imdel(auxname)
		print 'output: '+name
	print 'The normalized flat image filenames are:'
	print dicts.normflats

# This function divides each science image by the normalized flat with the same gr-angle.
def flattensciences(scienceimages=dicts.sciences,normalizedflats=dicts.normflats,indiv=False):
	# This prints an error message and quits function if either of input dictionaries don't exist
	try:
		len(normalizedflats),len(scienceimages)
	except:
		print 'An inputted dictionary does not exist: cannot proceed.'
		return
	# This runs imarith to divide a science image by the normalized flat with the same gr-angle.
	angles = scienceimages.keys()
	for angle in angles:
		science = scienceimages[angle] # this returns the list of all science images with GR-ANGLE=angle
		suffixlist = science.split('.')
		suffix = suffixlist[0][15:]
		# this stores the old headers (extension 0 and 1)
		hdulist = pyfits.open(science,mode='update')
		hdr0_old = hdulist[0].header.copy()
		hdr1_old = hdulist[1].header.copy()
		hdulist.close()
		normflat = normalizedflats.get(angle,-1)
		# if no normflat with same gr-angle found, continue onto next image to be flat-fielded
		if normflat == -1:
			print 'no normalized flat found for '+science+' with angle: '+angle
			print 'skipping flat-fielding of '+science
			continue
		# output names
		namesci = 'sci'+str(angle)+'flt'+suffix+'.fits'
		nameaux = 'sci'+str(angle)+'flt'+suffix+'AUX'+'.fits'
		# run iraf.images.imutil.imarith to simply divide science by normflat
		pyrafCalls.run_imarith(dividend=science,divisor=normflat+'[1]',quotient=nameaux,customRun=indiv)
		# This copies the original file to the final file (to preserve file structure), and replaces data
		pyrafCalls.run_imcopy(inputname=science,outputname=namesci)
		hdulist = pyfits.open(nameaux)
		datnew = hdulist[0].data.copy()
		hdulist.close()
		hdunew = pyfits.open(namesci,mode='update')
		hdunew[1].data = datnew
		hdunew.flush()
		hdunew.close()
		pipeHistory.updatePipeKeys(inputname=namesci,imagetype='science',procChar='f')
		# This adds the newly flat-fielded science image's name to the dictionary flatsciences.
		dicts.flatsciences[angle] = namesci
		# this deletes the aux image file
		pyrafCalls.run_imdel(nameaux)
		# updating user about successful flat-fielding
		print 'output: '+namesci
	print 'The flat-fielded science image filenames are:'
	print dicts.flatsciences

# This function divides each arc image by the normalized flat with the same gr-angle.
def flattenarcs(arcimages=dicts.arcs,normalizedflats=dicts.normflats,indiv=False):
	# This prints an error message and quits the function if one of the inputted dictionary doesn't exist
	try:
		len(normalizedflats),len(arcimages)
	except:
		print 'An inputted dictionary does not exist: cannot proceed'
		return
	# This runs imarith to divide an arc image by the normalized flat with the same gr-angle.
	angles = arcimages.keys()
	for angle in angles:
		arc = arcimages[angle]
		# this copies (and stores) the old (original) arc image's header
		suffixlist = arc.split('.')
		suffix = suffixlist[0][15:]
		hdulist = pyfits.open(arc,mode='update')
		hdr0_old = hdulist[0].header.copy()
		hdr1_old = hdulist[1].header.copy()
		hdulist.close()
		normflat = normalizedflats.get(angle,-1)
		# if no normflat with same gr-angle found, continue onto next image to be flat-fielded
		if normflat == -1:
			print 'WARNING: no normalized flat found for '+arc+' with angle: '+angle
			print 'Skipping flat-fielding of arc: '+arc
			continue
		# output name
		namearc = 'arc'+str(angle)+'flt'+suffix+'.fits'
		nameaux = 'arc'+str(angle)+'flt'+suffix+'AUX'+'.fits'
		# run iraf.images.imutil.imarith to simply divide science by normflat
		pyrafCalls.run_imarith(dividend=arc,divisor=normflat+'[1]',quotient=nameaux,customRun=indiv)
		# This copies the original file to the final file (to preserve file structure), and replaces data
		pyrafCalls.run_imcopy(inputname=arc,outputname=namearc)
		hdulist = pyfits.open(nameaux)
		datnew = hdulist[0].data.copy()
		hdulist.close()
		hdunew = pyfits.open(namearc,mode='update')
		hdunew[1].data = datnew
		hdunew.flush()
		hdunew.close()
		pipeHistory.updatePipeKeys(inputname=namearc,imagetype='arc',procChar='f')
		# This adds the newly flat-fielded arc image's name to the dictionary flatarcs.
		dicts.flatarcs[angle] = namearc
		# this deletes the aux image files
		pyrafCalls.run_imdel(nameaux)
		print 'output: '+namearc
	print 'The flat-fielded arc image filenames are:'
	print dicts.flatarcs

# This function divides each standard star image by the normalized flat with the same gr-angle.
def flattenstandards(standardimages=dicts.standards,normalizedflats=dicts.normflats,indiv=False):
	# This prints an error message and quits the function if an inputted dictionary doesn't exist
	try:
		len(normalizedflats),len(standardimages)
	except:
		print 'An inputted dictionary does not exist: cannot proceed.'
		return
	# This runs imarith to divide a standard star image by the normalized flat with the same gr-angle.
	angles = standardimages.keys()
	for angle in angles:
		standardlist = standardimages[angle] # this returns the list of all standard images with GR-ANGLE=angle
		if type(standardlist) == str:
			standardsWereCombined = False
		elif type(standardlist)==list and len(standardlist)==1:
			standardimages[angle]=standardlist[0] # sets standardimages[angle] equal to the filename of the sole standard image
			standardsWereCombined = False
		else:
			print("MULTIPLE STANDARDS for angle "+str(angle),standardlist)
			avgAIRMASS = 0 # will find the average AIRMASS
			suffix = ''
			for img in standardlist: # create suffix sequence of standard image numbers and avgAIRMASS
				suffixlist = img.split('.')
				suffix = suffix+suffixlist[0][15:]
				hdr = pyfits.getheader(img,0)
				avgAIRMASS = avgAIRMASS + abs(hdr['AIRMASS'])
				saltgain = hdr.get('GAIN',2.443)
			avgAIRMASS = avgAIRMASS/len(standardlist) # dividing sum of airmasses by number of images
			combstdname = 'std'+str(angle)+'cmb'+suffix+'.fits'
			auxname = 'std'+str(angle)+'cmb'+suffix+'AUX.fits'
			pyrafCalls.run_imcopy(inputname=standardlist[0],outputname=combstdname)
			for i,v in enumerate(standardlist): # need to append [1] to each filename since imcombine will need that
				standardlist[i] = v+'[1]'
			standardseq = ','.join(standardlist)
			pyrafCalls.run_imcombine(imagestocombine=standardseq,combimgname=auxname,commongain=saltgain,customRun=indiv)
			hdulist = pyfits.open(auxname)
			hdr = hdulist[0].header
			datnew = hdulist[0].data.copy()
			avgEXPTIME = hdr['EXPTIME']
			hdulist.close()
			hdunew = pyfits.open(combstdname,mode='update')
			hdr0 = hdunew[0].header
			hdr1 = hdunew[1].header
			hdunew[1].data = datnew
			hdr0['EXPTIME'] = avgEXPTIME
			hdr1['EXPTIME'] = avgEXPTIME
			hdr0['AIRMASS'] = avgAIRMASS
			hdr1['AIRMASS'] = avgAIRMASS
			hdunew.flush()
			hdunew.close()
			standardimages[angle] = combstdname # replaces the list of standards with the combined standard
			standardsWereCombined = True
			# this deletes the aux image file
			pyrafCalls.run_imdel(auxname)
		standard = standardimages[angle]
		suffixlist = standard.split('.')
		if standardsWereCombined == False:
			suffix = suffixlist[0][15:]
		elif standardsWereCombined == True:
			suffix = suffixlist[1][5:]
		hdulist = pyfits.open(standard,mode='update')
		hdr0_old = hdulist[0].header.copy()
		hdr1_old = hdulist[1].header.copy()
		hdulist.close()
		# choose normflat with same rounded gr-angle (string format)
		normflat = normalizedflats.get(angle,-1)
		# if no normflat with same gr-angle found, continue onto next image to be flat-fielded
		if normflat == -1:
			print 'no normalized flat found for '+standard+' with angle: '+angle
			continue
		# output name
		namestd = 'std'+str(angle)+'flt'+suffix+'.fits'
		nameaux = 'std'+str(angle)+'flt'+suffix+'AUX'+'.fits'
		# run iraf.images.imutil.imarith to simply divide standard by normflat
		pyrafCalls.run_imarith(dividend=standard,divisor=normflat+'[1]',quotient=nameaux,customRun=indiv)
		# This copies the original file to the final file (to preserve file structure), and replaces data
		pyrafCalls.run_imcopy(inputname=standard,outputname=namestd)
		hdulist = pyfits.open(nameaux)
		datnew = hdulist[0].data.copy()
		hdulist.close()
		hdunew = pyfits.open(namestd,mode='update')
		hdunew[1].data = datnew
		hdunew.flush()
		hdunew.close()
		pipeHistory.updatePipeKeys(inputname=namestd,imagetype='standard',procChar='f')
		# This adds the newly flat-fielded standard star image's name to the dictionary flatstandards.
		dicts.flatstandards[angle] = namestd
		# this deletes the aux image files
		pyrafCalls.run_imdel(nameaux)
		print 'output: '+namestd
	print 'The flat-fielded standard star image filenames are:'
	print dicts.flatstandards

# This function creates the bad pixel masks (chip gaps) for all (flat-fielded) images.
def maskimages(scienceimages=dicts.flatsciences,arcimages=dicts.flatarcs,standardimages=dicts.flatstandards,indiv=False):
	# This prints an error message and quits function if either of input dictionaries don't exist
	try:
		len(scienceimages),len(standardimages),len(arcimages)
	except:
		print 'An inputted dictionary does not exist: cannot proceed.'
		return
	if scienceimages == {}:
		print "ERROR: there are no flat-fielded science images."
		print "Skipping creation of bad pixel masks for science images."		
	if arcimages == {}:
		print "ERROR: there are no flat-fielded arc images."
		print "Skipping creation of bad pixel masks for arc images."
	if standardimages == {}:
		print "ERROR: there are no flat-fielded standard star images."
		print "Skipping creation of bad pixel masks for standard star images."
	# This runs imexpr on each science, standard star, and arc image for all gr-angles. 
	angles = scienceimages.keys() # flat-fielded arcs and standards should have same angles as science images
	for angle in angles:
		flatscience = scienceimages.get(angle,'')
		if flatscience == '': # no flatscience => no flat for this angle (usually)
			print "There is no flat-fielded science image for angle: "+angle
			print "Skipping creation of science bad pixel mask for this angle."
		else:
			suffixlist = flatscience.split('.')
			suffix = suffixlist[1][5:]
			namebpmsci = 'sci'+str(angle)+'bpm'+suffix+'.fits'
			imutil.imexpr.unlearn()
			imutil.imexpr(expr='a==0',output=namebpmsci,a=flatscience+'[1]')
			pipeHistory.updatePipeKeys(inputname=namebpmsci,imagetype='bpm_sci',procChar='o')
			hdu = pyfits.open(flatscience,mode='update')
			hdr = hdu[0].header
			hdr['RUBPM'] = namebpmsci
			hdr['GR-ANGLE'] = float(angle)
			hdu.flush()
			hdu.close()
			hdu=pyfits.open(namebpmsci,mode='update')
			hdr=hdu[0].header
			hdr['GR-ANGLE'] = float(angle)
			hdu.flush()
			hdu.close()
			dicts.bpmsciences[angle] = namebpmsci
		flatarc = arcimages.get(angle,'')
		if flatarc == '': # no flatarc => no flat for this angle
			print "There is no flat-fielded arc image for angle: "+angle
			print "Skipping creation of arc bad pixel mask for this angle."
		else:
			suffixlist = flatarc.split('.')
			suffix = suffixlist[1][5:]
			namebpmarc = 'arc'+str(angle)+'bpm'+suffix+'.fits'
			imutil.imexpr.unlearn()
			imutil.imexpr(expr='a==0',output=namebpmarc,a=flatarc+'[1]')
			pipeHistory.updatePipeKeys(inputname=namebpmarc,imagetype='bpm_arc',procChar='o')
			hdu = pyfits.open(flatarc,mode='update')
			hdr = hdu[0].header
			hdr['RUBPM'] = namebpmarc
			hdr['GR-ANGLE'] = float(angle)
			hdu.flush()
			hdu.close()
			hdu=pyfits.open(namebpmarc,mode='update')
			hdr=hdu[0].header
			hdr['GR-ANGLE'] = float(angle)
			hdu.flush()
			hdu.close()
			dicts.bpmarcs[angle] = namebpmarc
		flatstandard = standardimages.get(angle,'')
		if flatstandard == '': # no flatstandard => no flat for this angle or no standard for this angle
			print "There is no flat-fielded standard star image for angle: "+angle
			print "Skipping creation of standard star bad pixel mask for this angle."
		else:
			suffixlist = flatstandard.split('.')
			suffix = suffixlist[1][5:]
			namebpmstd = 'std'+str(angle)+'bpm'+suffix+'.fits'
			imutil.imexpr.unlearn()
			imutil.imexpr(expr='a==0',output=namebpmstd,a=flatstandard+'[1]')
			pipeHistory.updatePipeKeys(inputname=namebpmstd,imagetype='bpm_std',procChar='o')
			hdu = pyfits.open(flatstandard,mode='update')
			hdr = hdu[0].header
			hdr['RUBPM'] = namebpmstd
			hdr['GR-ANGLE'] = float(angle)
			hdu.flush()
			hdu.close()
			hdu=pyfits.open(namebpmstd,mode='update')
			hdr=hdu[0].header
			hdr['GR-ANGLE'] = float(angle)
			hdu.flush()
			hdu.close()
			dicts.bpmstandards[angle] = namebpmstd
			

# This function runs lacosmicx on (flat-fielded) science images to remove cosmic rays.
def LAxsciences(scienceimages=dicts.flatsciences,bpmasks=dicts.bpmsciences,indiv=False):
	# This prints an error message and quits function if either of input dictionaries don't exist
	try:
		len(scienceimages),len(bpmasks)
	except:
		print 'An inputted dictionary does not exist: cannot proceed.'
		return
	if scienceimages == {}:
		print "ERROR: there are no flat-fielded science images."
		print "Skipping Lacosmicx for science images."
		return
	# This runs lacosmicx on each science image
	angles = scienceimages.keys()
	for angle in angles:
		flatscience = scienceimages.get(angle,'')
		if flatscience == '': # no flatsience => no flat for this angle
			print "ERROR: there is no flat-fielded science image for angle: "+angle
			print "Skipping Lacosmicx for this angle."
			continue
		hdu = pyfits.open(flatscience)
		hdr = hdu[0].header
		bpmscience = hdr.get('RUBPM','')
		if bpmscience == '':
			print "WARNING: no bad pixel mask for science angle: "+angle
		hdu.close()
		suffixlist = flatscience.split('.')
		suffix = suffixlist[1][5:]
		namesci = 'sci'+str(angle)+'lax'+suffix+'.fits'
		namebpmlaxsci = 'sci'+str(angle)+'cpm'+suffix+'.fits' # name for the bad pixel mask that lacosmicx will produce
		# this stores the NumPy data arrays into variables to feed to lacosmicx.run()
		hduimg = pyfits.open(flatscience)
		dataimg = hduimg[1].data.copy()
		mode = imutil.imstat(images=flatscience+'[1]',fields='mode',format='no',Stdout=1) # returns list with mode of flatscience as string
		mode = float(mode[0])
		hdr0 = hduimg[0].header
		saltgain = hdr0.get('GAIN',2.443) # get the gain reported by salt or added by us in dictionaries()
		hduimg.close()
		hdubpm = pyfits.open(bpmscience)
		databpm = hdubpm[0].data.copy() # pipeline doesn't bother with the file structure (FITS extensions) of bad pixel masks
		hdubpm.close()
		# this calls the run() function in the lacosmicx package and stores the output NumPy array into a variable
		# there is no run_lacosmicx function in this pipeline; this command includes the optimal parameters for lacosmicx SNe spectra
		'''
		December 7, 2013: added robust=True; prevents nuking of sky emission lines in 4x4 images, leaves some cosmic rays in 2x4 images
		check for evidence of leftover cosmic rays, nuked galaxy and sky emission lines in extracted (and 2-D) spectra
		'''
		datalax = lacosmicx.run(inmat=dataimg,inmask=databpm,outmaskfile=namebpmlaxsci,sigclip=6.0,objlim=3.0,sigfrac=0.1,gain=saltgain,pssl=mode,robust=True)
		# This copies the original file to the final file (to preserve file structure), and replaces data
		pyrafCalls.run_imcopy(inputname=flatscience,outputname=namesci)
		hdunew = pyfits.open(namesci,mode='update')
		hdunew[1].data = datalax
		hdr = hdunew[0].header
		hdr['RUCPM'] = namebpmlaxsci # name of cosmic (ray) pixel mask
		hdr['GR-ANGLE'] = float(angle)
		hdunew.flush()
		hdunew.close()
		hdu=pyfits.open(namebpmlaxsci,mode='update')
		hdr=hdu[0].header
		hdr['GR-ANGLE'] = float(angle)
		hdu.flush()
		hdu.close()
		pipeHistory.updatePipeKeys(inputname=namesci,imagetype='science',procChar='l')
		# this adds the lacosmicx-corrected science image to the dictionary laxsciences
		dicts.laxsciences[angle] = namesci
		# this adds the 'RUIMGTYP' and 'RUPIPE' keywords to the lacosmicx pixel mask and adds it to the dictionary bpmlaxsciences
		pipeHistory.updatePipeKeys(inputname=namebpmlaxsci,imagetype='bpm_sci',procChar='l')
		dicts.bpmlaxsciences[angle] = namebpmlaxsci
	print "The LAcosmicx-corrected science image filenames are:"
	print dicts.laxsciences


# This function runs saltspec.specidentify to produce wavelength solution files.
# This doesn't yet create 2 separate dictionaries -- 1 for the arcs associated with science images, and 1 for
# the arcs associated with standard star images -- partly because SALT doesn't do that (yet).
# 2014-03-18: this function now also rectifies the identified arc, opens it in ds9, lets user re-run specidentify if bad rectification
def specidentifyarcs(arcimages=dicts.flatarcs,indiv=False):
	if len(arcimages) == 0:
		print "WARNING: there are no flat-fielded arcs. Skipping specidentify step."
		return
	# this prepares for and runs saltspec.specidentify (from pysalt)
	angles = dicts.arcs.keys() # use the keys from the product dictionary since we want to use all arcs (even non-flat-fielded ones)
	if indiv==True: # if running this as an individual task, only use keys in arcimages dictionary
		angles = arcimages.keys()
	for angle in angles:
		try: # error handling if the arc for this angle was not flat-fielded (for example, no flat for this angle)
			flatarc = arcimages[angle]
		except: # if no flat arc exists for this angle, use the product arc
			print "WARNING: No flat-fielded arc image exists for this angle: "+angle
			print "Skipping line identification."
			continue
		suffixlist = flatarc.split('.')
		suffix = suffixlist[1][5:]
		# this gets the arc lamp type so the correct linelist can be used by saltspec.specidentify
		hdulist = pyfits.open(flatarc,mode='update')
		hdr = hdulist[0].header
		lamp = hdr['LAMPID']
		hdr0_old = hdulist[0].header.copy()
		hdr1_old = hdulist[1].header.copy()
		hdulist.close()
		# this sets the arguments for run_specidentify function below
		linelistpath = params.lineListPath
		# picking the linelist file based on header LAMPID reference variable lamp
		if lamp == 'Th Ar':
			print 'the lamp is '+lamp+' for angle '+angle
			linelistpath = params.ThAr
		elif lamp == 'Xe':
			print 'the lamp is '+lamp+' for angle '+angle
			linelistpath = params.Xe
		elif lamp == 'Ne': # might instead need to use NeAr.salt
			print 'the lamp is '+lamp+' for angle '+angle
			linelistpath = params.Ne
		elif lamp == 'Cu Ar':
			print 'the lamp is '+lamp+' for angle '+angle
			linelistpath = params.CuAr
		elif lamp == 'Ar':
			print 'the lamp is '+lamp+' for angle '+angle
			linelistpath = params.Ar
		elif lamp == 'Hg Ar':
			print 'the lamp is '+lamp+' for angle '+angle
			linelistpath = params.HgAr
		else:
			print 'Could not find the proper linelist for '+lamp+' lamp.'
			continue
		while True:
			# this sets the name for the output wavelength solution file which will be produced by saltspec.specidentify
			idfilesci = 'arc'+str(angle)+'sol'+suffix+'.db'
			# idfilestd = 'arcid'+str(angle)+'_std.db'
			# this calls the run_specidentify function found further below
			pyrafCalls.run_specidentify(arcimage=flatarc,lamplines=linelistpath,idfile=idfilesci,customRun=indiv)
			# this updates the header of flatarc with 'RUPIPE' and 'RUIMGTYP' and idfile name
			hdulist = pyfits.open(flatarc,mode='update')
			hdr = hdulist[0].header
			hdr['RUSPECID'] = idfilesci
			hdulist.flush()
			hdulist.close()
			pipeHistory.updatePipeKeys(inputname=flatarc,imagetype='arc',procChar='a')
			# this adds the newly created wavelength solution file to the wavesols dictionary
			dicts.wavesols[angle] = idfilesci
			# run specrectify on flatarc
			# this sets the output filename for the 2-D wavelength-corrected arc image
			namearc = 'arc'+str(angle)+'wav'+suffix+'.fits'
			# this calls the run_specrectify function found further below
			pyrafCalls.run_specrectify(input=flatarc,output=namearc,idfile=idfilesci,customRun=indiv)		
			# this adds the 'RUPIPE' and 'RUIMGTYP' keywords to namearc's 0-header
			pipeHistory.updatePipeKeys(inputname=namearc,imagetype='arc',procChar='r')
			# this opens the rectified arc in ds9, asks user to re-run or proceed to next arc
			print namearc+' has been opened in ds9 so you can verify that rectification was successful.'
			d = ds9.ds9(start=True)
			try:
				d.set('file '+namearc)
				d.set('zoom to fit')
				d.set('zscale')
			except: # Error most likely because user closed ds9 instance.
				d = ds9.ds9(start=True)
				d.set('file '+namearc)
				d.set('zoom to fit')
				d.set('zscale')
			while True:
				answer = raw_input("Enter 0 to redo specidentify, or 1 to move onto next arc: ")
				if answer == '0' or answer == '1':
					break
				else:
					print "ERROR: you must enter either 0 or 1."
			if answer == '0': # delete wavelength solution and rectified arc before specidentify is redone 
				os.remove(idfilesci)
				os.remove(namearc)
				print "Wavelength solution and rectified arc deleted; re-running specidentify."
			elif answer == '1': # add rectified arc to relevant global dictionary, continue onto next angle (break out of while loop)
				dicts.wavearcs[angle] = namearc
				print 'output: '+namearc
				break
	print 'The specidentify wavelength solution filenames are:'
	print dicts.wavesols
	print 'The wavelength-calibrated arc image filenames are:'
	print dicts.wavearcs

# This function calls saltspec.specrectify to apply the wavelength solution to arc images.
def wavecalarc(arcimages=dicts.flatarcs,sols=dicts.wavesols,indiv=False):
	# This checks for the existence of input dictionaries by checking their length; if one is non-existent, quits function
	try:
		len(arcimages)
	except:
		print 'The inputted dictionary of flat-fielded arc images does not exist.'
		return
	try:
		len(sols)
	except:
		print 'The inputted dictionary of wavelength solutions does not exist.'
		return
	# This runs saltspec.specrectify on each flat-fielded arc image.
	angles = arcimages.keys() #all arcs should have been flat-fielded up to this point
	for angle in angles:
		# this stores the filenames of an arc image and its associated wavelength solution
		flatarc = arcimages.get(angle,'')
		if flatarc == '':
			print "WARNING: No flat-fielded arc image exists for this angle: "+angle
			print "Skipping rectification of arc image for this angle."
			continue
		suffixlist = flatarc.split('.')
		suffix = suffixlist[1][5:]
		hdu = pyfits.open(flatarc)
		hdr = hdu[0].header
		sol = hdr['RUSPECID'] # wavelength sol name
		hdu.close()
		# this stores the old headers of flatarc to copy over to the new image
		hdulist = pyfits.open(flatarc,mode='update')
		hdr0_old = hdulist[0].header.copy()
		hdr1_old = hdulist[1].header.copy()
		hdulist.close()
		# this sets the output filename for the 2-D wavelength-corrected arc image
		namearc = 'arc'+str(angle)+'wav'+suffix+'.fits'
		# this calls the run_specrectify function found further below
		pyrafCalls.run_specrectify(input=flatarc,output=namearc,idfile=sol,customRun=indiv)
		# this adds the 'RUPIPE' and 'RUIMGTYP' keywords to namearc's 0-header
		pipeHistory.updatePipeKeys(inputname=namearc,imagetype='arc',procChar='r')
		# this adds the filename of the new 2-D wavelength-corrected image to proper dictionary
		dicts.wavearcs[angle] = namearc
		# this notifies the user that an output was created by printing its name
		print 'output: '+namearc
	print 'The wavelength-calibrated arc image filenames are:'
	print dicts.wavearcs
		
		
# This function calls saltspec.specrectify to apply the wavelength solution to lacosmicx-corrected science images.
def wavecalsci(scienceimages=dicts.laxsciences,sols=dicts.wavesols,arcimages=dicts.flatarcs,indiv=False):
	# This checks for the existence of input dictionaries by checking their length; if one is non-existent, quits function
	try:
		len(scienceimages)
	except:
		print 'The inputted dictionary of flat-fielded science images does not exist.'
		return
	try:
		len(sols)
	except:
		print 'The inputted dictionary of wavelength solutions does not exist.'
		return
	angles = scienceimages.keys() #all science images should have been flat-fielded and run through lacosmicx, even with old flats if necessary
	for angle in angles:
		laxscience = scienceimages.get(angle,'')
		if laxscience == '':
			print "WARNING: No LAcosmicx-corrected science image exists for this angle: "+angle
			print "Skipping rectification of science image at this angle."
			continue
		suffixlist = laxscience.split('.')
		suffix = suffixlist[1][5:]
		arc = arcimages.get(angle,'')
		hdu = pyfits.open(arc)
		hdr = hdu[0].header
		sol = hdr['RUSPECID'] # wavelength sol name
		hdu.close()
		# this stores the old header of laxscience to copy over to the new image
		hdulist = pyfits.open(laxscience,mode='update')
		hdr0_old = hdulist[0].header.copy()
		hdr1_old = hdulist[1].header.copy()
		hdulist.close()
		# this sets the output filename for the 2-D wavelength-corrected science image
		namescience = 'sci'+str(angle)+'wav'+suffix+'.fits'
		# this calls the run_specrectify function found further below
		pyrafCalls.run_specrectify(input=laxscience,output=namescience,idfile=sol,customRun=indiv)
		# this displays the image in ds9 so the user can verify that the sky lines are straight
		print namescience+' has been opened in ds9 so you can verify that rectification was successful.'
		d = ds9.ds9(start=True)
		try:
			d.set('file '+namescience)
			d.set('zoom to fit')
			d.set('zscale')
		except: # Error most likely because user closed ds9 instance.
			d = ds9.ds9(start=True)
			d.set('file '+namescience)
			d.set('zoom to fit')
			d.set('zscale')
		# this adds the 'RUPIPE' and 'RUIMGTYP' keywords to namescience's 0-header
		pipeHistory.updatePipeKeys(inputname=namescience,imagetype='science',procChar='r')
		# this adds the filename of the new 2-D wavelength-corrected image to proper dictionary
		dicts.wavesciences[angle] = namescience
		# this notifies the user that an output was created by printing its name
		print 'output: '+namescience	
	print 'The wavelength-calibrated science image filenames are:'
	print dicts.wavesciences

# This function calls saltspec.specrectify to apply the wavelength solution to standard star images.
def wavecalstd(standardimages=dicts.flatstandards,sols=dicts.wavesols,arcimages=dicts.flatarcs,indiv=False):
	# This checks for the existence of input dictionaries by checking their length; if one is non-existent, quits function
	try:
		len(standardimages)
	except:
		print 'The inputted dictionary of flat-fielded standard star images does not exist.'
		return
	try:
		len(sols)
	except:
		print 'The inputted dictionary of wavelength solutions does not exist.'
		return
	angles = standardimages.keys() #all standards, if existent, should have been flat-fielded up to this point, even with old flats if necessary
	for angle in angles:
		flatstandard = standardimages.get(angle,'')
		if flatstandard == '':
			print "WARNING: No flat-fielded standard star image exists for this angle: "+angle
			print "Skipping rectification of standard star image at this angle."
			continue
		suffixlist = flatstandard.split('.')
		suffix = suffixlist[1][5:]
		arc = arcimages.get(angle,'')
		hdu = pyfits.open(arc)
		hdr = hdu[0].header
		sol = hdr['RUSPECID'] # wavelength sol name
		hdu.close()
		# this stores the old header of flatstandard to copy over to the new image
		hdulist = pyfits.open(flatstandard,mode='update')
		hdr0_old = hdulist[0].header.copy()
		hdr1_old = hdulist[1].header.copy()
		hdulist.close()
		# this sets the output name for the 2-D wavelength-corrected standard star image
		namestandard = 'std'+str(angle)+'wav'+suffix+'.fits'
		# this calls the run_specrectify function found further below
		pyrafCalls.run_specrectify(input=flatstandard,output=namestandard,idfile=sol,customRun=indiv)
		# this adds the 'RUPIPE' and 'RUIMGTYP' keywords to namestandard's 0-header
		pipeHistory.updatePipeKeys(inputname=namestandard,imagetype='standard',procChar='r')
		# this adds the name of the new 2-D wavelength-corrected image to proper dictionary
		dicts.wavestandards[angle] = namestandard
		# this notifies the user that an output was created by printing its name
		print 'output: '+namestandard
	print 'The wavelength-calibrated standard star image filenames are:'
	print dicts.wavestandards

# This function subtracts the background (sky lines essentially) from the two-dimensional science images.
def subtractbackground(scienceimages=dicts.wavesciences,indiv=False):
	# This checks if the inputted dictionary of images exists; if not, it quits the function
	try:
		len(scienceimages)
	except:
		print 'The inputted dictionary of wavelength-corrected 2D science images does not exist.'
		return
	# this sets the DISPAXIS=1 in the header and calls the run_background function further below
	angles = scienceimages.keys()
	for angle in angles:
		wavescience = scienceimages[angle]
		suffixlist = wavescience.split('.')
		suffix = suffixlist[1][5:]
		hdulist = pyfits.open(wavescience,mode='update')
		hdr = hdulist[1].header
		hdr.update('DISPAXIS',1)
		hdulist.flush()
		hdulist.close()
		# output filenames
		name = 'sci'+str(angle)+'bkg'+suffix+'.fits'
		auxname = 'sci'+str(angle)+'bkg'+suffix+'AUX'+'.fits'
		# this displays the image in ds9 to ask user if background should be run interactively (for faint spectra)
		print "You will now use the PyRAF task BACKGROUND to background-subtract the image: "+wavescience
		print "For your convenience, the image has also been opened in ds9 so you can find:"
		print '1. The central row of the supernova spectrum. (Consult your acquisition images if needed.)'
		print '2. Whether the supernova spectrum is bright or faint (barely distinguishable from the background).'
		print 'It is recommended to run BACKGROUND interactively if the spectrum is faint (enter 1).'
		# opens current image in ds9
		d = ds9.ds9(start=True)
		try:
			d.set('file '+wavescience)
			d.set('zoom to fit')
			d.set('zscale')
		except: # Error most likely because user closed ds9 instance.
			d = ds9.ds9(start=True)
			d.set('file '+wavescience)
			d.set('zoom to fit')
			d.set('zscale')
		# this asks the user if the spectral line looks very faint (look in ds9)
		# if the user says the line is very faint, then apall is run with a different set of parameters
		while True:
			answer = raw_input("Please enter 0 if the spectrum looks bright, or 1 if it looks faint: ")
			if answer == '0' or answer == '1':
				break
			else:
				print "Invalid input. You must enter either 0 or 1."
		# this calls the run_background function
		print "Running background subtraction on: "+wavescience
		pyrafCalls.run_background(twodimage=wavescience,newimage=auxname,faint=answer,customRun=indiv)
		# this copies the old images to preserve file structure, and transfers the extracted data+header
		pyrafCalls.run_imcopy(inputname=wavescience,outputname=name)
		hdulist = pyfits.open(auxname)
		datnew = hdulist[0].data.copy()
		hdr = hdulist[0].header
		hdr.update('EXTNAME','SCI')
		hdr.update('EXTVER',1)
		hdrnew = hdr.copy()
		hdulist.close()
		hdunew = pyfits.open(name,mode='update')
		hdunew[1].data = datnew
		hdunew[1].header = hdrnew
		hdr = hdunew[0].header
		hdr['RUFAINT'] = answer # necessary for automatic apsum extraction of non-bkg-subtracted sci images in extractSciSigma()
		hdunew.flush()
		hdunew.close()
		pipeHistory.updatePipeKeys(inputname=name,imagetype='science',procChar='b')
		# Add background-subtracted science image filename to dictionary backgroundsciences indexed by gr-angle.
		dicts.backgroundsciences[angle] = name
		# this deletes the aux image file
		pyrafCalls.run_imdel(auxname)
		# notify user of successful background-subtraction
		print 'output: '+name
	print 'The background-subtracted science image filenames are:'
	print dicts.backgroundsciences
		

# This function runs apall on the wavelength-corrected, background-subtracted science images to extract the SN spectra.
def extractsciences(scienceimages=dicts.backgroundsciences,wavescispectra=dicts.wavesciences,indiv=False):
	# This checks if the inputted dictionary of images exists; if not, it quits the function
	try:
		len(scienceimages)
	except:
		print 'The inputted dictionary of wavelength-corrected 2D images does not exist.'
		return	
	# This adds the DISPAXIS keyword to each image's header with value 1 and then runs apall.
	angles = scienceimages.keys()
	for angle in angles:
		bkgscience = scienceimages[angle]
		suffixlist = bkgscience.split('.')
		suffix = suffixlist[1][5:]
		hdulist = pyfits.open(bkgscience,mode='update')
		hdr = hdulist[1].header
		hdr.update('DISPAXIS',1)
		hdr0 = hdulist[0].header
		airmass = hdr0['AIRMASS']
		exptime = hdr0['EXPTIME']
		object = hdr0['OBJECT']
		thegain = hdr0.get('GAIN',2.443)
		hdr.update('DISPAXIS',1)
		hdr.update('AIRMASS',airmass)
		hdr.update('EXPTIME',exptime)
		hdr.update('OBJECT',object)
		hdulist.flush()
		hdr0_old = hdulist[0].header.copy()
		hdr1_old = hdulist[1].header.copy()
		hdulist.close()
		# this displays the image in ds9 for the user's convenience
		print "You will now use the PyRAF task APALL to extract the spectrum in the image: "+bkgscience
		print "For your convenience, the image has also been opened in ds9 so you can find:"
		print '1. The central row of the supernova spectrum. (Consult your acquisition images if needed.)'
		print '2. Whether the supernova spectrum is bright or faint (barely distinguishable from the background).'
		# opens current image in ds9
		d = ds9.ds9(start=True)
		hdulist = pyfits.open(bkgscience)
		try:
			d.set('file '+bkgscience)
			d.set('zoom to fit')
			d.set('zscale')
		except: # Error most likely because user closed ds9 instance.
			d = ds9.ds9(start=True)
			d.set('file '+bkgscience)
			d.set('zoom to fit')
			d.set('zscale')
		# this asks the user if the spectral line looks very faint (look in ds9)
		# if the user says the line is very faint, then apall is run with a different set of parameters
		print "1: Spectrum looks BRIGHT. (No local background subtraction.)"
		print "2: Spectrum looks BRIGHT but has extended galaxy emission (local background subtraction)."
		print "3: Spectrum looks FAINT. (No local background subtraction; good for host galaxy.)"
		print "4: Spectrum looks FAINT but has extended galaxy emission (local background subtraction)."
		print "For 4, a pre-sky-subtracted level will be re-added to the background-subtracted image."
		while True:
			answer = raw_input("Please enter 1, 2, 3 or 4 based on the above options: ")
			if answer == '1' or answer == '2' or answer == '3' or answer == '4':
				break
			else:
				print "Invalid input. You must enter either 1, 2, 3 or 4."
		# declares output filenames and calls run_apall function
		name = 'sci'+str(angle)+'ext'+suffix+'.fits'
		auxname = 'sci'+str(angle)+'ext'+suffix+'AUX'+'.fits'
		if answer == '4':
			# Ask user for an estimate of the median local background around the spectrum in the non-background-subtracted image.
			wavescience = wavescispectra[angle] # non-background-subtracted science image
			d = ds9.ds9(start=True)
			hdulist = pyfits.open(wavescience)
			try:
				d.set('file '+wavescience)
				d.set('zoom to fit')
				d.set('zscale')
			except: # Error most likely because user closed ds9 instance.
				d = ds9.ds9(start=True)
				d.set('file '+wavescience)
				d.set('zoom to fit')
				d.set('zscale')
			print "You chose option 4. The non-background-subtracted image has been opened in ds9."
			while True:
				pssl = raw_input("Please enter an estimate of the median local background level around the spectrum: ")
				if pssl != '': # need actual error condition using pssl
					break
				else:
					print "Invalid input. Use ds9 to estimate the median local background level around the spectrum."
			# add pssl to background-subtracted image (overwrite it)
			pyrafCalls.run_imarithGeneral(op1=bkgscience+'[1]',op2=pssl,oper='+',outname='tempapall.fits',customRun=False)
			hdu = pyfits.open('tempapall.fits')
			datnew = hdu[0].data.copy()
			hdu.close()
			hdu = pyfits.open(bkgscience,mode='update')
			hdu[1].data = datnew
			hdu.flush()
			hdu.close()
			os.remove('tempapall.fits')
		pyrafCalls.run_apall(twodimage=bkgscience,spectrum=auxname,saltgain=thegain,faint=answer,customRun=indiv)
		# sets output filenames and calls run_apall function
		pyrafCalls.run_imcopy(inputname=bkgscience,outputname=name)
		hdulist = pyfits.open(auxname)
		datnew = hdulist[0].data.copy()
		hdr = hdulist[0].header
		hdr.update('EXTNAME','SCI')
		hdr.update('EXTVER',1)
		hdrnew = hdr.copy()
		hdulist.close()
		hdunew = pyfits.open(name,mode='update')
		hdunew[1].data = datnew
		hdunew[1].header = hdrnew
		hdunew.flush()
		if answer == '0' or answer == '2': # for spectra without background, add a 4th band with temporary data which will be replaced in extractSciSigma(...)
			dat4 = np.array([[datnew[2][0]]])
			datnew = np.concatenate((datnew,dat4))
			hdunew[1].data = datnew
			hdunew.flush()
		hdunew.close()
		pipeHistory.updatePipeKeys(inputname=name,imagetype='science',procChar='e')
		# Add extracted science image filename to dictionary extractedsciences indexed by gr-angle.
		dicts.extractedsciences[angle] = name
		# this deletes aux files
		pyrafCalls.run_imdel(auxname)
		# notify user of successful extraction
		print 'output: '+name
	print 'The extracted science image filenames are:'
	print dicts.extractedsciences
	

# This function runs apall on the wavelength-corrected standard star images to extract the standard star spectra.
def extractstandards(standardimages=dicts.wavestandards,indiv=False):
	# This checks if the inputted dictionary of images exists; if not, it quits the function
	try:
		len(standardimages)
	except:
		print 'The inputted dictionary of wavelength-corrected 2D images does not exist.'
		return
	# This adds the DISPAXIS keyword to each image's header with value 1 and then runs apall.
	angles = standardimages.keys()
	for angle in angles:
		wavestandard = standardimages[angle]
		suffixlist = wavestandard.split('.')
		suffix = suffixlist[1][5:]
		hdulist = pyfits.open(wavestandard,mode='update')
		hdr = hdulist[1].header
		hdr.update('DISPAXIS',1)
		hdr0 = hdulist[0].header
		airmass = hdr0['AIRMASS']
		exptime = hdr0['EXPTIME']
		object = hdr0['OBJECT']
		thegain = hdr0.get('GAIN',2.443)
		hdr.update('DISPAXIS',1)
		hdr.update('AIRMASS',airmass)
		hdr.update('EXPTIME',exptime)
		hdr.update('OBJECT',object)
		hdulist.flush()
		hdr0_old = hdulist[0].header.copy()
		hdr1_old = hdulist[1].header.copy()
		hdulist.close()
		# setting answer =2 (bright spectrum with local bkg subtraction) since standards will be extracted non-interactively
		answer = '2'
		# output filenames
		name = 'std'+str(angle)+'ext'+suffix+'.fits'
		auxname = 'std'+str(angle)+'ext'+suffix+'AUX'+'.fits'
		pyrafCalls.run_apall(twodimage=wavestandard,spectrum=auxname,saltgain=thegain,faint=answer,customRun=indiv)
		# this copies the old images to preserve file structure, and transfers the extracted data+header
		pyrafCalls.run_imcopy(inputname=wavestandard,outputname=name)
		hdulist = pyfits.open(auxname)
		datnew = hdulist[0].data.copy()
		hdr = hdulist[0].header
		hdr.update('EXTNAME','SCI')
		hdr.update('EXTVER',1)
		hdrnew = hdr.copy()
		hdulist.close()
		hdunew = pyfits.open(name,mode='update')
		hdunew[1].data = datnew
		hdunew[1].header = hdrnew
		hdunew.flush()
		hdunew.close()
		pipeHistory.updatePipeKeys(inputname=name,imagetype='standard',procChar='e')
		# Add extracted standard star image filename to dictionary extractedstandards indexed by gr-angle.
		dicts.extractedstandards[angle] = name
		# this deletes aux files
		pyrafCalls.run_imdel(auxname)
		# notify user of successful extraction
		print 'output: '+name
	print 'The extracted standard star image filenames are:'
	print dicts.extractedstandards


# This function runs apsum on the wavelength-corrected arc images to extract the arc spectra.
# will be run twice for each arc -- once for science image, and once for standard star image
def extractarcs(arcimages=dicts.wavearcs,refsciences=dicts.backgroundsciences,refstandards=dicts.wavestandards,indiv=False):
	# This checks if the inputted dictionary of images exists; if not, quits the function
	try:
		len(arcimages)
	except:
		print 'The inputted dictionary of images does not exist.'
		return
	# This runs apsum twice for each arc image (for science-arcs and standard-arcs)
	angles = arcimages.keys()
	for angle in angles:
		wavearc = arcimages[angle]
		suffixlist = wavearc.split('.')
		suffix = suffixlist[1][5:]
		# this sets 'DISPAXIS'=1 and stores the old header so it can be copied over to the new images
		hdulist = pyfits.open(wavearc,mode='update')
		hdr = hdulist[1].header
		hdr0 = hdulist[0].header
		airmass = hdr0['AIRMASS']
		exptime = hdr0['EXPTIME']
		object = hdr0['OBJECT']
		hdr.update('DISPAXIS',1)
		hdr.update('AIRMASS',airmass)
		hdr.update('EXPTIME',exptime)
		hdr.update('OBJECT',object)
		hdulist.flush()
		hdr0_old = hdulist[0].header.copy()
		hdr1_old = hdulist[1].header.copy()
		hdulist.close()
		
		# this sets the output and reference image filename for the science-reference arc
		namesci = 'arc'+str(angle)+'ext'+suffix+'_sci.fits'
		auxnamesci = 'arc'+str(angle)+'ext'+suffix+'_sciAUX.fits'
		refsci = refsciences.get(angle,'')
		if refsci == '':
			print 'No reference science image found for angle: '+angle
			print 'Skipping extraction of sci-arc for this angle.'
			continue
		# this calls the run_apsum function found further below for the SCI-extracted arc
		pyrafCalls.run_apsum(twodimage=wavearc,refimage=refsci+'[1]',spectrum=auxnamesci,customRun=indiv)
		# this copies the old images to preserve file structure, and transfers the extracted data+header
		pyrafCalls.run_imcopy(inputname=wavearc,outputname=namesci)
		hdulist = pyfits.open(auxnamesci)
		datnew = hdulist[0].data.copy()
		hdr = hdulist[0].header
		hdr.update('EXTNAME','SCI')
		hdr.update('EXTVER',1)
		hdrnew = hdr.copy()
		hdulist.close()
		hdunew = pyfits.open(namesci,mode='update')
		hdunew[1].data = datnew
		hdunew[1].header = hdrnew # this is the extracted spectrum's header (with spectral WCS info, etc.)
		hdunew.flush()
		hdunew.close()
		pipeHistory.updatePipeKeys(inputname=namesci,imagetype='arc_sci',procChar='e')
		# Add extracted arc image to dictionary extractedarcs_sci indexed by gr-angle.
		dicts.extractedarcs_sci[angle] = namesci
		# this deletes aux files
		pyrafCalls.run_imdel(auxnamesci)
		# notify user of successful extraction
		print 'output: '+namesci
		
		if refstandards != {}: # only extract arcs for standards if we have standards
			# this sets the output and reference image names for the standard star-reference arc
			namestd = 'arc'+str(angle)+'ext'+suffix+'_std.fits'
			auxnamestd = 'arc'+str(angle)+'ext'+suffix+'_stdAUX.fits'
			try: 
				refstd = refstandards[angle]+'[1]'
			except:
				print "No standard star for angle: "+angle
				print "Cannot extract arc for this non-existent standard star."
				continue
			# this calls the run_apsum function found further below for the STD-extracted arc
			pyrafCalls.run_apsum(twodimage=wavearc,refimage=refstd,spectrum=auxnamestd,customRun=indiv)
			# this copies the old images to preserve file structure, and transfers the extracted data+header
			pyrafCalls.run_imcopy(inputname=wavearc,outputname=namestd)
			hdulist = pyfits.open(auxnamestd)
			datnew = hdulist[0].data.copy()
			hdr = hdulist[0].header
			hdr.update('EXTNAME','SCI')
			hdr.update('EXTVER',1)
			hdrnew = hdr.copy()
			hdulist.close()
			hdunew = pyfits.open(namestd,mode='update')
			hdunew[1].data = datnew
			hdunew[1].header = hdrnew
			hdunew.flush()
			hdunew.close()
			pipeHistory.updatePipeKeys(inputname=namestd,imagetype='arc_std',procChar='e')
			# Add STD-extracted arc image to dictionary extractedarcs_std indexed by gr-angle.
			dicts.extractedarcs_std[angle] = namestd
			# this deletes aux files
			pyrafCalls.run_imdel(auxnamestd)
			# notify user of successful extraction
			print 'output: '+namestd
	
	print 'The extracted arc (for science images) image filenames are:'
	print dicts.extractedarcs_sci
	if dicts.extractedarcs_std != {}: # only print extractedarcs_std dictionary if we have standards
		print 'The extracted arc (for standard star images) image filenames are:'
		print dicts.extractedarcs_std

		
# This function runs apsum on the wavelength-corrected science images.
# Its purpose is to acquire the sky spectrum (band 3) and sigma spectrum (band 4)
# from the non-2D-background-subtracted image since that should better approximate the background and error.
# The first background-extracted science image will be fed in as the reference image for automatic aperture finding.
def extractSciSigma(scienceimages=dicts.wavesciences,refsciences=dicts.backgroundsciences,mainsciences=dicts.extractedsciences,indiv=False):
	# This checks if the inputted dictionary of images exists; if not, quits the function
	try:
		len(scienceimages),len(refsciences),len(mainsciences)
	except:
		print 'Cannot proceed: an inputted dictionary of images does not exist.'
		return
	# This runs apsum on each wavescience image (with corresponding backgroundscience as reference)
	angles = scienceimages.keys()
	for angle in angles:
		wavesci = scienceimages[angle]
		# this checks for the existence of the main extracted science spectra (2D-background-subtracted extractions)
		mainsci = mainsciences.get(angle,'')
		if mainsci == '': # no need to extract original sky and sigma spectra if there's nowhere to insert them
			print "There does not exist a 2D-background-subtracted science extraction for angle: "+angle
			print "Skipping extraction and insertion of original sky and sigma spectra for image: "+wavesci
			continue
		# this checks the bkgsci 0-header for 'RUFAINT' keyword (1 => yes, faint spectrum; default = 0)
		bkgsci = refsciences.get(angle,'') # there should exist a background science for each wavescience
		# but if no bkgsci exists, skip this apsum extraction and just leave extsci's sky and sigma spectrum as they are
		if bkgsci == '':
			print "No background science image exists for this angle: "+angle
			print "Cannot acquire, with apsum, a sky spectrum and sigma spectrum for: "+wavesci
			print "Retaining the current extracted science spectrum's sky and sigma spectra."
			continue
		hdubkg =  pyfits.open(bkgsci) # else, find the 'RUFAINT' keyword
		hdrbkg = hdubkg[0].header
		faintSpectrum = hdrbkg.get('RUFAINT','0') # default value is '0', i.e., not faint
		hdubkg.close()
		# this gets the relevant information about the wavesci image, and updates the 0-header and 1-header
		suffixlist = wavesci.split('.')
		suffix = suffixlist[1][5:]
		# this sets 'DISPAXIS'=1 in wavesci's 1-header
		hdulist = pyfits.open(wavesci,mode='update')
		hdr = hdulist[1].header
		hdr0 = hdulist[0].header
		airmass = hdr0['AIRMASS']
		exptime = hdr0['EXPTIME']
		object = hdr0['OBJECT']
		hdr.update('DISPAXIS',1)
		hdr.update('AIRMASS',airmass)
		hdr.update('EXPTIME',exptime)
		hdr.update('OBJECT',object)
		hdr0.update('RUFAINT',faintSpectrum)
		hdulist.flush()
		hdr0_old = hdulist[0].header.copy()
		hdr1_old = hdulist[1].header.copy()
		hdulist.close()
		# this sets the output filename for the apsum-extracted science spectrum
		namesci = 'sci'+str(angle)+'aps'+suffix+'.fits'
		auxnamesci = 'sci'+str(angle)+'aps'+suffix+'AUX.fits' # output name for auxiliary science (with data)
		# this calls the run_apsumSci function found further below
		pyrafCalls.run_apsumSci(inputimage=wavesci+'[1]',refimage=bkgsci+'[1]',spectrum=auxnamesci,faint=faintSpectrum,customRun=indiv)
		# this copies the old images to preserve file structure, and transfers the extracted data+header
		pyrafCalls.run_imcopy(inputname=wavesci,outputname=namesci)
		hdulist = pyfits.open(auxnamesci)
		datnew = hdulist[0].data.copy()
		hdr = hdulist[0].header
		hdr.update('EXTNAME','SCI')
		hdr.update('EXTVER',1)
		hdrnew = hdr.copy()
		hdulist.close()
		hdunew = pyfits.open(namesci,mode='update')
		hdunew[1].data = datnew
		hdunew[1].header = hdrnew # this is the extracted spectrum's header (with spectral WCS info, etc.)
		hdunew.flush()
		hdunew.close()
		pipeHistory.updatePipeKeys(inputname=namesci,imagetype='science',procChar='j')
		# Add apsum-extracted science spectrum to dictionary apsumsciences indexed by gr-angle.
		dicts.apsumsciences[angle] = namesci
		# this deletes aux files
		pyrafCalls.run_imdel(auxnamesci)
		# insert namesci's sky and sigma extraction into the main 2D-background-extracted mainsci
		hduaps = pyfits.open(namesci,mode='update')
		dataps = hduaps[1].data.copy()
		hdra = hduaps[0].header
		hdra['RUAPSUM'] = mainsci # new keyword 'RUAPSUM' used to identify if mainsci's sky and sigma spectra were replaced with namesci's
		hduaps.flush()
		hduaps.close()
		newband3 = dataps[2] # 0-based NumPy array of 4 arrays (called 'spectral bands' in IRAF)
		newband4 = dataps[3] 
		hdumain = pyfits.open(mainsci,mode='update')
		datmain = hdumain[1].data.copy()
		datmain[2] = newband3 # this assumes that mainsci has well-defined sky and sigma numpy arrays, else index error will occur
		datmain[3] = newband4
		hdumain[1].data = datmain
		hdumain.flush()
		hdrm = hdumain[0].header
		hdrm['RUAPSUM'] = namesci # new keyword 'RUAPSUM' used to identify if a mainsci has had its sky and sigma data updated
		hdumain.close()		
		# notify user of successful sky and sigma spectra extraction and insertion
		print 'Successfully transferred sky and sigma spectrum of '+namesci+' to: '+mainsci
	print 'The apsum-extracted science spectrum filenames are:'
	print dicts.apsumsciences
		

# This function runs identify on the wavelength-corrected, SCI-extracted arc images to reduce errors further.
def identifyarcs(arcsci=dicts.extractedarcs_sci,scienceimages=dicts.extractedsciences,indiv=False):
	# This prints an error and quits the function if any of the inputted dictionaries don't exist
	try:
		len(arcsci),len(scienceimages)
	except:
		print 'An inputted dictionary does not exist: cannot proceed.'
		return
	# this prepares for and runs identify on each extracted arc (associated with a science image)
	angles = arcsci.keys()
	for angle in angles:
		# this stores the filenames of the extracted arc and its associated extracted science
		extarc = arcsci[angle]
		extsci = scienceimages[angle] # should exist if arcsci[angle] exists
		# this gets the arc lamp type so the correct linelist can be used by identify
		hdulist = pyfits.open(extarc,mode='update')
		hdr = hdulist[0].header
		lamp = hdr['LAMPID']
		hdulist.close()
		# this sets the arguments for run_identify function below
		linelistpath = params.lineListPath
		# picking the linelist file based on header LAMPID reference variable lamp
		if lamp == 'Th Ar':
			print 'the lamp is '+lamp+' for angle '+angle
			linelistpath = params.ThAr
		elif lamp == 'Xe':
			print 'the lamp is '+lamp+' for angle '+angle
			linelistpath = params.Xe
		elif lamp == 'Ne': # might instead need to use NeAr.salt
			print 'the lamp is '+lamp+' for angle '+angle
			linelistpath = params.Ne
		elif lamp == 'Cu Ar':
			print 'the lamp is '+lamp+' for angle '+angle
			linelistpath = params.CuAr
		elif lamp == 'Ar':
			print 'the lamp is '+lamp+' for angle '+angle
			linelistpath = params.Ar
		elif lamp == 'Hg Ar':
			print 'the lamp is '+lamp+' for angle '+angle
			linelistpath = params.HgAr
		else:
			print 'Could not find the proper linelist for '+lamp+' lamp.'
			continue
		# this calls the run_identify function found further below on extarc with data extension [1] ([SCI])
		extarc1 = extarc+'[1]'
		pyrafCalls.run_identify(arcsci=extarc1,lamplines=linelistpath,customRun=indiv)
		# this adds the 'REFSPEC1' keyword to the associated extracted science image for dispcor later
		hdulist = pyfits.open(extsci,mode='update')
		hdr = hdulist[0].header
		hdr.update('REFSPEC1',extarc)
		hdulist.flush()
		hdulist.close()
		# this adds the 'RUPIPE' and 'RUIMGTYP' keywords to the identified arc
		pipeHistory.updatePipeKeys(inputname=extarc,imagetype='arc_sci',procChar='i')
		# notifies user of successful identification 
		print 'emission lines in the file '+extarc+' have been identified'
		print 'and it is now linked to its science image: '+extsci
		

# This function runs reidentify on the remaining STD-extracted arcs.
def reidentifyarcs(arcstd=dicts.extractedarcs_std,arcsci=dicts.extractedarcs_sci,standardimages=dicts.extractedstandards,indiv=False):
	# This prints an error and quits the function if any of the inputted dictionaries don't exist
	try:
		len(arcstd),len(arcsci),len(standardimages)
	except:
		print 'An inputted dictionary does not exist: cannot proceed.'
		return
	# this prepares for and runs reidentify on each extracted arc (associated with a standard star image)
	angles = arcstd.keys()
	for angle in angles:
		# this stores the filenames of the extracted arc, reference sci-extracted arc, and associated extracted std star
		extarcstd = arcstd[angle] # current arc to be identified
		extarcsci = arcsci.get(angle,'') # reference arc already identified (with identify)
		if extarcsci == '':
			print 'No (reference) extracted arc (for science image) exists for angle: '+angle
			print 'Skipping reidentify for this extracted arc (for standard star).'
			continue
		extstd = standardimages[angle] # associated standard star image
		# this gets the arc lamp type so the correct linelist can be used by identify
		hdulist = pyfits.open(extarcstd,mode='update')
		hdr = hdulist[0].header
		lamp = hdr['LAMPID']
		hdulist.close()
		# this sets the arguments for run_reidentify function below
		linelistpath = params.lineListPath
		# picking the linelist file based on header LAMPID reference variable lamp
		if lamp == 'Th Ar':
			print 'the lamp is '+lamp+' for angle '+angle
			linelistpath = params.ThAr
		elif lamp == 'Xe':
			print 'the lamp is '+lamp+' for angle '+angle
			linelistpath = params.Xe
		elif lamp == 'Ne': # might instead need to use NeAr.salt
			print 'the lamp is '+lamp+' for angle '+angle
			linelistpath = params.Ne
		elif lamp == 'Cu Ar':
			print 'the lamp is '+lamp+' for angle '+angle
			linelistpath = params.CuAr
		elif lamp == 'Ar':
			print 'the lamp is '+lamp+' for angle '+angle
			linelistpath = params.Ar
		elif lamp == 'Hg Ar':
			print 'the lamp is '+lamp+' for angle '+angle
			linelistpath = params.HgAr
		else:
			print 'Could not find the proper linelist for '+lamp+' lamp.'
			continue
		# this calls the run_reidentify function found further below on extarcstd with data extension [1] ([SCI])
		extarcstd1 = extarcstd+'[1]'
		extarcsci1 = extarcsci+'[1]'
		pyrafCalls.run_reidentify(input=extarcstd1,ref=extarcsci1,lamplines=linelistpath,customRun=indiv)
		# this adds the 'REFSPEC1' keyword to the associated extracted standard star image for dispcor later
		hdulist = pyfits.open(extstd,mode='update')
		hdr = hdulist[0].header
		hdr.update('REFSPEC1',extarcstd)
		hdulist.flush()
		hdulist.close()
		# this adds the 'RUPIPE' and 'RUIMGTYP' keywords to the identified arc
		pipeHistory.updatePipeKeys(inputname=extarcstd,imagetype='arc_std',procChar='i')
		# notifies user of successful identification 
		print 'emission lines in the file '+extarcstd+' have been identified'
		print 'and it is now linked to its standard star image: '+extstd
		

# This function runs dispcor on the extracted science images to apply the lower-error wavelength solution.
# REFSPEC1 keyword was already added to images by the functions, identifyarcs and reidentifyarcs
def dispcorsci(scienceimages=dicts.extractedsciences,indiv=False):
	# This checks for the existence of input dictionaries by checking their length; if one is non-existent, quits function
	try:
		len(scienceimages)
	except:
		print 'The inputted dictionary does not exist: cannot proceed.'
		return
	# This runs onedspec.dispcor on each extracted science image.
	angles = scienceimages.keys() # each science image should have an extracted and identified sci-arc
	for angle in angles:
		# this stores the filenames of a science image
		extscience = scienceimages[angle]
		suffixlist = extscience.split('.')
		suffix = suffixlist[1][5:]
		# this stores the old header of extscience to copy over to the new image
		hdulist = pyfits.open(extscience,mode='update')
		hdr0_old = hdulist[0].header.copy()
		hdr1_old = hdulist[1].header.copy()
		hdulist.close()
		# this sets the output filename for the 1-D wavelength-corrected science image
		namescience = 'sci'+str(angle)+'dsp'+suffix+'.fits'
		auxnamescience = 'sci'+str(angle)+'dsp'+suffix+'AUX.fits'
		# this calls the run_dispcor function found further below on extscience with data extension [1] ([SCI])
		extscience1 = extscience+'[1]'
		pyrafCalls.run_dispcor(original=extscience1,corrected=auxnamescience,customRun=indiv)
		# this copies the old images to preserve file structure, and transfers the extracted data+header
		pyrafCalls.run_imcopy(inputname=extscience,outputname=namescience)
		hdulist = pyfits.open(auxnamescience)
		datnew = hdulist[0].data.copy()
		hdr = hdulist[0].header
		hdr.update('EXTNAME','SCI')
		hdr.update('EXTVER',1)
		hdrnew = hdr.copy()
		hdulist.close()
		hdunew = pyfits.open(namescience,mode='update')
		hdunew[1].data = datnew
		hdunew[1].header = hdrnew
		hdunew.flush()
		hdunew.close()
		pipeHistory.updatePipeKeys(inputname=namescience,imagetype='science',procChar='d')
		# this adds the filename of the new 1-D wavelength-corrected image to proper dictionary
		dicts.dispsciences[angle] = namescience
		# this deletes aux files
		pyrafCalls.run_imdel(auxnamescience)
		# this notifies the user that an output was created by printing its name
		print 'output: '+namescience
	print 'The dispersion-corrected science image filenames are:'
	print dicts.dispsciences


# This function runs dispcor on the extracted standard star images to apply the lower-error wavelength solution.
# REFSPEC1 keyword was already added to images by the functions, identifyarcs and reidentifyarcs
def dispcorstd(standardimages=dicts.extractedstandards,indiv=False):
	# This checks for the existence of input dictionaries by checking their length; if one is non-existent, quits function
	try:
		len(standardimages)
	except:
		print 'The inputted dictionary does not exist: cannot proceed.'
		return
	# This runs onedspec.dispcor on each extracted standard star image.
	angles = standardimages.keys() # each standard star should have an extracted and identified std-arc
	for angle in angles: # won't run if there are no standards for if not, then there are no angles either
		# this stores the filenames of a science image
		extstandard = standardimages[angle]
		suffixlist = extstandard.split('.')
		suffix = suffixlist[1][5:]
		# this stores the old header of extstandard to copy over to the new image
		hdulist = pyfits.open(extstandard,mode='update')
		hdr0_old = hdulist[0].header.copy()
		hdr1_old = hdulist[1].header.copy()
		hdulist.close()
		# this sets the output filename for the 2-D wavelength-corrected standard star image
		namestandard = 'std'+str(angle)+'dsp'+suffix+'.fits'
		auxnamestandard = 'std'+str(angle)+'dsp'+suffix+'AUX.fits'
		# this calls the run_dispcor function found further below on extstandard with data extension [1] ([SCI])
		extstandard1 = extstandard+'[1]'
		pyrafCalls.run_dispcor(original=extstandard1,corrected=auxnamestandard,customRun=indiv)
		# this copies the old images to preserve file structure, and transfers the extracted data+header
		pyrafCalls.run_imcopy(inputname=extstandard,outputname=namestandard)
		hdulist = pyfits.open(auxnamestandard)
		datnew = hdulist[0].data.copy()
		hdr = hdulist[0].header
		hdr.update('EXTNAME','SCI')
		hdr.update('EXTVER',1)
		hdrnew = hdr.copy()
		hdulist.close()
		hdunew = pyfits.open(namestandard,mode='update')
		hdunew[1].data = datnew
		hdunew[1].header = hdrnew
		hdunew.flush()
		hdunew.close()
		pipeHistory.updatePipeKeys(inputname=namestandard,imagetype='standard',procChar='d')
		# this adds the filename of the new 2-D wavelength-corrected image to proper dictionary
		dicts.dispstandards[angle] = namestandard
		# this deletes aux files
		pyrafCalls.run_imdel(auxnamestandard)
		# this notifies the user that an output was created by printing its name
		print 'output: '+namestandard
	print 'The dispersion-corrected standard star image filenames are:'
	print dicts.dispstandards


# This function creates the bad pixel ("wavelength") masks (chip gaps) for all dispersion-corrected science spectra.
# Bad wavelength masks created from dispersion-corrected science spectra should remain applicable to flux-calibrated science spectra.
# This function also creates bad wavelength masks for dispersion-corrected standard star spectra in case they come in handy for flux calibration.
def maskspectra(dispscispectra=dicts.dispsciences,dispstdspectra=dicts.dispstandards,indiv=False):
	# This prints an error message and quits function if dispersion-corrected science spectra do not exist.
	try:
		len(dispscispectra)
	except:
		print 'The inputted dictionary of dispersion-corrected science spectra does not exist: cannot proceed.'
		return
	# This runs imexpr on each dispersion-corrected science and standard star spectrum for all gr-angles. 
	angles = dispscispectra.keys() # there may not be standard stars, but there must be science spectra
	for angle in angles:
		dispscience = dispscispectra.get(angle,'')
		dispstandard = dispstdspectra.get(angle,'')
		if dispstandard == '': # no dispstandard => no standard star for this angle
			print "There is no dispersion-corrected standard star spectrum for angle: "+angle
			print "Skipping creation of bad wavelength mask for a standard star spectrum at this angle."
		if dispscience != '': # create BWM for science spectrum if science spectrum exists
			suffixlist = dispscience.split('.')
			suffix = suffixlist[1][5:]
			namebwmsci = 'sci'+str(angle)+'bwm'+suffix+'.fits'
			imutil.imexpr.unlearn()
			imutil.imexpr(expr='a==0',output=namebwmsci,a=dispscience+'[1]')
			# this adds the 'RUIMGTYP' and 'RUPIPE' keywords to the bad wavelength mask
			hdubwm = pyfits.open(namebwmsci,mode='update')
			# this manually sets all pixels to 1 within the chip gap regions (based on their ~invariant pixel numbers)
			# this is done to potentially correct for any erroneous identification of "good" pixels within chip gaps
			# figure out chip gap pixels based on CCDSUM binning header keyword
			hdusci = pyfits.open(dispscience)
			ccdsum = (hdusci[0].header)['CCDSUM']
			if ccdsum == '2 4': # 2x4 binning
				c1min,c1max = params.chipGapPix24[0]
				c2min,c2max = params.chipGapPix24[1]
			elif ccdsum == '4 4': # 4x4 binning
				c1min,c1max = params.chipGapPix44[0]
				c2min,c2max = params.chipGapPix44[1]
			hdusci.close()
			datbwm = hdubwm[0].data.copy()
			datbwm = datbwm[0][0] # spectrum bad wavelength mask is in "band 0" (rest is sky, noise, etc.)
			datbwm[c1min:c1max] = 1
			datbwm[c2min:c2max] = 1
			hdubwm[0].data = datbwm
			hdrbwm = hdubwm[0].header
			hdrbwm.update('GR-ANGLE',float(angle))
			hdubwm.flush()
			hdubwm.close()
			pipeHistory.updatePipeKeys(inputname=namebwmsci,imagetype='bwm_sci',procChar='o')
			hdudsp = pyfits.open(dispscience,mode='update') # adds BPM keyword with BWM filename for odcombine
			hdr = hdudsp[1].header
			hdr.update('BPM',namebwmsci)
			hdr0 = hdudsp[0].header
			hdr0.update('BPM',namebwmsci)
			hdudsp.flush()
			hdudsp.close()
			dicts.bwmsciences[angle] = namebwmsci
		if dispstandard != '': # create BWM for standard star spectrum if standard star spectrum exists
			suffixlist = dispstandard.split('.')
			suffix = suffixlist[1][5:]
			namebwmstd = 'std'+str(angle)+'bwm'+suffix+'.fits'
			imutil.imexpr.unlearn()
			imutil.imexpr(expr='a==0',output=namebwmstd,a=dispstandard+'[1]')
			# this adds the 'RUIMGTYP' and 'RUPIPE' keywords to the bad wavelength mask
			hdubwm = pyfits.open(namebwmstd,mode='update')
			# this manually sets all pixels to 1 within the chip gap regions (based on their ~invariant pixel numbers)
			# this is done to potentially correct for any erroneous identification of "good" pixels within chip gaps
			# figure out chip gap pixels based on CCDSUM binning header keyword
			hdustd = pyfits.open(dispstandard)
			ccdsum = (hdustd[0].header)['CCDSUM']
			if ccdsum == '2 4': # 2x4 binning
				c1min,c1max = params.chipGapPix24[0]
				c2min,c2max = params.chipGapPix24[1]
			elif ccdsum == '4 4': # 4x4 binning
				c1min,c1max = params.chipGapPix44[0]
				c2min,c2max = params.chipGapPix44[1]
			hdustd.close()
			datbwm = hdubwm[0].data.copy()
			datbwm = datbwm[0][0] # spectrum bad wavelength mask is in "band 0" (rest is sky, noise, etc.)
			datbwm[c1min:c1max] = 1
			datbwm[c2min:c2max] = 1
			hdubwm[0].data = datbwm
			hdrbwm = hdubwm[0].header
			hdrbwm.update('GR-ANGLE',float(angle))
			hdubwm.flush()
			hdubwm.close()
			pipeHistory.updatePipeKeys(inputname=namebwmstd,imagetype='bwm_std',procChar='o')
			hdudsp = pyfits.open(dispstandard,mode='update') # adds BPM keyword with BWM filename for odcombine
			hdr = hdudsp[1].header
			hdr.update('BPM',namebwmstd)
			hdr0 = hdudsp[0].header
			hdr0.update('BPM',namebwmstd)
			hdudsp.flush()
			hdudsp.close()
			dicts.bwmstandards[angle] = namebwmstd


# This function flux-calibrates the science spectra using the standard star spectra with standard, sensfunc, and calibrate.
# This function also flux-calibrates the standard star spectra since they'll be used later for telluric corrections.
def fluxcal(dspsci=dicts.dispsciences,dspstd=dicts.dispstandards,indiv=False):
	# This prints an error and quits the function if any of the inputted dictionaries don't exist
	try:
		len(dspsci),len(dspstd)
	except:
		print 'One of the inputted dictionaries does not exist: cannot proceed.'
		return
	# This runs standard, sensfunc for each standard star spectrum; and then calibrate on the associated science image
	angles = dspstd.keys() # use angles from dspstd dictionary since we need only do flux cal if there are standards
	for angle in angles: # won't run if no standards for if not, then there are no angles
		# this stores the filename for a science image, and its header
		dspscience = dspsci[angle]
		suffixlist = dspscience.split('.')
		suffix = suffixlist[1][5:]
		hdulist = pyfits.open(dspscience,mode='update')
		hdr0_old = hdulist[0].header.copy()
		hdr1_old = hdulist[1].header.copy()
		# figure out chip gap pixels based on CCDSUM binning header keyword
		ccdsum = (hdulist[0].header)['CCDSUM']
		if ccdsum == '2 4': # 2x4 binning
			c1min,c1max = params.chipGapPix24[0]
			c2min,c2max = params.chipGapPix24[1]
		elif ccdsum == '4 4': # 4x4 binning
			c1min,c1max = params.chipGapPix44[0]
			c2min,c2max = params.chipGapPix44[1]
		hdulist.close()
		# this stores the filename of the corresponding standard star image and the class of standard star it is
		dspstandard = dspstd[angle]
		hdulist = pyfits.open(dspstandard,mode='update')
		hdr = hdulist[1].header # identify and reidentify messes around with header 0, so added keys to header 1
		starclass = hdr['OBJECT']
		airmass = hdr['AIRMASS']
		exptime = hdr['EXPTIME']
		hdulist.close()
		# this sets the output names for the flux-calibrated science image, the std file, and the sens file
		namesci = 'sci'+str(angle)+'flx'+suffix+'.fits'
		auxnamesci = 'sci'+str(angle)+'flx'+suffix+'AUX.fits'
		# This checks if a sensfile exists from a previous std star reduction for this angle, if yes, just run calibrate
		if dicts.sensfiles.get(angle,'') != '':
			print "sensfunc from previous std star reduction found for angle "+angle
			stdfile = ''
			sensfile = dicts.sensfiles[angle]
			pyrafCalls.run_calibrate(scienceimage=dspscience,fluximage=auxnamesci,sensfilename=sensfile,customRun=indiv)
			# this copies the old images to preserve file structure, and transfers the extracted data+header
			pyrafCalls.run_imcopy(inputname=dspscience,outputname=namesci)
			hdulist = pyfits.open(auxnamesci)
			datnew = hdulist[0].data.copy()
			hdr = hdulist[0].header
			hdr.update('EXTNAME','SCI')
			hdr.update('EXTVER',1)
			hdrnew = hdr.copy()
			hdulist.close()
			hdunew = pyfits.open(namesci,mode='update')
			hdunew[1].data = datnew
			hdunew[1].header = hdrnew
			hdunew.flush() # this should cause pyfits to generate a message that it's fixing some stuff
			hdr = hdunew[0].header
			hdr['RUSTD'] = stdfile # for dictionaries()
			hdr['RUSENS'] = sensfile # for dictionaries()
			hdunew.flush()
			hdunew.close()
			pipeHistory.updatePipeKeys(inputname=namesci,imagetype='science',procChar='g')
			# this adds the flux-calibrated science image filename to dictionary fluxsciences indexed by gr-angle
			dicts.fluxsciences[angle] = namesci
			# this deletes aux files
			pyrafCalls.run_imdel(auxnamesci)
			# updates user on successful flux calibration
			print 'output: '+namesci
		# no existing sensfunc: this informs the user that standard will be run automatically
		else:
			print "Running standard automatically to define bandpasses for: "+dspstandard
			print "The pipeline will automatically exclude wavelengths of chip gaps and major telluric absorption features: "
			print "Telluric bands to be excluded: 6860A to 6890A, 7170A to 7350A, and 7600A to 7630A"
			stdfile = 'std'+str(angle)+'flx'+suffix
			sensfile = 'sens'+str(angle)+'flx'+suffix
			# this runs onedspec.standard to produce a unique std file
			# get star_name to feed into run_standard; star_name is the same as the 'OBJECT' 0-header keyword (based on aschere's SALT modifications)
			hdu = pyfits.open(dspstandard)
			hdr = hdu[0].header
			objName = (hdr['OBJECT']).lower()
			hdu.close()
			# check if objName matches one of the SALT standard stars
			# This creates a list of the SALT standard stars minus the '.dat' extension
			# /usr/local/astro64/iraf/extern/pysalt/data/standards/spectroscopic/
			cwd = os.getcwd()
			os.chdir(params.standardsPath)
			possiblestandards = glob('m*.dat')
			os.chdir(cwd)
			for n,s in enumerate(possiblestandards): # remove '.dat.' extension
				snew = s.split('.')
				possiblestandards[n] = snew[0].lower()
			# this calls run_standard only if 'm'+objName is in the caldir, else skips this flux calibration
			if 'm'+objName in possiblestandards:
				pyrafCalls.run_standard(standardimage=dspstandard,outputname=stdfile,exp=exptime,air=airmass,starName='m'+objName,customRun=indiv)
			else:
				print "There is no bandpass data file in the calibration directory for: "+objName
				print "Skipping flux calibration for: "+dspscience
				print "Skipping flux calibration for: "+dspstandard
				continue
			# std file was created non-interactively with bandpasses in chip gaps and telluric absorption lines
			# go through it and delete those lines whose first column number (wavelength) falls in a chip gap or telluric absorption line
			hdu = pyfits.open(dspstandard)
			dat = hdu[1].data
			dat = dat[0][0]
			hdr = hdu[1].header
			crval = hdr['CRVAL1']
			crpix = hdr['CRPIX1']
			cd1 = hdr['CD1_1']
			pix = np.arange(len(dat))
			wave = np.linspace(crval,crval+cd1*(len(dat)-crpix),len(dat))
			wavec1 = wave[np.logical_and(pix>=c1min,pix<=c1max)] # wavelengths of first chip gap (+/- epsilon)
			wavec2 = wave[np.logical_and(pix>=c2min,pix<=c2max)] # wavelengths of second chip gap
			wavec1min = np.min(wavec1)
			wavec1max = np.max(wavec1)
			wavec2min = np.min(wavec2)
			wavec2max = np.max(wavec2)
			hdu.close()
			# now go through the std file and remove rows whose first column # (bandpass central wavelength) falls 
			# between wavec1min and wavec1max, or between wavec2min and wavec2max
			# make sure you read in dtype same as wavec1 and wavec2 else comparisons might not be meaningful 
			f = open(stdfile,mode="r") # open stdfile for reading only
			lines = f.readlines()
			f.close()
			hdrline = lines[0]
			datalines = lines[1:] # first line is header, skip it
			f = open(stdfile,mode="w") # re-open stdfile for writing (deletes all previous content)
			f.write(hdrline)
			for l in datalines:
				w = (l.split())[0] # get the wavelength (first column)
				w = float(w) # typecast w from string to float
				# only write line if w not in chip gap or telluric
				if (w>=wavec1min and w<=wavec1max) or (w>=wavec2min and w<=wavec2max) or (w>=7600.0 and w<=7630.0) or (w>=6860.0 and w<=6890.0) or (w>=7170.0 and w<=7350.0):
					continue # don't write line because this is a bad bandpass
				else:
					f.write(l)
			f.close()
			# this adds the std file to the stdfiles dictionary
			dicts.stdfiles[angle] = stdfile
			# this runs onedspec.sensfunc to produce a unique sens file
			pyrafCalls.run_sensfunc(stddata=stdfile,sensname=sensfile,customRun=indiv)
			# this adds GR-ANGLE to sensfile header for dictionaries()
			hdusens = pyfits.open(sensfile+'.fits',mode='update') 
			hdrsens = hdusens[0].header
			hdrbwm.update('GR-ANGLE',float(angle))
			hdrsens['RUPIPE'] = 'g'
			hdrsens['RUIMGTYPE'] = 'sensfunc'
			hdusens.flush()
			hdusens.close()
			# this adds the sens file to the sensfiles dictionary
			dicts.sensfiles[angle] = sensfile+'.fits.' # sensfunc takes as input just the root name for the output, not the extension (.fits)
			# this runs onedspec.calibrate to calibrate science image using the just-produced sens file (on dspscience[SCI])
			dspscience1 = dspscience+'[SCI]'
			pyrafCalls.run_calibrate(scienceimage=dspscience,fluximage=auxnamesci,sensfilename=sensfile,customRun=indiv)
			# this copies the old images to preserve file structure, and transfers the extracted data+header
			pyrafCalls.run_imcopy(inputname=dspscience,outputname=namesci)
			hdulist = pyfits.open(auxnamesci)
			datnew = hdulist[0].data.copy()
			hdr = hdulist[0].header
			hdr.update('EXTNAME','SCI')
			hdr.update('EXTVER',1)
			hdrnew = hdr.copy()
			hdulist.close()
			hdunew = pyfits.open(namesci,mode='update')
			hdunew[1].data = datnew
			hdunew[1].header = hdrnew
			hdunew.flush() # this should cause pyfits to generate a message that it's fixing some stuff
			hdr = hdunew[0].header
			hdr['RUSTD'] = stdfile # for dictionaries()
			hdr['RUSENS'] = sensfile # for dictionaries()
			hdunew.flush()
			hdunew.close()
			pipeHistory.updatePipeKeys(inputname=namesci,imagetype='science',procChar='g')
			# this adds the flux-calibrated science image filename to dictionary fluxsciences indexed by gr-angle
			dicts.fluxsciences[angle] = namesci
			# this deletes aux files
			pyrafCalls.run_imdel(auxnamesci)
			# updates user on successful flux calibration
			print 'output: '+namesci
		######## use the exact same sensfunc file to flux-calibrate the dispersion-corrected STANDARD STAR
		# set the names for the output files
		suffixliststd = dspstandard.split('.')
		suffixstd = suffixliststd[1][5:]
		nameflxstd = 'std'+str(angle)+'flx'+suffixstd+'.fits'
		auxnameflxstd = 'std'+str(angle)+'flx'+suffixstd+'AUX.fits'
		pyrafCalls.run_calibrate(scienceimage=dspstandard,fluximage=auxnameflxstd,sensfilename=sensfile,customRun=indiv)
		# this copies the old images to preserve file structure, and transfers the extracted data+header
		pyrafCalls.run_imcopy(inputname=dspstandard,outputname=nameflxstd)
		hdulist = pyfits.open(auxnameflxstd)
		datnew = hdulist[0].data.copy()
		hdr = hdulist[0].header
		hdr.update('EXTNAME','SCI')
		hdr.update('EXTVER',1)
		hdrnew = hdr.copy()
		hdulist.close()
		hdunew = pyfits.open(nameflxstd,mode='update')
		hdunew[1].data = datnew
		hdunew[1].header = hdrnew
		hdunew.flush() # this should cause pyfits to generate a message that it's fixing some stuff
		hdr = hdunew[0].header
		hdr['RUSTD'] = stdfile # for dictionaries()
		hdr['RUSENS'] = sensfile # for dictionaries()
		hdunew.flush()
		hdunew.close()
		pipeHistory.updatePipeKeys(inputname=nameflxstd,imagetype='standard',procChar='g')
		# this adds the flux-calibrated science image filename to dictionary fluxsciences indexed by gr-angle
		dicts.fluxstandards[angle] = nameflxstd
		# this deletes aux files
		pyrafCalls.run_imdel(auxnameflxstd)
		# updates user on successful flux calibration
		print 'output: '+nameflxstd
	print 'The flux-calibrated science spectra filenames are:'
	print dicts.fluxsciences
	print 'The flux-calibrated standard star spectra filenames are:'
	print dicts.fluxstandards
	


# This function adds highly deviant pixels in the science and standard star spectra
# to their respective bad wavelength masks.
# Inspired by and adapted from astropy.stats.sigma_clip(...).
def sigmaClipSpectra(dispspectra=dicts.dispsciences,fluxspectra=dicts.fluxsciences,standardspectra=dicts.fluxstandards,masksciences=dicts.bwmsciences,maskstandards=dicts.bwmstandards,indiv=False):
	try:
		len(masksciences)
	except:
		print "One of the inputted dictionaries does not exist: cannot proceed."
		return
	angles = dispspectra.keys() # clean dispersion-corrected science spectra
	answer = -1
	for angle in angles:
		dspscience = dispspectra[angle]
		bwmscience = masksciences[angle]
		hdu = pyfits.open(dspscience)
		dat = hdu[1].data.copy()
		dat = dat[0][0]
		hdr = hdu[0].header
		faintKey = hdr.get('RUFAINT','1') # 3,4 => faint extraction, 1,2 => bright extraction
		hdu.close()
		if faintKey == '3' or faintKey == '4':
			print 'WARNING: This was a FAINT reduction, possibly a host galaxy spectrum.'
			print 'Sigma-clipping the spectrum might remove strong, narrow emission lines.'
			print 'This algorithm should be okay for actual supernovae.'
			while True:
				answer = raw_input("Please enter 0 to skip sigma-clipping, or 1 to continue: ")
				if answer == '0' or answer == '1':
					break
				else:
					print "Invalid input. You must enter either 0 or 1."
		if answer == '0':
			break
		datM = np.ma.array(dat,copy=True)
		for i in range(2):
			dev = datM - np.median(datM)
			datM.mask |= dev**2 > np.var(datM)*4**2 # 4-sigma clipping with 2 iterations (tested)
		badPix = np.where(datM.mask == True)[0]
		hdubwm = pyfits.open(bwmscience,mode='update')
		datbwm = hdubwm[0].data.copy()
		for p in badPix:
			datbwm[p] = 1
		hdubwm[0].data = datbwm
		hdubwm.flush()
		hdubwm.close()
		print 'deviant pixels identified: '+dspscience # updates user on successful cleaning
	angles = fluxspectra.keys() # clean flux-calibrated science spectra
	for angle in angles:
		flxscience = fluxspectra[angle]
		if answer == '0':
			break
		elif answer == -1:
			hdu = pyfits.open(flxscience)
			hdr = hdu[0].header
			faintKey = hdr.get('RUFAINT','1') # 3,4 => faint extraction, 1,2 => bright extraction
			hdu.close()
			if faintKey == '3' or faintKey == '4':
				print 'WARNING: This was a FAINT reduction, possibly a host galaxy spectrum.'
				print 'Sigma-clipping the spectrum might remove strong, narrow emission lines.'
				print 'This algorithm should be okay for actual supernovae.'
				while True:
					answer = raw_input("Please enter 0 to skip sigma-clipping, or 1 to continue: ")
					if answer == '0' or answer == '1':
						break
					else:
						print "Invalid input. You must enter either 0 or 1."
		if answer == '0':
			break
		bwmscience = masksciences[angle]
		hdu = pyfits.open(flxscience)
		dat = hdu[1].data.copy()
		dat = dat[0][0]
		hdu.close()
		datM = np.ma.array(dat,copy=True)
		for i in range(2):
			dev = datM - np.median(datM)
			datM.mask |= dev**2 > np.var(datM)*4**2 # 4-sigma clipping with 2 iterations (tested)
		badPix = np.where(datM.mask == True)[0]
		hdubwm = pyfits.open(bwmscience,mode='update')
		datbwm = hdubwm[0].data.copy()
		for p in badPix:
			datbwm[p] = 1
		hdubwm[0].data = datbwm
		hdubwm.flush()
		hdubwm.close()	
		print 'deviant pixels identified: '+flxscience # updates user on successful cleaning
	angles = standardspectra.keys() # clean flux-calibrated standard star spectra
	for angle in angles:
		flxstandard = standardspectra[angle]
		bwmstandard = maskstandards[angle]
		hdu = pyfits.open(flxstandard)
		dat = hdu[1].data.copy()
		dat = dat[0][0]
		hdu.close()
		datM = np.ma.array(dat,copy=True)
		for i in range(2):
			dev = datM - np.median(datM)
			datM.mask |= dev**2 > np.var(datM)*4**2 # 4-sigma clipping with 2 iterations (tested)
		badPix = np.where(datM.mask == True)[0]
		hdubwm = pyfits.open(bwmstandard,mode='update')
		datbwm = hdubwm[0].data.copy()
		for p in badPix:
			datbwm[p] = 1
		hdubwm[0].data = datbwm
		hdubwm.flush()
		hdubwm.close()	
		print 'deviant pixels identified: '+flxstandard # updates user on successful cleaning
	

# This function scales flux-calibrated science spectra so that they have similar flux values at overlapping wavelengths.
def scaleFluxScienceSpectra(sciencespectra=dicts.fluxsciences,indiv=False):
	# This prints an error and quits the function if any of the inputted dictionaries don't exist
	try:
		len(sciencespectra)
	except:
		print 'One of the inputted dictionaries does not exist: cannot proceed.'
		return
	# This begins the flux-scaling process for the flux-calibrated science spectra.
	angles = sciencespectra.keys()
	if len(angles) == 1: # if only one flux-calibrated spectrum, no other spectra to compare to
		print 'There is only one flux-calibrated spectrum -- the one with angle: '+angles[0]
		print 'Cannot scale a single spectrum.'
		print 'Skipping flux scaling.'
		return
	elif len(angles) == 2: # only two flux-calibrated spectra => compare these two only
		print 'There are only two flux-calibrated spectra to scale.'
		print 'Running this task non-interactively.'
		angle0 = angles[0] # default assumption: 1st element of angles is the lower angle number
		angle1 = angles[1] 
		angOne = float(angles[0]) # angles are strings currently
		angTwo = float(angles[1]) # converting to floats to check the default assumption
		if angTwo < angOne: # need to make this comparison to create appropriate common wavelength array
			angle0 = angles[1]
			angle1 = angles[0]
		sci0 = sciencespectra[angle0] # science spectrum 0 (with lower gr-angle) -- this will be scaled
		sci1 = sciencespectra[angle1] # science spectrum 1 (with greater/equal gr-angle) -- this is the reference
		# retrieving necessary information for science spectrum 0
		hdu0 = pyfits.open(sci0)
		hdr0 = hdu0[1].header
		dat0 = hdu0[1].data.copy()
		dat0 = dat0[0][0] # data for science spectrum itself is the first "band" (rest are noise, sky, etc.)
		crpix0 = hdr0['CRPIX1'] # dispersion (WCS) information for science spectrum 0
		cd0 = hdr0['CD1_1']
		crval0 = hdr0['CRVAL1']
		crpix0 = hdr0['CRPIX1']
		hdu0.close()
		wave0 = np.linspace(crval0,crval0+cd0*(len(dat0)-crpix0),len(dat0)) # wavelength array for science spectrum 0
		# retrieving necessary information for science spectrum 1
		hdu1 = pyfits.open(sci1)
		hdr1 = hdu1[1].header
		dat1 = hdu1[1].data.copy()
		dat1 = dat1[0][0] # data for science spectrum itself is the first "band" (rest are noise, sky, etc.)
		crpix1 = hdr1['CRPIX1'] # dispersion (WCS) information for science spectrum 1
		cd1 = hdr1['CD1_1']
		crval1 = hdr1['CRVAL1']
		# figure out chip gap pixels based on CCDSUM binning header keyword (same for dat0 and dat1)
		ccdsum = (hdu1[0].header)['CCDSUM']
		if ccdsum == '2 4': # 2x4 binning
			c1min,c1max = params.chipGapPix24[0]
			c2min,c2max = params.chipGapPix24[1]
		elif ccdsum == '4 4': # 4x4 binning
			c1min,c1max = params.chipGapPix44[0]
			c2min,c2max = params.chipGapPix44[1]
		hdu1.close()
		wave1 = np.linspace(crval1,crval1+cd1*(len(dat1)-crpix1),len(dat1)) # wavelength array for science spectrum 1
		# get rid of the chip gaps in spectrum 0
		pix0 = np.arange(len(dat0))+1 # 1-based pixel array for dat0
		dat0good = dat0[np.logical_or(pix0<=c1min,pix0>=c1max)] # remove flux values of first chip gap for spectrum 0
		wave0good = wave0[np.logical_or(pix0<=c1min,pix0>=c1max)] # remove wavelengths of first chip gap for spectrum 0
		pix0good = pix0[np.logical_or(pix0<=c1min,pix0>=c1max)] # remove pixels of first chip gap for spectrum 0
		dat0good = dat0good[np.logical_or(pix0good<=c2min,pix0good>=c2max)] # remove flux values of second chip gap for spectrum 0
		wave0good = wave0good[np.logical_or(pix0good<=c2min,pix0good>=c2max)] # remove wavelengths of second chip gap for spectrum 0
		pix0good = pix0good[np.logical_or(pix0good<=c2min,pix0good>=c2max)] # remove pixels of second chip gap for spectrum 0
		# get rid of the chip gaps in spectrum 1
		pix1 = np.arange(len(dat1))+1 # 1-based pixel array for dat1
		dat1good = dat1[np.logical_or(pix1<=c1min,pix1>=c1max)] # remove flux values of first chip gap for spectrum 1
		wave1good = wave1[np.logical_or(pix1<=c1min,pix1>=c1max)] # remove wavelengths of first chip gap for spectrum 1
		pix1good = pix1[np.logical_or(pix1<=c1min,pix1>=c1max)] # remove pixels of first chip gap for spectrum 1
		dat1good = dat1good[np.logical_or(pix1good<=c2min,pix1good>=c2max)] # remove flux values of second chip gap for spectrum 1
		wave1good = wave1good[np.logical_or(pix1good<=c2min,pix1good>=c2max)] # remove wavelengths of second chip gap for spectrum 1
		pix1good = pix1good[np.logical_or(pix1good<=c2min,pix1good>=c2max)] # remove pixels of second chip gap for spectrum 1
		# construct common wavelength array
		minwave0g = min(wave0good)
		maxwave0g = max(wave0good)
		minwave1g = min(wave1good)
		maxwave1g = max(wave1good)
		wavecom = wave1good[np.logical_and(wave1good>=minwave1g,wave1good<=maxwave0g)] # common wavelength array constructed from wave1
		# find flux values of spectrum 1 at common wavelengths
		dat1com = dat1good[np.logical_and(wave1good>=minwave1g,wave1good<=maxwave0g)] # common-wavelength flux values of spectrum 1
		# find interpolant of entire chip-gap-corrected wavelength values of spectrum 0
		f0 = interp1d(wave0good,dat0good)
		# use interpolant to get flux values of spectrum 0 at common wavelengths only
		dat0com = f0(wavecom)
		# use least-squares method to minimize distance between spectrum 0 and spectrum 1 flux values at common wavelengths
		print 'A second order scaling function will now be fit to reduce the error between the two spectra at common wavelengths: '
		print 'p[0] + p[1]*x + p[2]*x*x where x is the common wavelength NumPy array.'
		fitfunc = lambda p,x: (p[0] + p[1]*x + p[2]*x*x)*dat0com # this order-2 function will scale spectrum 0 to spectrum 1 flux values
		errfunc = lambda p,x,y: fitfunc(p,x) - y
		p0 = [1.0, 0.0, 0.0] # initial guess for parameters of second-order fitfunc is a constant function
		p1,success = leastsq(errfunc,p0,args=(wavecom,dat1com))
		print 'The parameters of the second-order scaling function are ([p0,p1,p2]): '
		print p1
		# apply the scaling function to entire original spectrum 0 data 
		dat0scaled = (p1[0]+p1[1]*wave0+p1[2]*wave0*wave0)*dat0
		# copy old file over and replace old data with new scaled data
		suffixlist = sci0.split('.')
		suffix = suffixlist[1][5:]
		namesci = 'sci'+angle0+'scl'+suffix+'.fits'
		pyrafCalls.run_imcopy(inputname=sci0,outputname=namesci)
		# we will replace the 1st spectral band (optimally-extracted spectrum) with the scaled one
		# band 2 will remain as the non-optimal, uncleaned, unweighted, non-scaled spectrum
		# band 3 will remain the non-scaled sky spectrum
		# band 4 will remain the non-scaled sigma spectrum
		hdu0 = pyfits.open(namesci,mode='update')
		datsci = hdu0[1].data.copy()
		datsci[0] = dat0scaled
		hdu0[1].data = datsci
		hdu0.flush()
		hdu0.close()
		pipeHistory.updatePipeKeys(inputname=namesci,imagetype='science',procChar='n')
		# add scaled science spectrum filename to dictionary scaledsciences
		dicts.scaledfluxsciences[angle0] = namesci
	elif len(angles) > 2: # for multiple flux-calibrated spectra, ask the user for a gr-angle and/or do a pairwise loop
		print "There are more than two flux-calibrated spectra."
		print "You have the option of picking one spectrum's GR-ANGLE."
		print "The flux values at overlapping wavelengths for all other spectra will be scaled to this spectrum."
		print "Here is the global dictionary of flux-calibrated science images:"
		print sciencespectra
		# for a non-interactive version, the median angle could just be chosen automatically
		# for an even number of spectra, the higher of the two median angles could be chosen (more flux than bluer angles)
		while True:
			refAngle = raw_input("Please enter a (preferably MEDIAN) GR-ANGLE: ")
			if refAngle in angles: 
				print "You chose to scale all other spectra with reference to: "+sciencespectra[refAngle]
				break
			else:
				print "Invalid input: you must enter an angle in the above dictionary to proceed."
		# all other spectra must now be scaled based on sciencespectra[refAngle]
		# all spectra should have overlapping wavelengths with each other
		# it is this common wavelength range, different for each pair, that will be used to find the optimal scaling
		# a second-order polynomial/function will be fit to each of the spectra that need to be scaled
		# first, gather and construct the necessary information about the reference spectrum
		refsci = sciencespectra[refAngle]
		hduref = pyfits.open(refsci)
		datr = hduref[1].data.copy() # reference spectrum data
		datr = datr[0][0] # science data is in "band 1" (rest is noise, sky, etc.)
		hdrr = hduref[1].header # reference spectrum header
		crpixr= hdrr['CRPIX1'] # reference spectrum WCS information
		cdr = hdrr['CD1_1']
		crvalr = hdrr['CRVAL1']
		hduref.close()
		# next, use a for loop to create a scaled spectrum for each non-refsci spectrum
		for angle in angles:
			if angle == refAngle: 
				continue # we don't want to scale the reference spectrum to itself, so skip it
			waver = np.linspace(crvalr,crvalr+cdr*(len(datr)-crpixr),len(datr)) # reference spectrum wavelength array
			# retrieve and construct the necessary information about the current spectrum to be scaled
			cursci = sciencespectra[angle] # current science spectrum
			print "The current spectrum being scaled is: "+cursci
			hducur = pyfits.open(cursci)
			datc = hducur[1].data.copy() # current science spectrum data
			datc = datc[0][0] # science data is "band 1" (rest is noise, sky, etc.)
			hdrc = hducur[1].header # current science spectrum header
			crpixc = hdrc['CRPIX1'] # current science spectrum WCS information
			cdc = hdrc['CD1_1']
			crvalc = hdrc['CRVAL1']
			# figure out chip gap pixels based on CCDSUM binning header keyword
			ccdsum = (hducur[0].header)['CCDSUM']
			if ccdsum == '2 4': # 2x4 binning
				c1min,c1max = params.chipGapPix24[0]
				c2min,c2max = params.chipGapPix24[1]
			elif ccdsum == '4 4': # 4x4 binning
				c1min,c1max = params.chipGapPix44[0]
				c2min,c2max = params.chipGapPix44[1]
			hducur.close()
			wavec = np.linspace(crvalc,crvalc+cdc*(len(datc)-crpixc),len(datc)) # current science spectrum wavelength array
			# get rid of the chip gaps (based on their ~invariant pixel numbers) for reference spectrum
			pixr = np.arange(len(datr))+1 # 1-based pixel array for reference spectrum
			datrgood = datr[np.logical_or(pixr<=c1min,pixr>=c1max)] # remove first chip gap's flux values
			wavergood = waver[np.logical_or(pixr<=c1min,pixr>=c1max)] # remove first chip gap's wavelengths
			pixrgood = pixr[np.logical_or(pixr<=c1min,pixr>=c1max)] # remove first chip gap's pixels
			datrgood = datrgood[np.logical_or(pixrgood<=c2min,pixrgood>=c2max)] # remove second chip gap's flux values
			wavergood = wavergood[np.logical_or(pixrgood<=c2min,pixrgood>=c2max)] # remove second chip gap's wavelengths
			pixrgood = pixrgood[np.logical_or(pixrgood<=c2min,pixrgood>=c2max)] # remove second chip gap's pixels
			# get rid of the chip gaps (based on their ~invariant pixel numbers) for current spectrum
			pixc = np.arange(len(datc))+1 # 1-based pixel array for current spectrum
			datcgood = datc[np.logical_or(pixc<=c1min,pixc>=c1max)] # remove first chip gap's flux values
			wavecgood = wavec[np.logical_or(pixc<=c1min,pixc>=c1max)] # remove first chip gap's wavelengths
			pixcgood = pixc[np.logical_or(pixc<=c1min,pixc>=c1max)] # remove first chip gap's pixels
			datcgood = datcgood[np.logical_or(pixcgood<=c2min,pixcgood>=c2max)] # remove second chip gap's flux values
			wavecgood = wavecgood[np.logical_or(pixcgood<=c2min,pixcgood>=c2max)] # remove second chip gap's wavelengths
			pixcgood = pixcgood[np.logical_or(pixcgood<=c2min,pixcgood>=c2max)] # remove second chip gap's pixels
			# construct a common wavelength array
			minwavergood = min(wavergood) # minimum wavelength of chip-corrected reference spectrum
			maxwavergood = max(wavergood) # maximum wavelength of chip-corrected ref spectrum
			minwavecgood = min(wavecgood) # minimum wavelength of chip-corrected current spectrum
			maxwavecgood = max(wavecgood) # maximum wavelength of chip-corrected cur spectrum
			# construct common wavelength array based on whether angle or refAngle is lower
			# and then find flux values of reference spectrum at these common wavelengths
			if angle <= refAngle: # the two angles should not really ever be equal (for this step)
				wavecom = wavergood[np.logical_and(wavergood>=minwavergood,wavergood<=maxwavecgood)]
				datrcom = datrgood[np.logical_and(wavergood>=minwavergood,wavergood<=maxwavecgood)]
			elif angle > refAngle:
				wavecom = wavergood[np.logical_and(wavergood>=minwavecgood,wavergood<=maxwavergood)]
				datrcom = datrgood[np.logical_and(wavergood>=minwavecgood,wavergood<=maxwavergood)]
			# find interpolant of entire chip-gap-corrected wavelength values of current spectrum
			fc = interp1d(wavecgood,datcgood)
			# use interpolant to get flux values of current spectrum at common wavelengths only
			datccom = fc(wavecom)
			# use least-squares method to minimize distance between ref and cur spectra at common wavelengths
			print 'A first order scaling function will now be fit to reduce the error between the two spectra at common wavelengths: '
			print 'p[0] + p[1]*x where x is the common wavelength NumPy array.'
			fitfunc = lambda p,x: (p[0] + p[1]*x)*datccom # this order-1 function will scale cur spectrum
			errfunc = lambda p,x,y: fitfunc(p,x) - y
			p0 = [1.0, 0.0] # initial guess for parameters of first-order fitfunc is that it is a constant function
			p1,success = leastsq(errfunc,p0,args=(wavecom,datrcom))
			print 'The parameters of the second-order scaling function are ([p0,p1,p2]): '
			print p1
			# apply the scaling function to entire original current spectrum data 
			datcscaled = (p1[0]+p1[1]*wavec)*datc
			# copy old file over and replace old data with new scaled data
			suffixlist = cursci.split('.')
			suffix = suffixlist[1][5:]
			namesci = 'sci'+angle+'scl'+suffix+'.fits'
			pyrafCalls.run_imcopy(inputname=cursci,outputname=namesci)
			hdunew = pyfits.open(namesci,mode='update')
			datsci = hdunew[1].data
			# we will replace the 1st spectral band (optimally-extracted spectrum) with the scaled one
			# band 2 will remain as the non-optimal, uncleaned, unweighted, non-scaled spectrum
			# band 3 will remain the non-scaled sky spectrum
			# band 4 will remain the non-scaled sigma spectrum
			datsci[0] = datcscaled
			hdunew[1].data = datsci
			hdunew.flush()
			hdunew.close()
			pipeHistory.updatePipeKeys(inputname=namesci,imagetype='science',procChar='n')
			# add scaled science spectrum filename to dictionary scaledfluxsciences
			dicts.scaledfluxsciences[angle] = namesci 
		# begin second round (r2) of scaling for the extreme gr-angles which may have had insufficient wavelength overlap with refAngle
		rescaledAngles = []
		r2angles = dicts.scaledfluxsciences.keys() # added after the fact: these angles don't contain refAngle due to round 1
		floatAngles = []
		floatAngMaps = {} # will be needed to access keys (string format) of scaledsciences
		for i,v in enumerate(r2angles): # convert angles from strings to floats for numerical comparison
			floatAngles.append(float(v))
			floatAngMaps[float(v)] = v
		refAngleFloat = float(refAngle) # convert refAngle from string to float for numerical comparison
		for r2ang1 in floatAngles:
			if r2ang1 == refAngleFloat:
				continue # this is pairwise rescaling based on the angles with refAngle excluded
			if r2ang1 in rescaledAngles:
				continue # each time an angle is rescaled, it is added to a list so that it is not rescaled again
			# check if r2ang1 spectrum is far enough from refAngle spectrum such that there was insufficient wavelength overlap in round 1
			# a difference in gr-angle of 2 degrees is arbitrary
			# I wonder if it can be made more dynamic or derived through a simple statistical value (median, mean, etc.)
			if abs(r2ang1-refAngleFloat) > 2.0:
				for r2ang2 in floatAngles: # nested for loop to get the other non-refAngle spectrum to scale (or become reference)
					if r2ang2 == refAngleFloat:
						continue # again, don't want to rescale refAngle spectrum. it is the basis for our automatic scaling procedure
					# r2ang2 can be in list rescaledAngles since we may need some of them for reference
					if abs(r2ang2 - r2ang1) < 2.0: # need r2ang1 and r2ang2 spectra to have significant wavelength overlap for meaningful rescaling
						# we want to use the spectrum whose angle is closer to refAngle (preferably median) as the new reference spectrum
						# i.e., we will scale the spectrum with the more distant angle (relative to refAngle) based on the other angle
						distAng1 = abs(r2ang1 - refAngleFloat)
						distAng2 = abs(r2ang2 - refAngleFloat)
						if distAng1 <= distAng2: # if distance between r2ang1 and refAngle is shorter, use r2ang1 as the new refAngle
							newRefAngle = r2ang1
							rescaleAngle = r2ang2
						else:
							newRefAngle = r2ang2
							rescaleAngle = r2ang1
						# begin the rescaling process in the same fashion as the scaling process from round 1
						print 'Beginning RESCALING process for spectrum with GR-ANGLE: '+str(rescaleAngle)
						refsci = dicts.scaledfluxsciences[floatAngMaps[newRefAngle]] # using scaled spectra now for the rescaling process
						cursci = dicts.scaledfluxsciences[floatAngMaps[rescaleAngle]]
						# get data and relevant header information
						hduref = pyfits.open(refsci)
						datr = hduref[1].data.copy() # reference spectrum data
						datr = datr[0][0]
						hdrr = hduref[1].header # reference spectrum header
						crpixr= hdrr['CRPIX1'] # reference spectrum WCS information
						cdr = hdrr['CD1_1']
						crvalr = hdrr['CRVAL1']
						hduref.close()
						waver = np.linspace(crvalr,crvalr+cdr*(len(datr)-crpixr),len(datr)) # reference spectrum wavelength array
						# retrieve and construct the necessary information about the current spectrum to be scaled
						print "The current spectrum being rescaled is: "+cursci
						hducur = pyfits.open(cursci)
						datc = hducur[1].data.copy() # current science spectrum data
						datc = datc[0][0]
						hdrc = hducur[1].header # current science spectrum header
						crpixc = hdrc['CRPIX1'] # current science spectrum WCS information
						cdc = hdrc['CD1_1']
						crvalc = hdrc['CRVAL1']
						# figure out chip gap pixels based on CCDSUM binning header keyword
						ccdsum = (hducur[0].header)['CCDSUM']
						if ccdsum == '2 4': # 2x4 binning
							c1min,c1max = params.chipGapPix24[0]
							c2min,c2max = params.chipGapPix24[1]
						elif ccdsum == '4 4': # 4x4 binning
							c1min,c1max = params.chipGapPix44[0]
							c2min,c2max = params.chipGapPix44[1]
						hducur.close()
						wavec = np.linspace(crvalc,crvalc+cdc*(len(datc)-crpixc),len(datc)) # current science spectrum wavelength array
						# even though these are scaled spectra, the chip gaps are not gone so follow similar procedure as in round 1
						# get rid of the chip gaps (based on their ~invariant pixel numbers) for reference spectrum
						pixr = np.arange(len(datr))+1 # 1-based pixel array for reference spectrum
						datrgood = datr[np.logical_or(pixr<=c1min,pixr>=c1max)] # remove first chip gap's flux values
						wavergood = waver[np.logical_or(pixr<=c1min,pixr>=c1max)] # remove first chip gap's wavelengths
						pixrgood = pixr[np.logical_or(pixr<=c1min,pixr>=c1max)] # remove first chip gap's pixels
						datrgood = datrgood[np.logical_or(pixrgood<=c2min,pixrgood>=c2max)] # remove second chip gap's flux values
						wavergood = wavergood[np.logical_or(pixrgood<=c2min,pixrgood>=c2max)] # remove second chip gap's wavelengths
						pixrgood = pixrgood[np.logical_or(pixrgood<=c2min,pixrgood>=c2max)] # remove second chip gap's pixels
						# get rid of the chip gaps (based on their ~invariant pixel numbers) for current spectrum
						pixc = np.arange(len(datc))+1 # 1-based pixel array for current spectrum
						datcgood = datc[np.logical_or(pixc<=c1min,pixc>=c1max)] # remove first chip gap's flux values
						wavecgood = wavec[np.logical_or(pixc<=c1min,pixc>=c1max)] # remove first chip gap's wavelengths
						pixcgood = pixc[np.logical_or(pixc<=c1min,pixc>=c1max)] # remove first chip gap's pixels
						datcgood = datcgood[np.logical_or(pixcgood<=c2min,pixcgood>=c2max)] # remove second chip gap's flux values
						wavecgood = wavecgood[np.logical_or(pixcgood<=c2min,pixcgood>=c2max)] # remove second chip gap's wavelengths
						pixcgood = pixcgood[np.logical_or(pixcgood<=c2min,pixcgood>=c2max)] # remove second chip gap's pixels
						# construct a common wavelength array
						minwavergood = min(wavergood) # minimum wavelength of chip-corrected reference spectrum
						maxwavergood = max(wavergood) # maximum wavelength of chip-corrected ref spectrum
						minwavecgood = min(wavecgood) # minimum wavelength of chip-corrected current spectrum
						maxwavecgood = max(wavecgood) # maximum wavelength of chip-corrected cur spectrum
						# construct common wavelength array based on whether rescaleAngle or newRefAngle is lower
						# and then find flux values of reference spectrum at these common wavelengths
						if rescaleAngle <= newRefAngle: # the two angles should not really ever be equal (for this step)
							wavecom = wavergood[np.logical_and(wavergood>=minwavergood,wavergood<=maxwavecgood)]
							datrcom = datrgood[np.logical_and(wavergood>=minwavergood,wavergood<=maxwavecgood)]
						elif rescaleAngle > newRefAngle:
							wavecom = wavergood[np.logical_and(wavergood>=minwavecgood,wavergood<=maxwavergood)]
							datrcom = datrgood[np.logical_and(wavergood>=minwavecgood,wavergood<=maxwavergood)]
						# find interpolant of entire chip-gap-corrected wavelength values of current spectrum
						fc = interp1d(wavecgood,datcgood)
						# use interpolant to get flux values of current spectrum at common wavelengths only
						datccom = fc(wavecom)
						# use least-squares method to minimize distance between ref and cur spectra at common wavelengths
						print 'Beginning optimization process for RESCALING.'
						print 'A second order scaling function will now be fit to reduce the error between the two spectra at common wavelengths: '
						print 'p[0] + p[1]*x + p[2]*x*x where x is the common wavelength NumPy array.'
						fitfunc = lambda p,x: (p[0] + p[1]*x + p[2]*x*x)*datccom # this order-2 function will scale cur spectrum
						errfunc = lambda p,x,y: fitfunc(p,x) - y
						p0 = [1.0, 0.0, 0.0] # initial guess for parameters of second-order fitfunc is that it is a constant function
						p1,success = leastsq(errfunc,p0,args=(wavecom,datrcom))
						print 'The parameters of the second-order scaling function are ([p0,p1,p2]): '
						print p1
						# apply the scaling function to entire original (already scaled) current spectrum data 
						datcscaled = (p1[0]+p1[1]*wavec+p1[2]*wavec*wavec)*datc
						# replace old already-scaled data with new RESCALED data
						hdunew = pyfits.open(cursci,mode='update')
						datsci = hdunew[1].data.copy()
						# we will replace the 1st spectral band (optimally-extracted spectrum) with the scaled one
						# band 2 will remain as the non-optimal, uncleaned, unweighted, non-scaled spectrum
						# band 3 will remain the non-scaled sky spectrum
						# band 4 will remain the non-scaled sigma spectrum
						datsci[0] = datcscaled
						hdunew[1].data = datsci
						hdunew.flush()
						# add rescaled spectrum's angle to list rescaledAngles
						rescaledAngles.append(rescaleAngle)
						print 'Rescaling complete for spectrum with GR-ANGLE: '+str(rescaleAngle)
		print "The angles of the rescaled spectra are: "
		print rescaledAngles
	print 'The scaled flux-calibrated science spectra filenames are: '
	print dicts.scaledfluxsciences
	
	
# This function scales flux-calibrated standard star spectra so that they have similar flux values at overlapping wavelengths.
# This is necessary to make sure telluric corrections work as intended.
def scaleStdStarSpectra(stdspectra=dicts.fluxstandards,indiv=False):
	# This prints an error and quits the function if any of the inputted dictionaries don't exist
	try:
		len(stdspectra)
	except:
		print 'One of the inputted dictionaries does not exist: cannot proceed.'
		return
	# This begins the flux-scaling process for the flux-calibrated standard star spectra.
	angles = stdspectra.keys()
	if len(angles) == 1: # if only one flux-calibrated spectrum, no other spectra to compare to
		print 'There is only one flux-calibrated spectrum -- the one with angle: '+angle
		print 'Cannot scale a single spectrum.'
		print 'Skipping flux scaling.'
		return
	elif len(angles) == 2: # only two flux-calibrated spectra => compare these two only
		print 'There are only two flux-calibrated spectra to scale.'
		print 'Running this task non-interactively.'
		angle0 = angles[0] # default assumption: 1st element of angles is the lower angle number
		angle1 = angles[1] 
		angOne = float(angles[0]) # angles are strings currently
		angTwo = float(angles[1]) # converting to floats to check the default assumption
		if angTwo < angOne: # need to make this comparison to create appropriate common wavelength array
			angle0 = angles[1]
			angle1 = angles[0]
		std0 = stdspectra[angle0] # std spectrum 0 (with lower gr-angle) -- this will be scaled
		std1 = stdspectra[angle1] # std spectrum 1 (with greater/equal gr-angle) -- this is the reference
		# retrieving necessary information for std spectrum 0
		hdu0 = pyfits.open(std0)
		hdr0 = hdu0[1].header
		dat0 = hdu0[1].data.copy()
		dat0 = dat0[0][0] # data for std spectrum itself is the first "band" (rest are noise, sky, etc.)
		crpix0 = hdr0['CRPIX1'] # dispersion (WCS) information for std spectrum 0
		cd0 = hdr0['CD1_1']
		crval0 = hdr0['CRVAL1']
		hdu0.close()
		wave0 = np.linspace(crval0,crval0+cd0*(len(dat0)-crpix0),len(dat0)) # wavelength array for std spectrum 0
		# retrieving necessary information for std spectrum 1
		hdu1 = pyfits.open(std1)
		hdr1 = hdu1[1].header
		dat1 = hdu1[1].data.copy()
		dat1 = dat1[0][0] # data for std spectrum itself is the first "band" (rest are noise, sky, etc.)
		crpix1 = hdr1['CRPIX1'] # dispersion (WCS) information for std spectrum 1
		cd1 = hdr1['CD1_1']
		crval1 = hdr1['CRVAL1']
		# figure out chip gap pixels based on CCDSUM binning header keyword
		ccdsum = (hdu1[0].header)['CCDSUM']
		if ccdsum == '2 4': # 2x4 binning
			c1min,c1max = params.chipGapPix24[0]
			c2min,c2max = params.chipGapPix24[1]
		elif ccdsum == '4 4': # 4x4 binning
			c1min,c1max = params.chipGapPix44[0]
			c2min,c2max = params.chipGapPix44[1]
		hdu1.close()
		wave1 = np.linspace(crval1,crval1+cd1*(len(dat1)-crpix1),len(dat1)) # wavelength array for std spectrum 1
		# get rid of the chip gaps in spectrum 0
		pix0 = np.arange(len(dat0))+1 # 1-based pixel array for dat0
		dat0good = dat0[np.logical_or(pix0<=c1min,pix0>=c1max)] # remove flux values of first chip gap for spectrum 0
		wave0good = wave0[np.logical_or(pix0<=c1min,pix0>=c1max)] # remove wavelengths of first chip gap for spectrum 0
		pix0good = pix0[np.logical_or(pix0<=c1min,pix0>=c1max)] # remove pixels of first chip gap for spectrum 0
		dat0good = dat0good[np.logical_or(pix0good<=c2min,pix0good>=c2max)] # remove flux values of second chip gap for spectrum 0
		wave0good = wave0good[np.logical_or(pix0good<=c2min,pix0good>=c2max)] # remove wavelengths of second chip gap for spectrum 0
		pix0good = pix0good[np.logical_or(pix0good<=c2min,pix0good>=c2max)] # remove pixels of second chip gap for spectrum 0
		# get rid of the chip gaps in spectrum 1
		pix1 = np.arange(len(dat1))+1 # 1-based pixel array for dat1
		dat1good = dat1[np.logical_or(pix1<=c1min,pix1>=c1max)] # remove flux values of first chip gap for spectrum 1
		wave1good = wave1[np.logical_or(pix1<=c1min,pix1>=c1max)] # remove wavelengths of first chip gap for spectrum 1
		pix1good = pix1[np.logical_or(pix1<=c1min,pix1>=c1max)] # remove pixels of first chip gap for spectrum 1
		dat1good = dat1good[np.logical_or(pix1good<=c2min,pix1good>=c2max)] # remove flux values of second chip gap for spectrum 1
		wave1good = wave1good[np.logical_or(pix1good<=c2min,pix1good>=c2max)] # remove wavelengths of second chip gap for spectrum 1
		pix1good = pix1good[np.logical_or(pix1good<=c2min,pix1good>=c2max)] # remove pixels of second chip gap for spectrum 1
		# construct common wavelength array
		minwave0g = min(wave0good)
		maxwave0g = max(wave0good)
		minwave1g = min(wave1good)
		maxwave1g = max(wave1good)
		wavecom = wave1good[np.logical_and(wave1good>=minwave1g,wave1good<=maxwave0g)] # common wavelength array constructed from wave1
		# find flux values of spectrum 1 at common wavelengths
		dat1com = dat1good[np.logical_and(wave1good>=minwave1g,wave1good<=maxwave0g)] # common-wavelength flux values of spectrum 1
		# find interpolant of entire chip-gap-corrected wavelength values of spectrum 0
		f0 = interp1d(wave0good,dat0good)
		# use interpolant to get flux values of spectrum 0 at common wavelengths only
		dat0com = f0(wavecom)
		# use least-squares method to minimize distance between spectrum 0 and spectrum 1 flux values at common wavelengths
		print 'A second order scaling function will now be fit to reduce the error between the two spectra at common wavelengths: '
		print 'p[0] + p[1]*x + p[2]*x*x where x is the common wavelength NumPy array.'
		fitfunc = lambda p,x: (p[0] + p[1]*x + p[2]*x*x)*dat0com # this order-2 function will scale spectrum 0 to spectrum 1 flux values
		errfunc = lambda p,x,y: fitfunc(p,x) - y
		p0 = [1.0, 0.0, 0.0] # initial guess for parameters of second-order fitfunc is a constant function
		p1,success = leastsq(errfunc,p0,args=(wavecom,dat1com))
		print 'The parameters of the second-order scaling function are ([p0,p1,p2]): '
		print p1
		# apply the scaling function to entire original spectrum 0 data 
		dat0scaled = (p1[0]+p1[1]*wave0+p1[2]*wave0*wave0)*dat0
		# copy old file over and replace old data with new scaled data
		suffixlist = std0.split('.')
		suffix = suffixlist[1][5:]
		namestd = 'std'+angle0+'scl'+suffix+'.fits'
		pyrafCalls.run_imcopy(inputname=std0,outputname=namestd)
		hdu0 = pyfits.open(namestd,mode='update')
		datstd = hdu0[1].data.copy()
		# we will replace the 1st spectral band (optimally-extracted spectrum) with the scaled one
		# band 2 will remain as the non-optimal, uncleaned, unweighted, non-scaled spectrum
		# band 3 will remain the non-scaled sky spectrum
		# band 4 will remain the non-scaled sigma spectrum
		datstd[0] = dat0scaled
		hdu0[1].data = datstd
		hdu0.flush()
		hdu0.close()
		pipeHistory.updatePipeKeys(inputname=namestd,imagetype='standard',procChar='n')
		# add scaled std spectrum filename to dictionary scaledstandards
		dicts.scaledstandards[angle0] = namestd
	elif len(angles) > 2: # for multiple flux-calibrated spectra, ask the user for a gr-angle and/or do a pairwise loop
		print "There are more than two flux-calibrated spectra."
		print "You have the option of picking one spectrum's GR-ANGLE."
		print "The flux values at overlapping wavelengths for all other spectra will be scaled to this spectrum."
		print "Here is the global dictionary of flux-calibrated standard star spectra:"
		print stdspectra
		# for a non-interactive version, the median angle could just be chosen automatically
		# for an even number of spectra, the higher of the two median angles could be chosen (more flux than bluer angles)
		while True:
			refAngle = raw_input("Please enter a (preferably MEDIAN) GR-ANGLE: ")
			if refAngle in angles: 
				print "You chose to scale all other spectra with reference to: "+stdspectra[refAngle]
				break
			else:
				print "Invalid input: you must enter an angle in the above dictionary to proceed."
		# all other spectra must now be scaled based on stdspectra[refAngle]
		# all spectra should have overlapping wavelengths with each other
		# it is this common wavelength range, different for each pair, that will be used to find the optimal scaling
		# a second-order polynomial/function will be fit to each of the spectra that need to be scaled
		# first, gather and construct the necessary information about the reference spectrum
		refstd = stdspectra[refAngle]
		hduref = pyfits.open(refstd)
		datr = hduref[1].data.copy() # reference spectrum data
		datr = datr[0][0] # std data is in "band 1" (rest is noise, sky, etc.)
		hdrr = hduref[1].header # reference spectrum header
		crpixr= hdrr['CRPIX1'] # reference spectrum WCS information
		cdr = hdrr['CD1_1']
		crvalr = hdrr['CRVAL1']
		hduref.close()
		# next, use a for loop to create a scaled spectrum for each non-refstd spectrum
		for angle in angles:
			if angle == refAngle: 
				continue # we don't want to scale the reference spectrum to itself, so skip it
			waver = np.linspace(crvalr,crvalr+cdr*(len(datr)-crpixr),len(datr)) # reference spectrum wavelength array
			# retrieve and construct the necessary information about the current spectrum to be scaled
			curstd = stdspectra[angle] # current std spectrum
			print "The current spectrum being scaled is: "+curstd
			hducur = pyfits.open(curstd)
			datc = hducur[1].data.copy() # current std spectrum data
			datc = datc[0][0] # std data is "band 1" (rest is noise, sky, etc.)
			hdrc = hducur[1].header # current std spectrum header
			crpixc = hdrc['CRPIX1'] # current std spectrum WCS information
			cdc = hdrc['CD1_1']
			crvalc = hdrc['CRVAL1']
			# figure out chip gap pixels based on CCDSUM binning header keyword
			ccdsum = (hducur[0].header)['CCDSUM']
			if ccdsum == '2 4': # 2x4 binning
				c1min,c1max = params.chipGapPix24[0]
				c2min,c2max = params.chipGapPix24[1]
			elif ccdsum == '4 4': # 4x4 binning
				c1min,c1max = params.chipGapPix44[0]
				c2min,c2max = params.chipGapPix44[1]
			hducur.close()
			wavec = np.linspace(crvalc,crvalc+cdc*(len(datc)-crpixc),len(datc)) # current std spectrum wavelength array
			# get rid of the chip gaps (based on their ~invariant pixel numbers) for reference spectrum
			pixr = np.arange(len(datr))+1 # 1-based pixel array for reference spectrum
			datrgood = datr[np.logical_or(pixr<=c1min,pixr>=c1max)] # remove first chip gap's flux values
			wavergood = waver[np.logical_or(pixr<=c1min,pixr>=c1max)] # remove first chip gap's wavelengths
			pixrgood = pixr[np.logical_or(pixr<=c1min,pixr>=c1max)] # remove first chip gap's pixels
			datrgood = datrgood[np.logical_or(pixrgood<=c2min,pixrgood>=c2max)] # remove second chip gap's flux values
			wavergood = wavergood[np.logical_or(pixrgood<=c2min,pixrgood>=c2max)] # remove second chip gap's wavelengths
			pixrgood = pixrgood[np.logical_or(pixrgood<=c2min,pixrgood>=c2max)] # remove second chip gap's pixels
			# get rid of the chip gaps (based on their ~invariant pixel numbers) for current spectrum
			pixc = np.arange(len(datc))+1 # 1-based pixel array for current spectrum
			datcgood = datc[np.logical_or(pixc<=c1min,pixc>=c1max)] # remove first chip gap's flux values
			wavecgood = wavec[np.logical_or(pixc<=c1min,pixc>=c1max)] # remove first chip gap's wavelengths
			pixcgood = pixc[np.logical_or(pixc<=c1min,pixc>=c1max)] # remove first chip gap's pixels
			datcgood = datcgood[np.logical_or(pixcgood<=c2min,pixcgood>=c2max)] # remove second chip gap's flux values
			wavecgood = wavecgood[np.logical_or(pixcgood<=c2min,pixcgood>=c2max)] # remove second chip gap's wavelengths
			pixcgood = pixcgood[np.logical_or(pixcgood<=c2min,pixcgood>=c2max)] # remove second chip gap's pixels
			# construct a common wavelength array
			minwavergood = min(wavergood) # minimum wavelength of chip-corrected reference spectrum
			maxwavergood = max(wavergood) # maximum wavelength of chip-corrected ref spectrum
			minwavecgood = min(wavecgood) # minimum wavelength of chip-corrected current spectrum
			maxwavecgood = max(wavecgood) # maximum wavelength of chip-corrected cur spectrum
			# construct common wavelength array based on whether angle or refAngle is lower
			# and then find flux values of reference spectrum at these common wavelengths
			if angle <= refAngle: # the two angles should not really ever be equal (for this step)
				wavecom = wavergood[np.logical_and(wavergood>=minwavergood,wavergood<=maxwavecgood)]
				datrcom = datrgood[np.logical_and(wavergood>=minwavergood,wavergood<=maxwavecgood)]
			elif angle > refAngle:
				wavecom = wavergood[np.logical_and(wavergood>=minwavecgood,wavergood<=maxwavergood)]
				datrcom = datrgood[np.logical_and(wavergood>=minwavecgood,wavergood<=maxwavergood)]
			# find interpolant of entire chip-gap-corrected wavelength values of current spectrum
			fc = interp1d(wavecgood,datcgood)
			# use interpolant to get flux values of current spectrum at common wavelengths only
			datccom = fc(wavecom)
			# use least-squares method to minimize distance between ref and cur spectra at common wavelengths
			print 'A second order scaling function will now be fit to reduce the error between the two spectra at common wavelengths: '
			print 'p[0] + p[1]*x + p[2]*x*x where x is the common wavelength NumPy array.'
			fitfunc = lambda p,x: (p[0] + p[1]*x + p[2]*x*x)*datccom # this order-2 function will scale cur spectrum
			errfunc = lambda p,x,y: fitfunc(p,x) - y
			p0 = [1.0, 0.0, 0.0] # initial guess for parameters of second-order fitfunc is that it is a constant function
			p1,success = leastsq(errfunc,p0,args=(wavecom,datrcom))
			print 'The parameters of the second-order scaling function are ([p0,p1,p2]): '
			print p1
			# apply the scaling function to entire original current spectrum data 
			datcscaled = (p1[0]+p1[1]*wavec+p1[2]*wavec*wavec)*datc
			# copy old file over and replace old data with new scaled data
			suffixlist = curstd.split('.')
			suffix = suffixlist[1][5:]
			namestd = 'std'+angle+'scl'+suffix+'.fits'
			pyrafCalls.run_imcopy(inputname=curstd,outputname=namestd)
			hdunew = pyfits.open(namestd,mode='update')
			datstd = hdunew[1].data.copy()
			# we will replace the 1st spectral band (optimally-extracted spectrum) with the scaled one
			# band 2 will remain as the non-optimal, uncleaned, unweighted, non-scaled spectrum
			# band 3 will remain the non-scaled sky spectrum
			# band 4 will remain the non-scaled sigma spectrum
			datstd[0] = datcscaled
			hdunew[1].data = datstd
			hdunew.flush()
			hdunew.close()
			pipeHistory.updatePipeKeys(inputname=namestd,imagetype='standard',procChar='n')
			# add scaled std spectrum filename to dictionary scaledstandards
			dicts.scaledstandards[angle] = namestd
		# begin second round (r2) of scaling for the extreme gr-angles which may have had insufficient wavelength overlap with refAngle
		rescaledAngles = []
		r2angles = dicts.scaledstandards.keys() # added after the fact: these angles don't contain refAngle due to round 1
		floatAngles = []
		floatAngMaps = {} # will be needed to access keys (string format) of scaledstandards
		for i,v in enumerate(r2angles): # convert angles from strings to floats for numerical comparison
			floatAngles.append(float(v))
			floatAngMaps[float(v)] = v
		refAngleFloat = float(refAngle) # convert refAngle from string to float for numerical comparison
		for r2ang1 in floatAngles:
			if r2ang1 == refAngleFloat:
				continue # this is pairwise rescaling based on the angles with refAngle excluded
			if r2ang1 in rescaledAngles:
				continue # each time an angle is rescaled, it is added to a list so that it is not rescaled again
			# check if r2ang1 spectrum is far enough from refAngle spectrum such that there was insufficient wavelength overlap in round 1
			# a difference in gr-angle of 2 degrees is arbitrary
			# I wonder if it can be made more dynamic or derived through a simple statistical value (median, mean, etc.)
			if abs(r2ang1-refAngleFloat) > 2.0:
				for r2ang2 in floatAngles: # nested for loop to get the other non-refAngle spectrum to scale (or become reference)
					if r2ang2 == refAngleFloat:
						continue # again, don't want to rescale refAngle spectrum. it is the basis for our automatic scaling procedure
					# r2ang2 can be in list rescaledAngles since we may need some of them for reference
					if abs(r2ang2 - r2ang1) < 2.0: # need r2ang1 and r2ang2 spectra to have significant wavelength overlap for meaningful rescaling
						# we want to use the spectrum whose angle is closer to refAngle (preferably median) as the new reference spectrum
						# i.e., we will scale the spectrum with the more distant angle (relative to refAngle) based on the other angle
						distAng1 = abs(r2ang1 - refAngleFloat)
						distAng2 = abs(r2ang2 - refAngleFloat)
						if distAng1 <= distAng2: # if distance between r2ang1 and refAngle is shorter, use r2ang1 as the new refAngle
							newRefAngle = r2ang1
							rescaleAngle = r2ang2
						else:
							newRefAngle = r2ang2
							rescaleAngle = r2ang1
						# begin the rescaling process in the same fashion as the scaling process from round 1
						print 'Beginning RESCALING process for spectrum with GR-ANGLE: '+str(rescaleAngle)
						refstd = dicts.scaledstandards[floatAngMaps[newRefAngle]] # using scaled spectra now for the rescaling process
						curstd = dicts.scaledstandards[floatAngMaps[rescaleAngle]]
						# get data and relevant header information
						hduref = pyfits.open(refstd)
						datr = hduref[1].data.copy() # reference spectrum data
						datr = datr[0][0]
						hdrr = hduref[1].header # reference spectrum header
						crpixr= hdrr['CRPIX1'] # reference spectrum WCS information
						cdr = hdrr['CD1_1']
						crvalr = hdrr['CRVAL1']
						hduref.close()
						waver = np.linspace(crvalr,crvalr+cdr*(len(datr)-crpixr),len(datr)) # reference spectrum wavelength array
						# retrieve and construct the necessary information about the current spectrum to be scaled
						print "The current spectrum being rescaled is: "+curstd
						hducur = pyfits.open(curstd)
						datc = hducur[1].data.copy() # current std spectrum data
						datc = datc[0][0]
						hdrc = hducur[1].header # current std spectrum header
						crpixc = hdrc['CRPIX1'] # current std spectrum WCS information
						cdc = hdrc['CD1_1']
						crvalc = hdrc['CRVAL1']
						# figure out chip gap pixels based on CCDSUM binning header keyword
						ccdsum = (hducur[0].header)['CCDSUM']
						if ccdsum == '2 4': # 2x4 binning
							c1min,c1max = params.chipGapPix24[0]
							c2min,c2max = params.chipGapPix24[1]
						elif ccdsum == '4 4': # 4x4 binning
							c1min,c1max = params.chipGapPix44[0]
							c2min,c2max = params.chipGapPix44[1]
						hducur.close()
						wavec = np.linspace(crvalc,crvalc+cdc*(len(datc)-crpixc),len(datc)) # current std spectrum wavelength array
						# even though these are scaled spectra, the chip gaps are not gone so follow similar procedure as in round 1
						# get rid of the chip gaps (based on their ~invariant pixel numbers) for reference spectrum
						pixr = np.arange(len(datr))+1 # 1-based pixel array for reference spectrum
						datrgood = datr[np.logical_or(pixr<=c1min,pixr>=c1max)] # remove first chip gap's flux values
						wavergood = waver[np.logical_or(pixr<=c1min,pixr>=c1max)] # remove first chip gap's wavelengths
						pixrgood = pixr[np.logical_or(pixr<=c1min,pixr>=c1max)] # remove first chip gap's pixels
						datrgood = datrgood[np.logical_or(pixrgood<=c2min,pixrgood>=c2max)] # remove second chip gap's flux values
						wavergood = wavergood[np.logical_or(pixrgood<=c2min,pixrgood>=c2max)] # remove second chip gap's wavelengths
						pixrgood = pixrgood[np.logical_or(pixrgood<=c2min,pixrgood>=c2max)] # remove second chip gap's pixels
						# get rid of the chip gaps (based on their ~invariant pixel numbers) for current spectrum
						pixc = np.arange(len(datc))+1 # 1-based pixel array for current spectrum
						datcgood = datc[np.logical_or(pixc<=c1min,pixc>=c1max)] # remove first chip gap's flux values
						wavecgood = wavec[np.logical_or(pixc<=c1min,pixc>=c1max)] # remove first chip gap's wavelengths
						pixcgood = pixc[np.logical_or(pixc<=c1min,pixc>=c1max)] # remove first chip gap's pixels
						datcgood = datcgood[np.logical_or(pixcgood<=c2min,pixcgood>=c2max)] # remove second chip gap's flux values
						wavecgood = wavecgood[np.logical_or(pixcgood<=c2min,pixcgood>=c2max)] # remove second chip gap's wavelengths
						pixcgood = pixcgood[np.logical_or(pixcgood<=c2min,pixcgood>=c2max)] # remove second chip gap's pixels
						# construct a common wavelength array
						minwavergood = min(wavergood) # minimum wavelength of chip-corrected reference spectrum
						maxwavergood = max(wavergood) # maximum wavelength of chip-corrected ref spectrum
						minwavecgood = min(wavecgood) # minimum wavelength of chip-corrected current spectrum
						maxwavecgood = max(wavecgood) # maximum wavelength of chip-corrected cur spectrum
						# construct common wavelength array based on whether rescaleAngle or newRefAngle is lower
						# and then find flux values of reference spectrum at these common wavelengths
						if rescaleAngle <= newRefAngle: # the two angles should not really ever be equal (for this step)
							wavecom = wavergood[np.logical_and(wavergood>=minwavergood,wavergood<=maxwavecgood)]
							datrcom = datrgood[np.logical_and(wavergood>=minwavergood,wavergood<=maxwavecgood)]
						elif rescaleAngle > newRefAngle:
							wavecom = wavergood[np.logical_and(wavergood>=minwavecgood,wavergood<=maxwavergood)]
							datrcom = datrgood[np.logical_and(wavergood>=minwavecgood,wavergood<=maxwavergood)]
						# find interpolant of entire chip-gap-corrected wavelength values of current spectrum
						fc = interp1d(wavecgood,datcgood)
						# use interpolant to get flux values of current spectrum at common wavelengths only
						datccom = fc(wavecom)
						# use least-squares method to minimize distance between ref and cur spectra at common wavelengths
						print 'Beginning optimization process for RESCALING.'
						print 'A second order scaling function will now be fit to reduce the error between the two spectra at common wavelengths: '
						print 'p[0] + p[1]*x + p[2]*x*x where x is the common wavelength NumPy array.'
						fitfunc = lambda p,x: (p[0] + p[1]*x + p[2]*x*x)*datccom # this order-2 function will scale cur spectrum
						errfunc = lambda p,x,y: fitfunc(p,x) - y
						p0 = [1.0, 0.0, 0.0] # initial guess for parameters of second-order fitfunc is that it is a constant function
						p1,success = leastsq(errfunc,p0,args=(wavecom,datrcom))
						print 'The parameters of the second-order scaling function are ([p0,p1,p2]): '
						print p1
						# apply the scaling function to entire original (already scaled) current spectrum data 
						datcscaled = (p1[0]+p1[1]*wavec+p1[2]*wavec*wavec)*datc
						# replace old already-scaled data with new RESCALED data
						hdunew = pyfits.open(curstd,mode='update')
						datstd = hdunew[1].data.copy()
						# we will replace the 1st spectral band (optimally-extracted spectrum) with the scaled one
						# band 2 will remain as the non-optimal, uncleaned, unweighted, non-scaled spectrum
						# band 3 will remain the non-scaled sky spectrum
						# band 4 will remain the non-scaled sigma spectrum
						datstd[0] = datcscaled
						hdunew[1].data = datstd
						hdunew.flush()
						# add rescaled spectrum's angle to list rescaledAngles
						rescaledAngles.append(rescaleAngle)
						print 'Rescaling complete for spectrum with GR-ANGLE: '+str(rescaleAngle)
		print "The angles of the rescaled spectra are: "
		print rescaledAngles
	print 'The scaled standard star spectra filenames are: '
	print dicts.scaledstandards
	
	
# This function creates separate sky and sigma spectra from a flux-calibrated science spectrum's 3rd and 4th spectral bands. 
# These new spectra will be combined in combinespectra() to produce combined sky and sigma spectra.
# Since scaleScienceSpectra() doesn't affect bands 3 or 4, this function checks only for the existence of flux-calibrated spectra.
def copySciSkySig(fluxspectra=dicts.fluxsciences,indiv=False):
	# This prints an error message and quits function if inputted dictionaries don't exist.
	try:
		len(fluxspectra)
	except:
		print 'The inputted dictionary does not exist: cannot proceed.'
		return
	# This runs scopy on each flux-calibrated science image to extract bands 3 (sky) and 4 (sigma)
	angles = fluxspectra.keys() # if there are no flux-calibrated spectra, there are no angles, and thus the following loop won't run
	for angle in angles: # creation of individual sky spectra
		fluxsci = fluxspectra[angle]
		suffixlist = fluxsci.split('.')
		suffix = suffixlist[1][5:]
		# need to check for the existence of sky and sigma spectral bands (3 and 4) in case user has custom images
		hduflx = pyfits.open(fluxsci)
		datflx = hduflx[1].data
		if len(datflx) == 4: # by IRAF standards, this must be an array of 4 arrays (these are the bands)
			skyband = '3'
		elif len(datflx) == 5: # why/how? 2014-03-15
			skyband = '4'
		else:
			print "Cannot proceed: the file structure of "+fluxsci+" is not correct."
			print "The data is not of the form: one NumPy array containing four NumPy arrays (the spectral bands)."
			print "Skipping extraction of sky spectral band for: "+fluxsci
			continue
		hduflx.close()
		# define name for sci-sky spectrum and call onedspec.scopy
		namesky = 'sci'+str(angle)+'sky'+suffix+'.fits'
		nameskyAUX = 'sci'+str(angle)+'sky'+suffix+'AUX.fits'
		pyrafCalls.run_scopy(inputname=fluxsci,outputname=nameskyAUX,band=skyband) # band 3 is sky background
		# copy previous spectra to maintain file structure
		pyrafCalls.run_imcopy(inputname=fluxsci,outputname=namesky)
		# add 'RUSKY' keyword to fluxsci 0-header to link it up with the sky spectrum, and get 'BPM' name
		hduflx = pyfits.open(fluxsci,mode='update')
		hdrflx = hduflx[0].header
		hdrflx['RUSKY'] = namesky
		thegain = hdrflx.get('GAIN',2.443) # will be needed in odcombine
		hduflx.flush()
		hdrflx1 = hduflx[1].header
		maskname = hdrflx1.get('BPM','')
		hduflx.close()
		# I must copy the auxiliary data
		hduaux = pyfits.open(nameskyAUX)
		datnew = hduaux[0].data.copy()
		hduaux.close()
		# I must add the data and the 'RUPIPE' and 'RUIMGTYP' and 'BPM' header keywords
		hdusky = pyfits.open(namesky,mode='update')
		hdusky[1].data = datnew
		hdrsky1 = hdusky[1].header
		hdrsky1['BPM'] = maskname
		hdrsky = hdusky[0].header
		hdrsky['GAIN'] = thegain
		hdrsky['RUSKY'] = fluxsci # link this sky spectrum with its flux-calibrated science spectrum
		hdusky.flush()
		hdusky.close()
		pipeHistory.updatePipeKeys(inputname=namesky,imagetype='science',procChar='y')
		os.remove(nameskyAUX)
		# add to skysciences global dictionary and notify user of successful extraction
		dicts.skysciences[angle] = namesky
		print "Sky spectrum successfully extracted: "+namesky+" from "+fluxsci
	# do same thing over gain but now for sigma spectra
	for angle in angles: # creation of individual sigma spectra
		fluxsci = fluxspectra[angle]
		suffixlist = fluxsci.split('.')
		suffix = suffixlist[1][5:]
		# need to check for the existence of sky and sigma spectral bands (3 and 4) in case user has custom images
		hduflx = pyfits.open(fluxsci)
		datflx = hduflx[1].data
		if len(datflx) == 4: # by IRAF standards, this must be an array of 4 arrays (these are the bands)
			sigband = '4'
		elif len(datflx) == 5: # why/how? 2014-03-15
			sigband = '5'
		else:
			print "Cannot proceed: the file structure of "+fluxsci+" is not correct."
			print "The data is not of the form: one NumPy array containing four NumPy arrays (the spectral bands)."
			print "Skipping extraction of sky spectral band for: "+fluxsci
			continue
		hduflx.close()
		# define name for sci-sig spectrum and call onedspec.scopy
		namesig = 'sci'+str(angle)+'sig'+suffix+'.fits'
		namesigAUX = 'sci'+str(angle)+'sig'+suffix+'AUX.fits'
		pyrafCalls.run_scopy(inputname=fluxsci,outputname=namesigAUX,band=sigband) # band 4 is sigma
		# copy previous spectra to maintain file structure
		pyrafCalls.run_imcopy(inputname=fluxsci,outputname=namesig)
		# add 'RUSIG' keyword to fluxsci 0-header to link it up with the sigma spectrum, and get 'BPM' name
		hduflx = pyfits.open(fluxsci,mode='update')
		hdrflx = hduflx[0].header
		hdrflx['RUSIGMA'] = namesig
		thegain = hdrflx.get('GAIN',2.443) # will be needed in odcombine
		hduflx.flush()
		hdrflx1 = hduflx[1].header
		maskname = hdrflx1.get('BPM','')
		hduflx.close()
		# copy auxiliary data
		hduaux = pyfits.open(namesigAUX)
		datnew = hduaux[0].data.copy()
		hduaux.close()
		# I must copy the data and add the 'RUPIPE' and 'RUIMGTYP' and 'BPM' header keywords
		hdusig = pyfits.open(namesig,mode='update')
		hdusig[1].data = datnew
		hdrsig1 = hdusig[1].header
		hdrsig1['BPM'] = maskname
		hdrsig = hdusig[0].header
		hdrsig['GAIN'] = thegain
		hdrsig['RUSIGMA'] = fluxsci # link this sigma spectrum with its flux-calibrated science spectrum
		hdusig.flush()
		hdusig.close()
		pipeHistory.updatePipeKeys(inputname=namesig,imagetype='science',procChar='z')
		os.remove(namesigAUX)
		# add to sigmasciences global dictionary and notify user of successful extraction
		dicts.sigmasciences[angle] = namesig
		print "Sigma spectrum successfully extracted: "+namesig+" from "+fluxsci
	print 'The sky (from science spectra) spectra filenames are:'
	print dicts.skysciences
	print 'The sigma (from science spectra) spectra filenames are:'
	print dicts.sigmasciences
	
	
# This function scales dispersion-corrected science spectra so that they have similar count/sec values at overlapping wavelengths.
# This is for creating the combined count-based spectrum.
# 2014-03-16: changed to scaling counts, not counts/EXPTIME, because that doesn't work if each spectrum has a different exptime.
def scaleDispScienceSpectra(sciencespectra=dicts.dispsciences,indiv=False):
	# This prints an error and quits the function if any of the inputted dictionaries don't exist
	try:
		len(sciencespectra)
	except:
		print 'One of the inputted dictionaries does not exist: cannot proceed.'
		return
	# This begins the count-scaling process for the dispersion-corrected science spectra.
	angles = sciencespectra.keys()
	if len(angles) == 1: # if only one dispersion-corrected spectrum, no other spectra to compare to
		print 'There is only one dispersion-corrected spectrum -- the one with angle: '+angle
		print 'Cannot scale a single spectrum.'
		print 'Skipping flux scaling.'
		return
	elif len(angles) == 2: # only two dispersion-corrected spectra => compare these two only
		print 'There are only two dispersion-corrected spectra to scale.'
		print 'Running this task non-interactively.'
		angle0 = angles[0] # default assumption: 1st element of angles is the lower angle number
		angle1 = angles[1] 
		angOne = float(angles[0]) # angles are strings currently
		angTwo = float(angles[1]) # converting to floats to check the default assumption
		if angTwo < angOne: # need to make this comparison to create appropriate common wavelength array
			angle0 = angles[1]
			angle1 = angles[0]
		sci0 = sciencespectra[angle0] # science spectrum 0 (with lower gr-angle) -- this will be scaled
		sci1 = sciencespectra[angle1] # science spectrum 1 (with greater/equal gr-angle) -- this is the reference
		# retrieving necessary information for science spectrum 0
		hdu0 = pyfits.open(sci0)
		hdr0 = hdu0[1].header
		dat0 = hdu0[1].data.copy()
		dat0 = dat0[0][0] # data for science spectrum itself is the first "band" (rest are noise, sky, etc.)
		crpix0 = hdr0['CRPIX1'] # dispersion (WCS) information for science spectrum 0
		cd0 = hdr0['CD1_1']
		crval0 = hdr0['CRVAL1']
		exptime0 = hdr0['EXPTIME']
		dat0 = dat0 / exptime0 # must scale counts/sec, not counts alone (unlike with flux)
		hdu0.close()
		wave0 = np.linspace(crval0,crval0+cd0*(len(dat0)-crpix0),len(dat0)) # wavelength array for science spectrum 0
		# retrieving necessary information for science spectrum 1
		hdu1 = pyfits.open(sci1)
		hdr1 = hdu1[1].header
		dat1 = hdu1[1].data.copy()
		dat1 = dat1[0][0] # data for science spectrum itself is the first "band" (rest are noise, sky, etc.)
		crpix1 = hdr1['CRPIX1'] # dispersion (WCS) information for science spectrum 1
		cd1 = hdr1['CD1_1']
		crval1 = hdr1['CRVAL1']
		exptime1 = hdr1['EXPTIME']
		#dat1 = dat1 / exptime1 # must scale counts/sec, not counts alone (unlike with flux)
		# figure out chip gap pixels based on CCDSUM binning header keyword
		ccdsum = (hdu1[0].header)['CCDSUM']
		if ccdsum == '2 4': # 2x4 binning
			c1min,c1max = params.chipGapPix24[0]
			c2min,c2max = params.chipGapPix24[1]
		elif ccdsum == '4 4': # 4x4 binning
			c1min,c1max = params.chipGapPix44[0]
			c2min,c2max = params.chipGapPix44[1]
		hdu1.close()
		wave1 = np.linspace(crval1,crval1+cd1*(len(dat1)-crpix1),len(dat1)) # wavelength array for science spectrum 1
		# get rid of the chip gaps in spectrum 0
		pix0 = np.arange(len(dat0))+1 # 1-based pixel array for dat0
		dat0good = dat0[np.logical_or(pix0<=c1min,pix0>=c1max)] # remove count/sec values of first chip gap for spectrum 0
		wave0good = wave0[np.logical_or(pix0<=c1min,pix0>=c1max)] # remove wavelengths of first chip gap for spectrum 0
		pix0good = pix0[np.logical_or(pix0<=c1min,pix0>=c1max)] # remove pixels of first chip gap for spectrum 0
		dat0good = dat0good[np.logical_or(pix0good<=c2min,pix0good>=c2max)] # remove count/sec values of second chip gap for spectrum 0
		wave0good = wave0good[np.logical_or(pix0good<=c2min,pix0good>=c2max)] # remove wavelengths of second chip gap for spectrum 0
		pix0good = pix0good[np.logical_or(pix0good<=c2min,pix0good>=c2max)] # remove pixels of second chip gap for spectrum 0
		# get rid of the chip gaps in spectrum 1
		pix1 = np.arange(len(dat1))+1 # 1-based pixel array for dat1
		dat1good = dat1[np.logical_or(pix1<=c1min,pix1>=c1max)] # remove count/sec values of first chip gap for spectrum 1
		wave1good = wave1[np.logical_or(pix1<=c1min,pix1>=c1max)] # remove wavelengths of first chip gap for spectrum 1
		pix1good = pix1[np.logical_or(pix1<=c1min,pix1>=c1max)] # remove pixels of first chip gap for spectrum 1
		dat1good = dat1good[np.logical_or(pix1good<=c2min,pix1good>=c2max)] # remove count/sec values of second chip gap for spectrum 1
		wave1good = wave1good[np.logical_or(pix1good<=c2min,pix1good>=c2max)] # remove wavelengths of second chip gap for spectrum 1
		pix1good = pix1good[np.logical_or(pix1good<=c2min,pix1good>=c2max)] # remove pixels of second chip gap for spectrum 1
		# construct common wavelength array
		minwave0g = min(wave0good)
		maxwave0g = max(wave0good)
		minwave1g = min(wave1good)
		maxwave1g = max(wave1good)
		wavecom = wave1good[np.logical_and(wave1good>=minwave1g,wave1good<=maxwave0g)] # common wavelength array constructed from wave1
		# find count/sec values of spectrum 1 at common wavelengths
		dat1com = dat1good[np.logical_and(wave1good>=minwave1g,wave1good<=maxwave0g)] # common-wavelength count/sec values of spectrum 1
		# find interpolant of entire chip-gap-corrected wavelength values of spectrum 0
		f0 = interp1d(wave0good,dat0good)
		# use interpolant to get count/sec values of spectrum 0 at common wavelengths only
		dat0com = f0(wavecom)
		# use least-squares method to minimize distance between spectrum 0 and spectrum 1 count/sec values at common wavelengths
		print 'A first order scaling function will now be fit to reduce the error between the two spectra at common wavelengths: '
		print 'p[0] + p[1]*x where x is the common wavelength NumPy array.'
		fitfunc = lambda p,x: (p[0] + p[1]*x)*dat0com # this order-1 function will scale spectrum 0 to spectrum 1 count/sec values
		errfunc = lambda p,x,y: fitfunc(p,x) - y
		p0 = [1.0, 0.0] # initial guess for parameters of first-order fitfunc is a constant function
		p1,success = leastsq(errfunc,p0,args=(wavecom,dat1com))
		print 'The optimal parameters of the first-order scaling function are ([p0,p1]): '
		print p1
		# apply the scaling function to entire original spectrum 0 data 
		dat0scaled = (p1[0]+p1[1]*wave0)*dat0
		#dat0scaled = dat0scaled * exptime0 # reverse mapping from counts/sec to counts
		# copy old file over and replace old data with new scaled data
		suffixlist = sci0.split('.')
		suffix = suffixlist[1][5:]
		namesci = 'sci'+angle0+'cnt'+suffix+'.fits'
		pyrafCalls.run_imcopy(inputname=sci0,outputname=namesci)
		hdu0 = pyfits.open(namesci,mode='update')
		datsci = hdu0[1].data.copy()
		# we will replace the 1st spectral band (optimally-extracted spectrum) with the scaled one
		# band 2 will remain as the non-optimal, uncleaned, unweighted, non-scaled spectrum
		# band 3 will remain the non-scaled sky spectrum
		# band 4 will remain the non-scaled sigma spectrum
		datsci[0] = dat0scaled
		hdu0[1].data = datsci
		hdu0.flush()
		hdu0.close()
		pipeHistory.updatePipeKeys(inputname=namesci,imagetype='science',procChar='c')
		# add scaled dispersion-corrected science spectrum filename to dictionary scaleddispsciences
		dicts.scaleddispsciences[angle0] = namesci
	elif len(angles) > 2: # for multiple dispersion-corrected spectra, ask the user for a gr-angle and/or do a pairwise loop
		print "There are more than two dispersion-corrected spectra."
		print "You have the option of picking one spectrum's GR-ANGLE."
		print "The count/sec values at overlapping wavelengths for all other spectra will be scaled to this spectrum."
		print "Here is the global dictionary of dispersion-corrected science images:"
		print sciencespectra
		# for a non-interactive version, the median angle could just be chosen automatically
		# for an even number of spectra, the higher of the two median angles could be chosen (more counts than bluer angles)
		while True:
			refAngle = raw_input("Please enter a (preferably MEDIAN) GR-ANGLE: ")
			if refAngle in angles: 
				print "You chose to scale all other spectra with reference to: "+sciencespectra[refAngle]
				break
			else:
				print "Invalid input: you must enter an angle in the above dictionary to proceed."
		# all other spectra must now be scaled based on sciencespectra[refAngle]
		# all spectra should have overlapping wavelengths with each other
		# it is this common wavelength range, different for each pair, that will be used to find the optimal scaling
		# a first-order polynomial/function will be fit to each of the spectra that need to be scaled
		# first, gather and construct the necessary information about the reference spectrum
		refsci = sciencespectra[refAngle]
		hduref = pyfits.open(refsci)
		datr = hduref[1].data.copy() # reference spectrum data
		datr = datr[0][0] # science data is in "band 1" (rest is noise, sky, etc.)
		hdrr = hduref[1].header # reference spectrum header
		crpixr= hdrr['CRPIX1'] # reference spectrum WCS information
		cdr = hdrr['CD1_1']
		crvalr = hdrr['CRVAL1']
		exptimer = hdrr['EXPTIME']
		#datr = datr / exptimer # must scale counts/sec, not counts alone (unlike with flux)
		hduref.close()
		# next, use a for loop to create a scaled spectrum for each non-refsci spectrum
		for angle in angles:
			if angle == refAngle: 
				continue # we don't want to scale the reference spectrum to itself, so skip it
			waver = np.linspace(crvalr,crvalr+cdr*(len(datr)-crpixr),len(datr)) # reference spectrum wavelength array
			# retrieve and construct the necessary information about the current spectrum to be scaled
			cursci = sciencespectra[angle] # current science spectrum
			print "The current spectrum being scaled is: "+cursci
			hducur = pyfits.open(cursci)
			datc = hducur[1].data.copy() # current science spectrum data
			datc = datc[0][0] # science data is "band 1" (rest is noise, sky, etc.)
			hdrc = hducur[1].header # current science spectrum header
			crpixc = hdrc['CRPIX1'] # current science spectrum WCS information
			cdc = hdrc['CD1_1']
			crvalc = hdrc['CRVAL1']
			exptimec = hdrc['EXPTIME']
			#datc = datc / exptimec # must scale counts/sec, not counts alone (unlike with flux)
			# figure out chip gap pixels based on CCDSUM binning header keyword
			ccdsum = (hducur[0].header)['CCDSUM']
			if ccdsum == '2 4': # 2x4 binning
				c1min,c1max = params.chipGapPix24[0]
				c2min,c2max = params.chipGapPix24[1]
			elif ccdsum == '4 4': # 4x4 binning
				c1min,c1max = params.chipGapPix44[0]
				c2min,c2max = params.chipGapPix44[1]
			hducur.close()
			wavec = np.linspace(crvalc,crvalc+cdc*(len(datc)-crpixc),len(datc)) # current science spectrum wavelength array
			# get rid of the chip gaps (based on their ~invariant pixel numbers) for reference spectrum
			pixr = np.arange(len(datr))+1 # 1-based pixel array for reference spectrum
			datrgood = datr[np.logical_or(pixr<=c1min,pixr>=c1max)] # remove first chip gap's count/sec values
			wavergood = waver[np.logical_or(pixr<=c1min,pixr>=c1max)] # remove first chip gap's wavelengths
			pixrgood = pixr[np.logical_or(pixr<=c1min,pixr>=c1max)] # remove first chip gap's pixels
			datrgood = datrgood[np.logical_or(pixrgood<=c2min,pixrgood>=c2max)] # remove second chip gap's count/sec values
			wavergood = wavergood[np.logical_or(pixrgood<=c2min,pixrgood>=c2max)] # remove second chip gap's wavelengths
			pixrgood = pixrgood[np.logical_or(pixrgood<=c2min,pixrgood>=c2max)] # remove second chip gap's pixels
			# get rid of the chip gaps (based on their ~invariant pixel numbers) for current spectrum
			pixc = np.arange(len(datc))+1 # 1-based pixel array for current spectrum
			datcgood = datc[np.logical_or(pixc<=c1min,pixc>=c1max)] # remove first chip gap's count/sec values
			wavecgood = wavec[np.logical_or(pixc<=c1min,pixc>=c1max)] # remove first chip gap's wavelengths
			pixcgood = pixc[np.logical_or(pixc<=c1min,pixc>=c1max)] # remove first chip gap's pixels
			datcgood = datcgood[np.logical_or(pixcgood<=c2min,pixcgood>=c2max)] # remove second chip gap's count/sec values
			wavecgood = wavecgood[np.logical_or(pixcgood<=c2min,pixcgood>=c2max)] # remove second chip gap's wavelengths
			pixcgood = pixcgood[np.logical_or(pixcgood<=c2min,pixcgood>=c2max)] # remove second chip gap's pixels
			# construct a common wavelength array
			minwavergood = min(wavergood) # minimum wavelength of chip-corrected reference spectrum
			maxwavergood = max(wavergood) # maximum wavelength of chip-corrected ref spectrum
			minwavecgood = min(wavecgood) # minimum wavelength of chip-corrected current spectrum
			maxwavecgood = max(wavecgood) # maximum wavelength of chip-corrected cur spectrum
			# construct common wavelength array based on whether angle or refAngle is lower
			# and then find count/sec values of reference spectrum at these common wavelengths
			if angle <= refAngle: # the two angles should not really ever be equal (for this step)
				wavecom = wavergood[np.logical_and(wavergood>=minwavergood,wavergood<=maxwavecgood)]
				datrcom = datrgood[np.logical_and(wavergood>=minwavergood,wavergood<=maxwavecgood)]
			elif angle > refAngle:
				wavecom = wavergood[np.logical_and(wavergood>=minwavecgood,wavergood<=maxwavergood)]
				datrcom = datrgood[np.logical_and(wavergood>=minwavecgood,wavergood<=maxwavergood)]
			# find interpolant of entire chip-gap-corrected wavelength values of current spectrum
			fc = interp1d(wavecgood,datcgood)
			# use interpolant to get count/sec values of current spectrum at common wavelengths only
			datccom = fc(wavecom)
			# use least-squares method to minimize distance between ref and cur spectra at common wavelengths
			print 'A first order scaling function will now be fit to reduce the error between the two spectra at common wavelengths: '
			print 'p[0] + p[1]*x where x is the common wavelength NumPy array.'
			fitfunc = lambda p,x: (p[0] + p[1]*x)*datccom # this order-1 function will scale cur spectrum
			errfunc = lambda p,x,y: fitfunc(p,x) - y
			p0 = [1.0, 0.0] # initial guess for parameters of first-order fitfunc is that it is a constant function
			p1,success = leastsq(errfunc,p0,args=(wavecom,datrcom))
			print 'The parameters of the first-order scaling function are ([p0,p1]): '
			print p1
			# apply the scaling function to entire original current spectrum data 
			datcscaled = (p1[0] + p1[1]*wavec)*datc
			#datcscaled = datcscaled * exptimec # reverse mapping from counts/sec to counts
			# copy old file over and replace old data with new scaled data
			suffixlist = cursci.split('.')
			suffix = suffixlist[1][5:]
			namesci = 'sci'+angle+'cnt'+suffix+'.fits'
			pyrafCalls.run_imcopy(inputname=cursci,outputname=namesci)
			hdunew = pyfits.open(namesci,mode='update')
			datsci = hdunew[1].data.copy()
			# we will replace the 1st spectral band (optimally-extracted spectrum) with the scaled one
			# band 2 will remain as the non-optimal, uncleaned, unweighted, non-scaled spectrum
			# band 3 will remain the non-scaled sky spectrum
			# band 4 will remain the non-scaled sigma spectrum
			datsci[0] = datcscaled
			hdunew[1].data = datsci
			hdunew.flush()
			hdunew.close()
			pipeHistory.updatePipeKeys(inputname=namesci,imagetype='science',procChar='c')
			# add scaled disp science spectrum filename to dictionary scaleddispsciences
			dicts.scaleddispsciences[angle] = namesci
	print 'The scaled dispersion-corrected science spectra filenames are: '
	print dicts.scaleddispsciences
	
	
# This function runs odcombine on the dispersion-corrected and flux-calibrated (or scaled) science spectra.
# It will also combine the sky and sigma spectra, and flux-calibrated standard star spectra.
def combinespectra(dispspectra=dicts.dispsciences,fluxspectra=dicts.fluxsciences,scaledfluxspectra=dicts.scaledfluxsciences,scaleddispspectra=dicts.scaleddispsciences,skyspectra=dicts.skysciences,sigmaspectra=dicts.sigmasciences,fluxstd=dicts.fluxstandards,sclstd=dicts.scaledstandards,indiv=False):
	# This prints an error message and quits function if inputted dictionaries don't exist.
	try:
		len(dispspectra),len(fluxspectra),len(scaledfluxspectra),len(scaleddispspectra)
	except:
		print 'An inputted dictionary does not exist: cannot proceed.'
		return
	# this combines the dispersion-corrected science spectra
	namescidisp = 'scifinal_dsp.fits'
	dispspectratocombine = []
	anglesdsp = dispspectra.keys() # if no disp spectra, then no angles, so the following would never be run
	for angle in anglesdsp:
		nextSpectrum = dispspectra[angle]+'[1]'
		scaledsci = scaleddispspectra.get(angle,'') # use scaled disp science spectrum instead if one is found
		if scaledsci != '':
			nextSpectrum = scaledsci+'[1]'
		dispspectratocombine.append(nextSpectrum)
	try:
		if len(dispspectratocombine)>1:
			img = dispspectratocombine[0][:-3] # gets an image filename without [1] to find 'GAIN'
			hdulist = pyfits.open(img)
			hdr0 = hdulist[0].header
			ccdgain = hdr0.get('GAIN',2.443) # reasonable assumption: all science spectra/images have same gain
			hdulist.close()
			odcombinedispseq = ','.join(dispspectratocombine)
			pyrafCalls.run_odcombine(inputlistseq=odcombinedispseq,finalspectrumname=namescidisp,saltgain=ccdgain,customRun=indiv)
			print 'The combined dispersion-corrected science image filename is:'
			print namescidisp
			# add 'RUPIPE' and 'RUIMGTYP' keywords to namesciflux's 0-header
			pipeHistory.updatePipeKeys(inputname=namescidisp,imagetype='science',procChar='u')
			# add dispersion-corrected spectrum filename to index 'dsp' of dictionary combinedspectra
			dicts.combinedspectra['dsp'] = namescidisp
		else:
			print 'There are not enough dispersion-corrected spectra: cannot run odcombine.'
	except:
		print 'Error running odcombine for dispersion-corrected science spectra: did not proceed.'
	# this combines the (scaled) flux-calibrated science spectra (if they exist)
	namesciflux = 'scifinal_flx.fits'
	fluxspectratocombine = []
	anglesflx = fluxspectra.keys() # if no flux spectra, then no angles, so the following would never be run
	for angle in anglesflx:
		nextSpectrum = fluxspectra[angle]+'[1]'
		scaledsci = scaledfluxspectra.get(angle,'') # use scaled science spectrum instead if one is found
		if scaledsci != '':
			nextSpectrum = scaledsci+'[1]'
		fluxspectratocombine.append(nextSpectrum)
	try:
		if len(fluxspectratocombine)>1:
			img = fluxspectratocombine[0][:-3] # gets an image filename without [1] to find 'GAIN'
			hdulist = pyfits.open(img)
			hdr0 = hdulist[0].header
			ccdgain = hdr0.get('GAIN',2.443) # reasonable assumption: all science spectra/images have same gain
			hdulist.close()
			odcombinefluxseq = ','.join(fluxspectratocombine)
			pyrafCalls.run_odcombine(inputlistseq=odcombinefluxseq,finalspectrumname=namesciflux,saltgain=ccdgain,customRun=indiv)
			print 'The combined flux-calibrated science image filename is:'
			print namesciflux
			# add 'RUPIPE' and 'RUIMGTYP' keywords to namesciflux's 0-header
			pipeHistory.updatePipeKeys(inputname=namesciflux,imagetype='science',procChar='v')
			# add flux-calibrated spectrum filename to index 'flx' of dictionary combinedspectra
			dicts.combinedspectra['flx'] = namesciflux
		else:
			print 'There are not enough flux-calibrated spectra: cannot run odcombine.'
	except:
		print 'Error running odcombine for flux-calibrated science spectra: did not proceed.'
	# if dispersion-corrected (and maybe flux-calibrated) combinations were done, also do the sky and sigma spectra combinations
	# later, add them as 2nd, 3rd extensions to scifinal_dsp.fits and scifinal_flx.fits
	if dicts.combinedspectra.get('dsp','') != '':
		anglessky = skyspectra.keys()
		skyspectratocombine = []
		for angle in anglessky:
			nextSkySpectrum = skyspectra[angle] # there are no 'SCI' FITS extensions in the sky spectra
			skyspectratocombine.append(nextSkySpectrum)
		if len(skyspectratocombine) > 1:
			# get gain from one sky spectrum's 0-header (reasonable assumption: all sky spectra have same gain)
			img = skyspectratocombine[0]
			hdulist = pyfits.open(img)
			hdr0 = hdulist[0].header
			ccdgain = hdr0.get('GAIN',2.443) # reasonable assumption: all science spectra/images have same gain
			hdulist.close()
			for i,v in enumerate(skyspectratocombine):
				skyspectratocombine[i] = v+'[1]' # data lies in 1-extension
			odcombineskyseq = ','.join(skyspectratocombine)
			namesky = 'scifinal_sky.fits'
			pyrafCalls.run_odcombine(inputlistseq=odcombineskyseq,finalspectrumname=namesky,saltgain=ccdgain,customRun=indiv)
			# add header keywords to scifinal_sky.fits
			pipeHistory.updatePipeKeys(inputname=namesky,imagetype='science',procChar='w')
			# add sky spectrum filename to index 'sky' of dictionary combinedspectra
			dicts.combinedspectra['sky'] = namesky
		# do same thing over again for sigma spectra
		anglessig = sigmaspectra.keys()
		sigspectratocombine = []
		for angle in anglessig:
			nextSigSpectrum = sigmaspectra[angle] # there are no 'SCI' FITS extensions in the sig spectra
			sigspectratocombine.append(nextSigSpectrum)
		if len(sigspectratocombine) > 1:
			# get gain from one sig spectrum's 0-header (reasonable assumption: all sig spectra have same gain)
			img = sigspectratocombine[0]
			hdulist = pyfits.open(img)
			hdr0 = hdulist[0].header
			ccdgain = hdr0.get('GAIN',2.443) # reasonable assumption: all science spectra/images have same gain
			hdulist.close()
			for i,v in enumerate(sigspectratocombine):
				sigspectratocombine[i] = v+'[1]' # data lies in 1-extension
			odcombinesigseq = ','.join(sigspectratocombine)
			namesig = 'scifinal_sig.fits'
			pyrafCalls.run_odcombine(inputlistseq=odcombinesigseq,finalspectrumname=namesig,saltgain=ccdgain,customRun=indiv)
			# add header keywords to scifinal_sig.fits
			pipeHistory.updatePipeKeys(inputname=namesig,imagetype='science',procChar='x')
			# add sig spectrum filename to index 'sig' of dictionary combinedspectra
			dicts.combinedspectra['sig'] = namesig
	# this combines the (scaled) flux-calibrated standard star spectra (if they exist)
	namestdflux = 'stdfinal_flx.fits'
	fluxspectratocombine = []
	anglesflx = fluxstd.keys() # if no flux spectra, then no angles, so the following would never be run
	for angle in anglesflx:
		nextSpectrum = fluxstd[angle]+'[1]'
		scaledstd = sclstd.get(angle,'') # use scaled science spectrum instead if one is found
		if scaledstd != '':
			nextSpectrum = scaledstd+'[1]'
		fluxspectratocombine.append(nextSpectrum)
	try:
		if len(fluxspectratocombine)>1:
			img = fluxspectratocombine[0][:-3] # gets an image filename without [1] to find 'GAIN'
			hdulist = pyfits.open(img)
			hdr0 = hdulist[0].header
			ccdgain = hdr0.get('GAIN',2.443) # reasonable assumption: all science spectra/images have same gain
			hdulist.close()
			odcombinefluxseq = ','.join(fluxspectratocombine)
			pyrafCalls.run_odcombine(inputlistseq=odcombinefluxseq,finalspectrumname=namestdflux,saltgain=ccdgain,customRun=indiv)
			print 'The combined flux-calibrated science image filename is:'
			print namestdflux
			# add 'RUPIPE' and 'RUIMGTYP' keywords to namestdflux's 0-header
			pipeHistory.updatePipeKeys(inputname=namestdflux,imagetype='standard',procChar='v')
			# add flux-calibrated spectrum filename to index 'flx' of dictionary combinedspectra
			dicts.combinedspectra['std'] = namestdflux
		else:
			print 'There are not enough flux-calibrated standard star spectra: cannot run odcombine.'
	except:
		print 'Error running odcombine for flux-calibrated standard star spectra: did not proceed.'
	
			
# This function operates on the combined flux-calibrated (and/or dispersion-corrected) science spectra.
# It interpolates over remaining chip gaps (that odcombine could not mask with non-existent replacement data)
# and interpolates over artificial deviations from the spectrum (background over-subtraction or cosmic rays).
def cleanCombinedSpectra(combspectra=dicts.combinedspectra,indiv=False):
	# This prints an error message and quits function if inputted dictionaries don't exist.
	try:
		len(combspectra)
	except:
		print 'The inputted dictionary does not exist: cannot proceed.'
		return
	uncleanSpectra = []
	if combspectra.get('flx','') != '':
		uncleanSpectra.append(combspectra.get('flx',''))
	if combspectra.get('std','') != '':
		uncleanSpectra.append(combspectra.get('std',''))
	for s in uncleanSpectra:
		# This reads in the FITS file header and data.
		hdu = pyfits.open(s)
		dat = hdu[0].data.copy() # combined spectra have only one FITS extension (0-extension)
		hdr = hdu[0].header
		cd = hdr['CD1_1'] # Dispersion (wavelength) information.
		crval = hdr['CRVAL1'] # Dispersion (wavelength information).
		crpix = hdr['CRPIX1']
		wave = np.linspace(crval,crval+cd*(len(dat)-crpix),len(dat)) # Wavelength array same size as dat.
		# 2 3 4 3 4 3 4 5 6 0 5 4 3 2 4 
		# Split dat and wave into subarrays -- each of which will be cleaned of leftover chip gaps.
		n = 10 # Number of subarrays. Can be made fancier for special cases.
		dats = np.array_split(dat,n) # NumPy array of n arrays containing the data.
		waves = np.array_split(wave,n) # NumPy array of n arrays containing the wavelengths.

		datFinal = [] # Final data array will be constructed by appending processed data subarrays to it.

		for i,d in enumerate(dats): # Filter each data subarray. 
			w = waves[i] # Corresponding wavelength subarray.
			goodInd = (np.where(d != 0))[0] # Indices where data subarray is not 0 (non-chip gap). 
			dGood = d[goodInd] # Non-chip gap data.
			wGood = w[goodInd] # Corresponding non-chip gap wavelengths.
			badInd = (np.where(d == 0))[0] # Indices where data subarray is 0 (chip gap).
			wBad = w[badInd] # Chip gap wavelengths for interpolation.
			print(w.min(),w.max(),len(dGood),len(wGood),len(wBad),badInd)
			dnew = d.copy()
			if len(wGood) != 0 and len(dGood) != 0: # If no chip gaps, don't try to interpolate.	
				f = interp1d(wGood,dGood)
				# v = np.empty_like(chipInd) # Corresponding interpolated values at chip gap indices.
				for p,q in enumerate(badInd): # Replace old values with interpolated values.
					try:
						v = f(wBad[p]) # Assuming one-to-one mapping is maintained between chipInd, wBad, and v elements.
					except: # In case bad wavelength is outside the good wavelength range (can't interpolate).
						v = np.median(dGood)
					dnew[q] = v
			datFinal = np.append(datFinal,dnew) # Append dnew to end of final data array.

		# Phase 2

		win = np.ones(5)
		datSmooth = np.convolve(win/win.sum(),datFinal,mode='same')
		datDiv = np.abs(datSmooth / datFinal)
		datBad = datDiv / np.std(datDiv) > 5
		(np.where(datBad == True)) [0]
		badBadPix = (np.where(datBad == True)) [0]
		datFinal2 = datFinal.copy()
		for i in badBadPix:
			datFinal2[i] = datSmooth[i]
		
		if s[0:3] == 'sci':
			namecln = 'scifinal_cln.fits'
			imgtype = 'science'
		elif s[0:3] == 'std':
			namecln = 'stdfinal_cln.fits'
			imgtype = 'standard'
		else:
			print "Cannot continue: could not classify spectrum as science or standard star."
			continue
		pyrafCalls.run_imcopy(inputname=s,outputname=namecln,numExt=1)
		hdunew = pyfits.open(namecln,mode='update')
		hdunew[0].data = datFinal2
		hdunew.flush()
		hdunew.close()
		pipeHistory.updatePipeKeys(inputname=namecln,imagetype=imgtype,procChar='k')
		if imgtype == 'science':
			dicts.combinedspectra['clnsci'] = namecln
		elif imgtype == 'standard':
			dicts.combinedspectra['clnstd'] = namecln
	

# This function telluric-corrects the combined, cleaned science spectrum.
def telluricCorrect(combspectra=dicts.combinedspectra,indiv=False):
	# This prints an error message and quits function if inputted dictionaries don't exist.
	try:
		len(combspectra)
	except:
		print 'The inputted dictionary does not exist: cannot proceed.'
		return
	# Importing science and standard star header info and data from FITS files.
	if combspectra.get('clnsci','')  == '' or combspectra.get('flx','') == '':
		print "Cannot proceed: no cleaned or combined flux-cal science spectrum found."
		return
	else:
		sci = combspectra.get('clnsci','')
	if combspectra.get('clnstd','')  == '' or combspectra.get('std','') == '':
		print "Cannot proceed: no cleaned or combined flux-cal standard star spectrum found."
		return
	else:
		std = combspectra.get('clnstd','')
	hdusci = pyfits.open(sci)
	hdrsci = hdusci[0].header
	datsciOrig = hdusci[0].data.copy()

	hdustd = pyfits.open(std) 
	hdrstd = hdustd[0].header
	datstdOrig = hdustd[0].data.copy()

	# Standardize data so that meaningful comparisons can be made.
	datsci = (datsciOrig-np.mean(datsciOrig))/np.std(datsciOrig)
	datstd = (datstdOrig-np.mean(datstdOrig))/np.std(datstdOrig)

	# Use WCS information from headers to create wavelength arrays.
	crpixstd = hdrstd['CRPIX1']
	cdstd = hdrstd['CD1_1']
	crvalstd = hdrstd['CRVAL1']
	wavestd = np.linspace(crvalstd,crvalstd+cdstd*(len(datstd)-crpixstd),len(datstd))

	crpixsci = hdrsci['CRPIX1']
	cdsci = hdrsci['CD1_1']
	crvalsci = hdrsci['CRVAL1']
	wavesci = np.linspace(crvalsci,crvalsci+cdsci*(len(datsci)-crpixsci),len(datsci))

	# Define and combine telluric band wavelength and interpolation ranges here.
	beginA = 7590 # Fraunhofer telluric A band (O_2)
	endA = 7685
	beginB = 6855 # Fraunhofer telluric B band (O_2)
	endB = 6935
	beginC = 6270 # Fraunhofer telluric a band (O_2)
	endC = 6310 # The a-band is not always present or properly resolved in the spectra => sometimes meaningless correction

	beginAinterp = 7500
	endAinterp = 7750
	beginBinterp = 6800
	endBinterp = 7000
	beginCinterp = 6200
	endCinterp = 6400

	telluricLines = [(beginA,endA,beginAinterp,endAinterp),(beginB,endB,beginBinterp,endBinterp)] # not doing a-band

	# Begin telluric correction process.
	corsci = datsci.copy()
	for beginT,endT,beginTinterp,endTinterp in telluricLines:
		wavestdTOrig = wavestd[np.where(np.logical_and(wavestd>beginT,wavestd<endT))] # For interactive comparison purposes.
		datstdTOrig = datstd[np.where(np.logical_and(wavestd>beginT,wavestd<endT))] # For interactive comparison purposes.
		print("Working in telluric region... "+str(beginT)+":"+str(endT)+" Angstroms.")
		# Find wavelengths and data subarray in extended telluric band surrounding.
		wavesciTFull = wavesci[np.where(np.logical_and(wavesci>beginTinterp,wavesci<endTinterp))] # Wavelengths of extended (surrounding) T band region.
		datsciTFull = datsci[np.where(np.logical_and(wavesci>beginTinterp,wavesci<endTinterp))] # Data of extended T band region.
		wavestdTFull = wavestd[np.where(np.logical_and(wavestd>beginTinterp,wavestd<endTinterp))] # Wavelengths of extended (surrounding) T band region.
		datstdTFull = datstd[np.where(np.logical_and(wavestd>beginTinterp,wavestd<endTinterp))] # Data of extended T band region.
		# Cross-correlate both sub-spectra to find wavelength shift of standard star.
		corrT = np.correlate(datstdTFull,datsciTFull,mode='full')
		xcorrT = np.arange(len(corrT))
		lagsT = xcorrT - len(datstdTFull) - 1
		dPerLagT = (wavestdTFull[-1] - wavestdTFull[0])/float(len(wavestdTFull))
		offsetsT = -lagsT * dPerLagT
		shiftT = offsetsT[np.argmax(corrT)]
		print("Standard star shifted by "+str(shiftT)+ " Angstroms.")
		# Apply the wavelength shift to the entire standard star spectrum and get the new data subarray.
		wavestdS = wavestd+shiftT
		# datstdTShifted = datstd[np.where(np.logical_and(wavestd>beginT,wavestd<endT))]
		# datstdT = datstd[np.where(np.logical_and(wavestdTemp>beginT,wavestdTemp<endT))]
		# Find wavelengths and data subarray in current telluric band.
		wavesciT = wavesci[np.where(np.logical_and(wavesci>beginT,wavesci<endT))]
		datsciT = datsci[np.where(np.logical_and(wavesci>beginT,wavesci<endT))]
		wavestdT = wavestdS[np.where(np.logical_and(wavestdS>beginT,wavestdS<endT))]
		datstdT = datstd[np.where(np.logical_and(wavestdS>beginT,wavestdS<endT))]
		# Error condition in case two data subarrays are not equal.
		if len(datstdT) > len(datsciT):
			print("Correcting for len(datSTD) > len(datasci) by "+str(np.abs(len(datstdT)-len(datsciT)))+"px")
			for i in np.arange(np.abs(len(datstdT)-len(datsciT)))+1: # remove last element
				wavestdT = wavestdT[0:-i]
				datstdT = datstdT[0:-i]
		elif len(datstdT) < len(datsciT):
			print("Correcting for len(datSTD) < len(datasci) by "+str(np.abs(len(datstdT)-len(datsciT)))+"px")
			for i in np.arange(np.abs(len(datstdT)-len(datsciT)))+1: # add element from shifted std star arrays
				maxInd = (np.where(np.logical_and(wavestdS>beginT,wavestdS<endT)))[0].max()
				wavestdT = np.append(wavestdT,wavestdS[maxInd+i]) # extend the smaller array to match the larger
				datstdT = np.append(datstdT,datstd[maxInd+i])
		# Scale, shift, and stretch standard star sub-spectrum to match science sub-spectrum.
		fitfunc = lambda p,x: (p[0]+p[1]*x)*datstdT
		errfunc = lambda p,x,y: fitfunc(p,x) - y
		p0 = [1.0,0.0]
		p1,success = leastsq(errfunc,p0,args=(wavesciT,datsciT))
		print("Standard star polynomial scaling parameters:")
		print(p1)
		datstdT = (p1[0] + p1[1]*wavestdT)*datstdT
		# Divide science sub-spectrum by transformed standard star sub-spectrum to get correction.
		corT = np.divide(datsciT,datstdT)
		# Find interpolation over telluric band to scale correction values to continuum flux values.
		wavesciTFull = wavesci[np.where(np.logical_and(wavesci>beginTinterp,wavesci<endTinterp))] # Wavelengths of extended (surrounding) T band region.
		datsciTFull = datsci[np.where(np.logical_and(wavesci>beginTinterp,wavesci<endTinterp))] # Data of extended T band region.
		wavesciTNonTell = wavesciTFull[np.where(np.logical_or(wavesciTFull<beginT,wavesciTFull>endT))] # Wavelengths of extended non-telluric T region (used to find interpolant).
		datsciTNonTell = datsciTFull[np.where(np.logical_or(wavesciTFull<beginT,wavesciTFull>endT))] # Data of extended non-telluric T region (used to find interpolant).
		fT = interp1d(wavesciTNonTell,datsciTNonTell) # Interpolating function in the extended non-telluric T region.
		intT = fT(wavesciT)
		# Scale corrected values to interpolated values to maintain continuum level.
		fitfunc = lambda p,x: (p[0]+p[1]*x+p[2]*x*x)*corT
		errfunc = lambda p,x,y: fitfunc(p,x) - y
		p0 = [1.0,0.0,0.0]
		p1,success = leastsq(errfunc,p0,args=(wavesciT,intT))
		print("Correction polynomial scaling parameters:")
		print p1
		corsclT = (p1[0]+p1[1]*wavesciT+p1[2]*wavesciT*wavesciT)*corT
		# Replace telluric band in entire spectrum with the corrected and scaled data values.
		corsci[np.where(np.logical_and(wavesci>beginT,wavesci<endT))] = corsclT.astype(corsci.dtype)
	
		print ""


	# De-regularize entire corrected science spectrum.
	sciFinal = corsci*np.std(datsciOrig)+np.mean(datsciOrig)
	# Copy telluric-corrected data to new file.
	nametel = 'scifinal_tel.fits'
	pyrafCalls.run_imcopy(inputname=sci,outputname=nametel,numExt=1)
	hdunew = pyfits.open(nametel,mode='update')
	hdunew[0].data = sciFinal
	hdunew.flush()
	hdunew.close()
	pipeHistory.updatePipeKeys(inputname=nametel,imagetype='science',procChar='t')
	dicts.combinedspectra['tel'] = nametel
	print "The name of the telluric-corrected science spectrum is:"
	print nametel

