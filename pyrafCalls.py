<<<<<<< HEAD
'''
>>> Rutgers SALT Supernova Spectral Reduction Pipeline <<<

This module contains functions which run a specific PyRAF task
corresponding to one step of the entire reduction process.

Each function sets the parameters for a PyRAF task and then
runs it on the data. These functions are primarily called by 
the core module containing the main pipeline functions: tasks.py.

Please refer to the documentation for more information about
how each PyRAF call affects the data.

*** Modifications ***
Sept. 26, 2013: Created module. -Viraj Pandya

'''

import pyfits

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

import dicts # These are the pipeline's global dictionaries.
import params # Customizable parameters for the pipeline.

# This function sets the parameters for the pyraf task ccdred.flatcombine and then runs it.
def run_flatcombine(flatstocombine,combflatname,customRun=False):
	# This clears the current parameter values for flatcombine, and sets general parameters again.
	ccdred.flatcombine.unlearn()
	ccdred.flatcombine.combine='average'
	ccdred.flatcombine.reject='avsigclip'
	ccdred.flatcombine.ccdtype=''
	ccdred.flatcombine.scale='mode'
	ccdred.flatcombine.process='no'
	ccdred.flatcombine.subsets='no'
	
	# appends '[1]' to end of each filename since that's the extension for the data to be combined
	flatstofeed = list(flatstocombine) # making copy of list (will need unmodified original for imcopy)
	for i,v in enumerate(flatstofeed):
		flatstofeed[i] = v+'[1]'
	
	# turn python list (of string filenames) into a string sequence
	flatcombineseq = ','.join(flatstofeed)
	# If running task individually, pop up epar window to let users change parameters if they want.
	if customRun == True:
		while True: 
			eparAnswer = raw_input("Do you want to further edit the parameters? 0 for no, 1 for yes.")
			if eparAnswer == '0' or eparAnswer == '1':
				break
			else:
				print "Invalid input: you must enter either 0 (no) or 1 (yes)."		
		if eparAnswer == '1':
			iraf.eparam(ccdred.flatcombine)
	
	ccdred.flatcombine(input=flatcombineseq,output=combflatname)

# This function sets the parameters for the pyraf task longslit.response and then runs it.
def run_response(combinedflat,normflatname,customRun=False):
	# this clears current parameter values for longslit.response and sets general parameters.
	longslit.response.unlearn()
	longslit.response.threshold=5.0
	longslit.response.function='spline3'
	longslit.response.order=19
	longslit.response.low_rej=3.0
	longslit.response.high_rej=3.0
	longslit.response.niterate=3
	longslit.response.grow=2
	longslit.response.interactive='no' # 2013-12-06: prevents bug in iraf that auto-inputs dispaxis=2
	# If running task individually, pop up epar window to let users change parameters if they want.
	if customRun == True:
		while True: 
			eparAnswer = raw_input("Do you want to further edit the parameters? 0 for no, 1 for yes.")
			if eparAnswer == '0' or eparAnswer == '1':
				break
			else:
				print "Invalid input: you must enter either 0 (no) or 1 (yes)."		
		if eparAnswer == '1':
			iraf.eparam(longslit.response)
	# This makes sure that the DISPAXIS is set correctly (if 'DISPAXIS' keyword not in header)
	longslit.dispaxis = 1
	# This runs longslit.response (interactively for now)
	longslit.response(calibration=combinedflat,normalization=combinedflat,response=normflatname)
	

# This function sets the parameters for the pyraf task immatch.imcombine and then runs it.
def run_imcombine(imagestocombine,combimgname,commongain,customRun=False):
	# This clears the current parameter values for imcombine, and sets general parameters again.
	immatch.imcombine.unlearn()
	immatch.imcombine.project='No'
	immatch.imcombine.combine='average'
	immatch.imcombine.reject='avsigclip'
	immatch.imcombine.mclip='yes'
	immatch.imcombine.nkeep=1
	immatch.imcombine.scale='mode'
	immatch.imcombine.grow=1.0
	immatch.imcombine.rdnoise=0.0
	immatch.imcombine.gain=commongain
	immatch.imcombine.snoise=0.0
	immatch.imcombine.lsigma=2.0
	immatch.imcombine.hsigma=2.0
	immatch.imcombine.blank=0.0
	immatch.imcombine.expname='EXPTIME' # average EXPTIME added to output image header
	immatch.imcombine.imcmb='$I,AIRMASS' # copy airmass to one of the IMCMBnnn keywords in output header
	
	# If running task individually, pop up epar window to let users change parameters if they want.
	if customRun == True:
		while True: 
			eparAnswer = raw_input("Do you want to further edit the parameters? 0 for no, 1 for yes.")
			if eparAnswer == '0' or eparAnswer == '1':
				break
			else:
				print "Invalid input: you must enter either 0 (no) or 1 (yes)."		
		if eparAnswer == '1':
			iraf.eparam(immatch.imcombine)
		
	immatch.imcombine(input=imagestocombine,output=combimgname)
	

# This function sets the parameters for the pyraf task imutil.imarith and then runs it.
def run_imarith(dividend,divisor,quotient,customRun=False):
	# This resets the parameters of imutil.imarith.
	imutil.imarith.unlearn()
	imutil.imarith.divzero=0.0
	# appends [1] to end of dividend image since that's where data lives
	dividend = dividend+'[1]'
	# If running task individually, pop up epar window to let users change parameters if they want.
	if customRun == True:
		while True: 
			eparAnswer = raw_input("Do you want to further edit the parameters? 0 for no, 1 for yes.")
			if eparAnswer == '0' or eparAnswer == '1':
				break
			else:
				print "Invalid input: you must enter either 0 (no) or 1 (yes)."		
		if eparAnswer == '1':
			iraf.eparam(imutil.imarith)
	# This runs imutil.imarith
	imutil.imarith(operand1=dividend,op='/',operand2=divisor,result=quotient)	
	
	
# This function sets the parameters for the pysalt task saltspec.specidentify and then runs it.
def run_specidentify(arcimage,lamplines,idfile,customRun=False):
	# this resets the parameters of saltspec.specidentify
	saltspec.specidentify.unlearn()
	saltspec.specidentify.guesstype='rss'
	saltspec.specidentify.automethod='Matchlines'
	saltspec.specidentify.function='legendre'
	saltspec.specidentify.order=3
	saltspec.specidentify.rstep=100
	saltspec.specidentify.rstart='middlerow'
	saltspec.specidentify.mdiff=5
	saltspec.specidentify.inter='yes' # interactive
	saltspec.specidentify.startext=1 # important because SALT data is in ext 1 (name: SCI)
	saltspec.specidentify.clobber='yes'
	saltspec.specidentify.logfile='pysalt.log'
	saltspec.specidentify.verbose='yes'
	# If running task individually, pop up epar window to let users change parameters if they want.
	if customRun == True:
		while True: 
			eparAnswer = raw_input("Do you want to further edit the parameters? 0 for no, 1 for yes.")
			if eparAnswer == '0' or eparAnswer == '1':
				break
			else:
				print "Invalid input: you must enter either 0 (no) or 1 (yes)."		
		if eparAnswer == '1':
			iraf.eparam(saltspec.specidentify)
	# this runs saltspec.specidentify
	saltspec.specidentify(images=arcimage,linelist=lamplines,outfile=idfile)


# This function sets the parameters for saltspec.specrectify and runs the task.	
# It outputs 2-D wavelength-corrected images (science, arc, and standard star).
def run_specrectify(input,output,idfile,customRun=False):
	# this resets the parameters of saltspec.specrectify
	saltspec.specrectify.unlearn()
	saltspec.specrectify.outpref=''
	saltspec.specrectify.caltype='line'
	saltspec.specrectify.function='legendre'
	saltspec.specrectify.order=3
	saltspec.specrectify.inttype='interp'
	saltspec.specrectify.clobber='yes'
	saltspec.specrectify.verbose='yes'
	# If running task individually, pop up epar window to let users change parameters if they want.
	if customRun == True:
		while True: 
			eparAnswer = raw_input("Do you want to further edit the parameters? 0 for no, 1 for yes.")
			if eparAnswer == '0' or eparAnswer == '1':
				break
			else:
				print "Invalid input: you must enter either 0 (no) or 1 (yes)."		
		if eparAnswer == '1':
			iraf.eparam(saltspec.specrectify)
	# this runs saltspec.specrectify
	saltspec.specrectify(images=input,outimages=output,solfile=idfile)

	
# This function sets the parameters for the pyraf task twodspec.longslit.background and then runs it.
def run_background(twodimage,newimage,faint='0',customRun=False):
	# this resets the parameters of longslit.background
	longslit.background.unlearn()
	longslit.background.axis='2'
	longslit.background.interactive='no'
	longslit.background.naverage='-100'
	longslit.background.function='legendre'
	longslit.background.order=2 # order 3 is not necessary, don't use it in interactive mode either
	longslit.background.low_reject=1.0
	longslit.background.high_reject=1.0
	longslit.background.niterate=10
	longslit.background.grow=0.0
	# since there are 2 extensions (0 and 1), need to specify data operation extension: 1
	twodimage = twodimage+'[1]'
	# runs background interactively if user indicated that spectrum is faint
	if faint=='1':
		longslit.background.interactive='yes'
	# If running task individually, pop up epar window to let users change parameters if they want.
	if customRun == True:
		while True: 
			eparAnswer = raw_input("Do you want to further edit the parameters? 0 for no, 1 for yes.")
			if eparAnswer == '0' or eparAnswer == '1':
				break
			else:
				print "Invalid input: you must enter either 0 (no) or 1 (yes)."		
		if eparAnswer == '1':		
			iraf.eparam(longslit.background)
	# This runs longslit.background
	longslit.background(input=twodimage,output=newimage)
	
# This function sets the parameters for the pyraf task apextract.apall and then runs it.
def run_apall(twodimage,spectrum,saltgain,faint,customRun=False):
	# this resets the parameters of apextract.apall.
	# with as many appropriate parameters set as possible:
	apextract.apall.unlearn()
	apextract.apall.format='multispec'
	apextract.apall.review='no'
	apextract.apall.nsum=-100 # set nsum equal to user-inputted colsum
	apextract.apall.line='INDEF' # set line equal to user-inputted col
	apextract.apall.lower=-3.5
	apextract.apall.upper=3.5
	apextract.apall.find='yes'
	apextract.apall.recenter='yes'
	apextract.apall.resize='yes'
	apextract.apall.edit='yes'
	apextract.apall.trace='yes'
	apextract.apall.fittrace='yes'
	apextract.apall.extract='yes'
	apextract.apall.extras='yes'
	apextract.apall.nfind=1 # setting this explicitly in call below so no query (for standard star)
	apextract.apall.b_function='legendre' # b. stuff = background, only used for standards
	apextract.apall.b_order=1
	apextract.apall.b_sample = '-15:-10,10:15' # clean up leftover skylines (longslit.background) around the spectrum
	apextract.apall.b_naverage = -5 # how many rows in each column to sum (positive) or median (negative)
	apextract.apall.b_low_reject=3.0
	apextract.apall.b_high_reject=3.0
	apextract.apall.b_niterate=0
	apextract.apall.b_grow=0
	apextract.apall.clean='yes' # because of leftover skylines from longslit.background
	apextract.apall.weights='variance' # because of leftover skylines from longslit.background
	apextract.apall.t_nsum=25 # increased from 15 to 25 (summing/tracing over more lines)
	apextract.apall.t_nlost=200
	apextract.apall.t_low_reject=1.0 # since there are so many data points (over 1000 generally)
	apextract.apall.t_high_reject=1.0 
	apextract.apall.t_step=20 # increased on 2013-10-14
	apextract.apall.t_niterate=2 # newly added
	apextract.apall.t_grow=1 # newly added; changed to 1 on 2013-10-14
	apextract.apall.t_function='legendre'
	apextract.apall.t_order=3
	apextract.apall.lsigma=2.0
	apextract.apall.usigma=2.0
	apextract.apall.background='fit' # because of leftover skylines from longslit.background
	apextract.apall.pfit = 'fit1d'
	apextract.apall.readnoise=3 # verify whether rdnoise is nearly the same for all images (SALT CCD)
	apextract.apall.gain=saltgain # taken from image header
	# interactive apall: yes for science; no for standards
	apextract.apall.interactive='yes'
	# changes some parameters again if user indicated that the line looks really faint
	if faint=='1': # bright without extended galaxy emission (no apall.background needed)
		apextract.apall.background='none'
	# for faint == '2' (bright with extended galaxy emission=>local bkg subtraction), the default parameters defined above are used.
	if faint=='3': # faint without extended galaxy emission (this can also be used for galaxies)
		apextract.apall.nsum=-1000
		apextract.apall.background='none'
		apextract.apall.t_nsum=50
		apextract.apall.t_step=15
		apextract.apall.t_nlost=100
		apextract.apall.t_niterate = 3 
		apextract.apall.t_naverage = 1
		apextract.apall.t_grow = 1.0
		# ask user to input a column where the SN is bright enough to start tracing from
		while True:
			colAnswer = raw_input("Please enter a column number where the SN is easily visible and can be traced: ")
			if int(colAnswer) > 0 or int(colAnswer) < 3100: # temporarily hard-coded upper limit for 2x4 binning and PG0900 grating
				break
			else:
				print "Invalid input: you must enter a column number between 0 and 3100."		
		apextract.apall.line=int(colAnswer)
	if faint=='4': # faint with extended galaxy emission (add local bkg counts; subtract local bkg)
		apextract.apall.nsum=-1000
		apextract.apall.background='fit' # to subtract local background for extended galaxy emission behind object
		apextract.apall.t_nsum=50
		apextract.apall.t_step=15
		apextract.apall.t_nlost=100
		apextract.apall.t_niterate = 3 
		apextract.apall.t_naverage = 1
		apextract.apall.t_grow = 1.0
		# ask user to input a column where the SN is bright enough to start tracing from
		while True:
			colAnswer = raw_input("Please enter a column number where the SN is easily visible and can be traced: ")
			if int(colAnswer) > 0 or int(colAnswer) < 3100: # temporarily hard-coded upper limit for 2x4 binning and PG0900 grating
				break
			else:
				print "Invalid input: you must enter a column number between 0 and 3100."		
		apextract.apall.line=int(colAnswer)
	first = twodimage[:3]
	if first == 'sci':
		apextract.apall.interactive='yes'
		print 'running apall interactively on science image: '+twodimage
	elif first=='std':
		apextract.apall.interactive='no'
		apextract.apall.background='fit' # since standards aren't background-subtracted; may not be necessary though
		print 'running apall non-interactively on standard star image: '+twodimage
	# since there are 2 extensions (0 and 1), need to specify extraction extension: 1
	twodimage = twodimage+'[1]'
	# If running task individually, pop up epar window to let users change parameters if they want.
	if customRun == True:
		while True: 
			eparAnswer = raw_input("Do you want to further edit the parameters? 0 for no, 1 for yes.")
			if eparAnswer == '0' or eparAnswer == '1':
				break
			else:
				print "Invalid input: you must enter either 0 (no) or 1 (yes)."		
		if eparAnswer == '1':		
			iraf.eparam(apextract.apall)
	# This runs apextract.apall
	apextract.apall(input=twodimage,output=spectrum,nfind=1,trace='yes',fittrace='yes',recenter='yes',resize='yes',edit='yes',extract='yes')


# This function sets the parameters for the pyraf task apextract.apsum and then runs it.  
def run_apsum(twodimage,refimage,spectrum,customRun=False):
	# this resets the parameters of apextract.apsum
	apextract.apsum.unlearn()
	apextract.apsum.interactive='no' # non-interactive
	apextract.apsum.review='no'
	apextract.apsum.background='none'
	apextract.apsum.format='multispec'
	# apextract.apsum.clean='yes'
	# apextract.apsum.weights='variance'
	apextract.apsum.nsum=50
	apextract.apsum.lsigma=2.0
	apextract.apsum.usigma=2.0
	# need to specify extraction extension (1, not 0) since there are 2
	twodimage = twodimage+'[1]'
	# If running task individually, pop up epar window to let users change parameters if they want.
	if customRun == True:
		while True: 
			eparAnswer = raw_input("Do you want to further edit the parameters? 0 for no, 1 for yes.")
			if eparAnswer == '0' or eparAnswer == '1':
				break
			else:
				print "Invalid input: you must enter either 0 (no) or 1 (yes)."		
		if eparAnswer == '1':
			iraf.eparam(apextract.apsum)
	# This runs apextract.apsum
	apextract.apsum(input=twodimage,output=spectrum,references=refimage)


# This function sets the parameters for the pyraf task apextract.apsum and then runs it
def run_apsumSci(inputimage,refimage,spectrum,faint,customRun=False):
	# this resets the parameters of apextract.apsum
	apextract.apsum.unlearn()
	apextract.apsum.interactive='no' # non-interactive
	apextract.apsum.review='no'
	apextract.apsum.background='fit' # we want the background subtracted for science
	apextract.apsum.extras = 'yes' # we want the additional spectral bands (sky (3), sigma (4) especially)
	apextract.apsum.format='multispec'
	apextract.apsum.clean='yes' # we want cosmic rays and bad pixels cleaned for science
	apextract.apsum.weights='variance' # we want a weighted extraction for science
	apextract.apsum.pfit='fit1d'
	apextract.apsum.nsum=50
	apextract.apsum.lsigma=2.0
	apextract.apsum.usigma=2.0
	# this makes additional parameter changes if spectrum is faint
	if faint == '1':
		apextract.apsum.nsum=-1000 # taken from apall for faint spectra above
	# this sets the gain parameter as specified in inputimage's 0-header
	hduin = pyfits.open(inputimage[:len(inputimage)-3]) # pyfits open fits file, not fits[1] "file"
	hdr = hduin[0].header
	thisgain = hdr.get('GAIN',2.146) # default for our usual SALT long-slit spectroscopy settings
	hduin.close()
	apextract.apsum.gain=thisgain
	# If running task individually, pop up epar window to let users change parameters if they want.
	if customRun == True:
		while True: 
			eparAnswer = raw_input("Do you want to further edit the parameters? 0 for no, 1 for yes.")
			if eparAnswer == '0' or eparAnswer == '1':
				break
			else:
				print "Invalid input: you must enter either 0 (no) or 1 (yes)."		
		if eparAnswer == '1':		
			iraf.eparam(apextract.apsum)
	# This runs apextract.apsum
	apextract.apsum(input=inputimage,output=spectrum,references=refimage)
	
	
# This function sets the parameters for and runs identify
def run_identify(arcsci,lamplines,customRun=False):
	# this resets the parameters of onedspec.identify
	onedspec.identify.unlearn()
	onedspec.identify.nsum=10
	onedspec.identify.match=-3.0
	onedspec.identify.maxfeatures=50
	onedspec.identify.zwidth=100.0
	onedspec.identify.ftype='emission'
	onedspec.identify.threshold=0.0
	onedspec.identify.function='spline3'
	onedspec.identify.order=1
	onedspec.identify.niterate=0
	onedspec.identify.low_reject=3.0
	onedspec.identify.high_reject=3.0
	onedspec.identify.grow=0.0
	onedspec.identify.autowrite='No'
	onedspec.identify.database='database'
	onedspec.identify.section='middle line'
	# If running task individually, pop up epar window to let users change parameters if they want.
	if customRun == True:
		while True: 
			eparAnswer = raw_input("Do you want to further edit the parameters? 0 for no, 1 for yes.")
			if eparAnswer == '0' or eparAnswer == '1':
				break
			else:
				print "Invalid input: you must enter either 0 (no) or 1 (yes)."		
		if eparAnswer == '1':			
			iraf.eparam(onedspec.identify)
	# this runs longslit.identify
	onedspec.identify(images=arcsci,coordlist=lamplines)
	

# This function sets the parameters for and runs the pyraf task reidentify.
def run_reidentify(input,ref,lamplines,customRun=False):
	# this resets the parameters of onedspec.identify
	onedspec.reidentify.unlearn()
	onedspec.reidentify.interactive='no'
	onedspec.reidentify.newaps='yes'
	onedspec.reidentify.override='no'
	onedspec.reidentify.refit='yes'
	onedspec.reidentify.trace='yes'
	onedspec.reidentify.step=10
	onedspec.reidentify.nsum=10
	onedspec.reidentify.shift=0.
	onedspec.reidentify.search=0.0
	onedspec.reidentify.nlost=0
	onedspec.reidentify.threshold=0.0
	onedspec.reidentify.addfeatures='no'
	onedspec.reidentify.match=-3.0
	onedspec.reidentify.maxfeatures=50
	onedspec.reidentify.verbose='yes'
	# If running task individually, pop up epar window to let users change parameters if they want.
	if customRun == True:
		while True: 
			eparAnswer = raw_input("Do you want to further edit the parameters? 0 for no, 1 for yes.")
			if eparAnswer == '0' or eparAnswer == '1':
				break
			else:
				print "Invalid input: you must enter either 0 (no) or 1 (yes)."		
		if eparAnswer == '1':		
			iraf.eparam(onedspec.reidentify)
	# this runs longslit.reidentify non-interactively
	onedspec.reidentify(reference=ref,images=input,coordlist=lamplines)


# This function sets the parameters for the pyraf task onedspec.dispcor and then runs it.
def run_dispcor(original,corrected,customRun=False):
	# this resets the parameters of onedspec.dispcor
	onedspec.dispcor.unlearn()
	onedspec.dispcor.verbose='yes'
	# If running task individually, pop up epar window to let users change parameters if they want.
	if customRun == True:
		while True: 
			eparAnswer = raw_input("Do you want to further edit the parameters? 0 for no, 1 for yes.")
			if eparAnswer == '0' or eparAnswer == '1':
				break
			else:
				print "Invalid input: you must enter either 0 (no) or 1 (yes)."		
		if eparAnswer == '1':		
			iraf.eparam(onedspec.dispcor)
	# This runs onedspec.dispcor (REFSPEC1 keyword addition was done by identifyarcs and reidentifyarcs)
	onedspec.dispcor(input=original,output=corrected)
	

# This function sets the parameters for the pyraf task onedspec.standard and then runs it.
# The user has already been alerted by fluxcal function above about the star_name.
def run_standard(standardimage,outputname,exp,air,starName,customRun=False):
	# this resets the parameters of onedspec.standard
	onedspec.standard.unlearn()
	onedspec.standard.caldir=params.standardsPath
	onedspec.standard.extinction=''
	# onedspec.standard.apertures='1' # only want to interactively choose bandpasses in aperture 1
	onedspec.standard.interact='NO' # non-interactive defining of bandpasses, need at least 15 bandpasses 
	onedspec.standard.airmass=air
	onedspec.standard.exptime=exp
	onedspec.standard.answer='NO'
	# magnitude of standard star (apparent/absolute) not given to us in SALT headers
	onedspec.standard.bandwidth = 50.0 # automatic definition of bandpasses for non-interactive flux calibration
	onedspec.standard.bandsep = 20.0
	# the star_name is the same as the root name of the dat file in the salt standards caldir
	onedspec.standard.star_name = starName
	# If running task individually, pop up epar window to let users change parameters if they want.
	if customRun == True:
		onedspec.standard.interact='YES' # give option to run interactively
		while True: 
			eparAnswer = raw_input("Do you want to further edit the parameters? 0 for no, 1 for yes.")
			if eparAnswer == '0' or eparAnswer == '1':
				break
			else:
				print "Invalid input: you must enter either 0 (no) or 1 (yes)."		
		if eparAnswer == '1':		
			iraf.eparam(onedspec.standard)
	# this runs onedspec.standard resulting in an std file named outputname 
	onedspec.standard(input=standardimage+'[1]',output=outputname)
	

# This function sets the parameters for the pyraf task onedspec.sensfunc and then runs it.
def run_sensfunc(stddata,sensname,customRun=False):
	# this resets the parameters of onedspec.sensfunc
	onedspec.sensfunc.unlearn()
	onedspec.sensfunc.apertures='1' # use only the first aperture (object data)
	onedspec.sensfunc.function='spline3'
	onedspec.sensfunc.order=2
	onedspec.sensfunc.extinction=''
	onedspec.sensfunc.ignoreaps='yes' # create only one sensfunc file
	onedspec.sensfunc.interactive='NO' # non-interactive fitting 
	onedspec.sensfunc.graphs='sri'
	onedspec.sensfunc.answer='NO' # user shouldn't have to input anything else, just check fit and quit
	# If running task individually, pop up epar window to let users change parameters if they want.
	if customRun == True:
		onedspec.sensfunc.interactive='YES' # give option to run interactively
		while True: 
			eparAnswer = raw_input("Do you want to further edit the parameters? 0 for no, 1 for yes.")
			if eparAnswer == '0' or eparAnswer == '1':
				break
			else:
				print "Invalid input: you must enter either 0 (no) or 1 (yes)."		
		if eparAnswer == '1':		
			iraf.eparam(onedspec.sensfunc)
	# runs onedspec.sensfunc resulting in a file named sensname
	onedspec.sensfunc(standards=stddata,sensitivity=sensname)
	

# This function sets the parameters for the pyraf task onedspec.calibrate and then runs it.
def run_calibrate(scienceimage,fluximage,sensfilename,customRun=False):
	# this resets the parameters of onedspec.calibrate
	onedspec.calibrate.unlearn()
	onedspec.calibrate.extinct='no'
	onedspec.calibrate.extinction=''
	onedspec.calibrate.sensitivity=sensfilename # the sensitivity function for flux calibration
	onedspec.calibrate.ignoreaps='yes' # look for and use sens*, not sens*0001, sens*0002, etc.
	# This stores some header information from the science image
	scihdr = pyfits.getheader(scienceimage,1)
	exptime = scihdr['EXPTIME']
	airmass = abs(scihdr['AIRMASS'])
	# this continues to set the necessary parameters of onedspec.standard
	onedspec.calibrate.airmass=airmass
	onedspec.calibrate.exptime=exptime
	
	# runs onedspec.calibrate resulting in a flux-calibrated spectrum with name fluximage (on [SCI] extension)
	scienceimage1 = scienceimage+'[1]'
	
	# If running task individually, pop up epar window to let users change parameters if they want.
	if customRun == True:
		while True: 
			eparAnswer = raw_input("Do you want to further edit the parameters? 0 for no, 1 for yes.")
			if eparAnswer == '0' or eparAnswer == '1':
				break
			else:
				print "Invalid input: you must enter either 0 (no) or 1 (yes)."		
		if eparAnswer == '1':		
			iraf.eparam(onedspec.calibrate)
	
	onedspec.calibrate(input=scienceimage1,output=fluximage)
	

# This function sets the parameters for the pyraf task onedspec.odcombine and then runs it.
def run_odcombine(inputlistseq,finalspectrumname,saltgain,customRun=False):
	# this resets the parameters of odcombine
	onedspec.odcombine.unlearn()
	onedspec.odcombine.group='all'
	onedspec.odcombine.combine='average'
	onedspec.odcombine.reject='avsigclip'
	onedspec.odcombine.apertures = 1
	onedspec.odcombine.outtype = 'real'
	onedspec.odcombine.smaskformat = 'bpmspectrum' # each input image has 'BPM' keyword linking it to its bpm
	onedspec.odcombine.smasktype = 'goodvalue'
	onedspec.odcombine.smaskvalue = 0.0
	onedspec.odcombine.blank = 0.0
	onedspec.odcombine.lsigma = 2.0
	onedspec.odcombine.hsigma = 2.0
	onedspec.odcombine.gain = saltgain
	onedspec.odcombine.masktype = 'goodvalue'
	onedspec.odcombine.maskvalue = 0.0
	# If running task individually, pop up epar window to let users change parameters if they want.
	if customRun == True:
		while True: 
			eparAnswer = raw_input("Do you want to further edit the parameters? 0 for no, 1 for yes.")
			if eparAnswer == '0' or eparAnswer == '1':
				break
			else:
				print "Invalid input: you must enter either 0 (no) or 1 (yes)."		
		if eparAnswer == '1':
			iraf.eparam(onedspec.odcombine)
	# running onedspec.odcombine resulting in the final wave- and flux-calibrated spectrum named finalspectrumname
	onedspec.odcombine(input=inputlistseq,output=finalspectrumname)
	

# This function creates a copy of a FITS file and preserves the 2-extension structure.
def run_imcopy(inputname,outputname,numExt=2,customRun=False):	
	if numExt == 1:
		imutil.imcopy(input=inputname+'[0]',output=outputname+'[append]')
	else:
		imutil.imcopy(input=inputname+'[0]',output=outputname+'[append]')
		imutil.imcopy(input=inputname+'[SCI]',output=outputname+'[append]')
		
		
# This function deletes a file.
def run_imdel(inputname,customRun=False):
	imutil.imdel(inputname)
	

# This function extracts bands from a spectrum stored in a FITS file.
def run_scopy(inputname,outputname,band,customRun=False):
	onedspec.scopy.unlearn()
	onedspec.scopy(input=inputname+'[1]',output=outputname,bands=band,format='multispec',clobber='yes',verbose='yes')

# This *generalized* function sets the parameters for the pyraf task imutil.imarith and then runs it.
def run_imarithGeneral(op1,op2,oper,outname,customRun=False):

	# This resets the parameters of imutil.imarith.
	imutil.imarith.unlearn()
	imutil.imarith.divzero=0.0
	# appends [1] to end of dividend image since that's where data lives
	# If running task individually, pop up epar window to let users change parameters if they want.
	if customRun == True:
		while True: 
			eparAnswer = raw_input("Do you want to further edit the parameters? 0 for no, 1 for yes.")
			if eparAnswer == '0' or eparAnswer == '1':
				break
			else:
				print "Invalid input: you must enter either 0 (no) or 1 (yes)."		
		if eparAnswer == '1':
			iraf.eparam(imutil.imarith)
	# This runs imutil.imarith
	imutil.imarith(operand1=op1,op=oper,operand2=op2,result=outname)	
	
	
	
	
			
=======
'''
>>> Rutgers SALT Supernova Spectral Reduction Pipeline <<<

This module contains functions which run a specific PyRAF task
corresponding to one step of the entire reduction process.

Each function sets the parameters for a PyRAF task and then
runs it on the data. These functions are primarily called by 
the core module containing the main pipeline functions: tasks.py.

Please refer to the documentation for more information about
how each PyRAF call affects the data.

*** Modifications ***
Sept. 26, 2013: Created module. -Viraj Pandya

'''

import pyfits

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

import dicts # These are the pipeline's global dictionaries.
import params # Customizable parameters for the pipeline.

# This function sets the parameters for the pyraf task ccdred.flatcombine and then runs it.
def run_flatcombine(flatstocombine,combflatname,customRun=False):
	# This clears the current parameter values for flatcombine, and sets general parameters again.
	ccdred.flatcombine.unlearn()
	ccdred.flatcombine.combine='average'
	ccdred.flatcombine.reject='avsigclip'
	ccdred.flatcombine.ccdtype=''
	ccdred.flatcombine.scale='mode'
	ccdred.flatcombine.process='no'
	ccdred.flatcombine.subsets='no'
	
	# appends '[1]' to end of each filename since that's the extension for the data to be combined
	flatstofeed = list(flatstocombine) # making copy of list (will need unmodified original for imcopy)
	for i,v in enumerate(flatstofeed):
		flatstofeed[i] = v+'[1]'
	
	# turn python list (of string filenames) into a string sequence
	flatcombineseq = ','.join(flatstofeed)
	# If running task individually, pop up epar window to let users change parameters if they want.
	if customRun == True:
		while True: 
			eparAnswer = raw_input("Do you want to further edit the parameters? 0 for no, 1 for yes.")
			if eparAnswer == '0' or eparAnswer == '1':
				break
			else:
				print "Invalid input: you must enter either 0 (no) or 1 (yes)."		
		if eparAnswer == '1':
			iraf.eparam(ccdred.flatcombine)
	
	ccdred.flatcombine(input=flatcombineseq,output=combflatname)

# This function sets the parameters for the pyraf task longslit.response and then runs it.
def run_response(combinedflat,normflatname,customRun=False):
	# this clears current parameter values for longslit.response and sets general parameters.
	longslit.response.unlearn()
	longslit.response.threshold=5.0
	longslit.response.function='spline3'
	longslit.response.order=19
	longslit.response.low_rej=3.0
	longslit.response.high_rej=3.0
	longslit.response.niterate=3
	longslit.response.grow=2
	longslit.response.interactive='no' # 2013-12-06: prevents bug in iraf that auto-inputs dispaxis=2
	# If running task individually, pop up epar window to let users change parameters if they want.
	if customRun == True:
		while True: 
			eparAnswer = raw_input("Do you want to further edit the parameters? 0 for no, 1 for yes.")
			if eparAnswer == '0' or eparAnswer == '1':
				break
			else:
				print "Invalid input: you must enter either 0 (no) or 1 (yes)."		
		if eparAnswer == '1':
			iraf.eparam(longslit.response)
	# This makes sure that the DISPAXIS is set correctly (if 'DISPAXIS' keyword not in header)
	longslit.dispaxis = 1
	# This runs longslit.response (interactively for now)
	longslit.response(calibration=combinedflat,normalization=combinedflat,response=normflatname)
	

# This function sets the parameters for the pyraf task immatch.imcombine and then runs it.
def run_imcombine(imagestocombine,combimgname,commongain,customRun=False):
	# This clears the current parameter values for imcombine, and sets general parameters again.
	immatch.imcombine.unlearn()
	immatch.imcombine.project='No'
	immatch.imcombine.combine='average'
	immatch.imcombine.reject='avsigclip'
	immatch.imcombine.mclip='yes'
	immatch.imcombine.nkeep=1
	immatch.imcombine.scale='mode'
	immatch.imcombine.grow=1.0
	immatch.imcombine.rdnoise=0.0
	immatch.imcombine.gain=commongain
	immatch.imcombine.snoise=0.0
	immatch.imcombine.lsigma=2.0
	immatch.imcombine.hsigma=2.0
	immatch.imcombine.blank=0.0
	immatch.imcombine.expname='EXPTIME' # average EXPTIME added to output image header
	immatch.imcombine.imcmb='$I,AIRMASS' # copy airmass to one of the IMCMBnnn keywords in output header
	
	# If running task individually, pop up epar window to let users change parameters if they want.
	if customRun == True:
		while True: 
			eparAnswer = raw_input("Do you want to further edit the parameters? 0 for no, 1 for yes.")
			if eparAnswer == '0' or eparAnswer == '1':
				break
			else:
				print "Invalid input: you must enter either 0 (no) or 1 (yes)."		
		if eparAnswer == '1':
			iraf.eparam(immatch.imcombine)
		
	immatch.imcombine(input=imagestocombine,output=combimgname)
	

# This function sets the parameters for the pyraf task imutil.imarith and then runs it.
def run_imarith(dividend,divisor,quotient,customRun=False):
	# This resets the parameters of imutil.imarith.
	imutil.imarith.unlearn()
	imutil.imarith.divzero=0.0
	# appends [1] to end of dividend image since that's where data lives
	dividend = dividend+'[1]'
	# If running task individually, pop up epar window to let users change parameters if they want.
	if customRun == True:
		while True: 
			eparAnswer = raw_input("Do you want to further edit the parameters? 0 for no, 1 for yes.")
			if eparAnswer == '0' or eparAnswer == '1':
				break
			else:
				print "Invalid input: you must enter either 0 (no) or 1 (yes)."		
		if eparAnswer == '1':
			iraf.eparam(imutil.imarith)
	# This runs imutil.imarith
	imutil.imarith(operand1=dividend,op='/',operand2=divisor,result=quotient)	
	
	
# This function sets the parameters for the pysalt task saltspec.specidentify and then runs it.
def run_specidentify(arcimage,lamplines,idfile,customRun=False):
	# this resets the parameters of saltspec.specidentify
	saltspec.specidentify.unlearn()
	saltspec.specidentify.guesstype='rss'
	saltspec.specidentify.automethod='Matchlines'
	saltspec.specidentify.function='legendre'
	saltspec.specidentify.order=3
	saltspec.specidentify.rstep=100
	saltspec.specidentify.rstart='middlerow'
	saltspec.specidentify.mdiff=5
	saltspec.specidentify.inter='yes' # interactive
	saltspec.specidentify.startext=1 # important because SALT data is in ext 1 (name: SCI)
	saltspec.specidentify.clobber='yes'
	saltspec.specidentify.logfile='pysalt.log'
	saltspec.specidentify.verbose='yes'
	# If running task individually, pop up epar window to let users change parameters if they want.
	if customRun == True:
		while True: 
			eparAnswer = raw_input("Do you want to further edit the parameters? 0 for no, 1 for yes.")
			if eparAnswer == '0' or eparAnswer == '1':
				break
			else:
				print "Invalid input: you must enter either 0 (no) or 1 (yes)."		
		if eparAnswer == '1':
			iraf.eparam(saltspec.specidentify)
	# this runs saltspec.specidentify
	saltspec.specidentify(images=arcimage,linelist=lamplines,outfile=idfile)


# This function sets the parameters for saltspec.specrectify and runs the task.	
# It outputs 2-D wavelength-corrected images (science, arc, and standard star).
def run_specrectify(input,output,idfile,customRun=False):
	# this resets the parameters of saltspec.specrectify
	saltspec.specrectify.unlearn()
	saltspec.specrectify.outpref=''
	saltspec.specrectify.caltype='line'
	saltspec.specrectify.function='legendre'
	saltspec.specrectify.order=3
	saltspec.specrectify.inttype='interp'
	saltspec.specrectify.clobber='yes'
	saltspec.specrectify.verbose='yes'
	# If running task individually, pop up epar window to let users change parameters if they want.
	if customRun == True:
		while True: 
			eparAnswer = raw_input("Do you want to further edit the parameters? 0 for no, 1 for yes.")
			if eparAnswer == '0' or eparAnswer == '1':
				break
			else:
				print "Invalid input: you must enter either 0 (no) or 1 (yes)."		
		if eparAnswer == '1':
			iraf.eparam(saltspec.specrectify)
	# this runs saltspec.specrectify
	saltspec.specrectify(images=input,outimages=output,solfile=idfile)

	
# This function sets the parameters for the pyraf task twodspec.longslit.background and then runs it.
def run_background(twodimage,newimage,faint='0',customRun=False):
	# this resets the parameters of longslit.background
	longslit.background.unlearn()
	longslit.background.axis='2'
	longslit.background.interactive='no'
	longslit.background.naverage='-100'
	longslit.background.function='legendre'
	longslit.background.order=2 # order 3 is not necessary, don't use it in interactive mode either
	longslit.background.low_reject=1.0
	longslit.background.high_reject=1.0
	longslit.background.niterate=10
	longslit.background.grow=0.0
	# since there are 2 extensions (0 and 1), need to specify data operation extension: 1
	twodimage = twodimage+'[1]'
	# runs background interactively if user indicated that spectrum is faint
	if faint=='1':
		longslit.background.interactive='yes'
	# If running task individually, pop up epar window to let users change parameters if they want.
	if customRun == True:
		while True: 
			eparAnswer = raw_input("Do you want to further edit the parameters? 0 for no, 1 for yes.")
			if eparAnswer == '0' or eparAnswer == '1':
				break
			else:
				print "Invalid input: you must enter either 0 (no) or 1 (yes)."		
		if eparAnswer == '1':		
			iraf.eparam(longslit.background)
	# This runs longslit.background
	longslit.background(input=twodimage,output=newimage)
	
# This function sets the parameters for the pyraf task apextract.apall and then runs it.
def run_apall(twodimage,spectrum,saltgain,faint,customRun=False):
	# this resets the parameters of apextract.apall.
	# with as many appropriate parameters set as possible:
	apextract.apall.unlearn()
	apextract.apall.format='multispec'
	apextract.apall.review='no'
	apextract.apall.nsum=-100 # set nsum equal to user-inputted colsum
	apextract.apall.line='INDEF' # set line equal to user-inputted col
	apextract.apall.lower=-3.5
	apextract.apall.upper=3.5
	apextract.apall.find='yes'
	apextract.apall.recenter='yes'
	apextract.apall.resize='yes'
	apextract.apall.edit='yes'
	apextract.apall.trace='yes'
	apextract.apall.fittrace='yes'
	apextract.apall.extract='yes'
	apextract.apall.extras='yes'
	apextract.apall.nfind=1 # setting this explicitly in call below so no query (for standard star)
	apextract.apall.b_function='legendre' # b. stuff = background, only used for standards
	apextract.apall.b_order=1
	apextract.apall.b_sample = '-15:-10,10:15' # clean up leftover skylines (longslit.background) around the spectrum
	apextract.apall.b_naverage = -5 # how many rows in each column to sum (positive) or median (negative)
	apextract.apall.b_low_reject=3.0
	apextract.apall.b_high_reject=3.0
	apextract.apall.b_niterate=0
	apextract.apall.b_grow=0
	apextract.apall.clean='yes' # because of leftover skylines from longslit.background
	apextract.apall.weights='variance' # because of leftover skylines from longslit.background
	apextract.apall.t_nsum=25 # increased from 15 to 25 (summing/tracing over more lines)
	apextract.apall.t_nlost=200
	apextract.apall.t_low_reject=1.0 # since there are so many data points (over 1000 generally)
	apextract.apall.t_high_reject=1.0 
	apextract.apall.t_step=20 # increased on 2013-10-14
	apextract.apall.t_niterate=2 # newly added
	apextract.apall.t_grow=1 # newly added; changed to 1 on 2013-10-14
	apextract.apall.t_function='legendre'
	apextract.apall.t_order=3
	apextract.apall.lsigma=2.0
	apextract.apall.usigma=2.0
	apextract.apall.background='fit' # because of leftover skylines from longslit.background
	apextract.apall.pfit = 'fit1d'
	apextract.apall.readnoise=3 # verify whether rdnoise is nearly the same for all images (SALT CCD)
	apextract.apall.gain=saltgain # taken from image header
	# interactive apall: yes for science; no for standards
	apextract.apall.interactive='yes'
	# changes some parameters again if user indicated that the line looks really faint
	if faint=='1': # bright without extended galaxy emission (no apall.background needed)
		apextract.apall.background='none'
	# for faint == '2' (bright with extended galaxy emission=>local bkg subtraction), the default parameters defined above are used.
	if faint=='3': # faint without extended galaxy emission (this can also be used for galaxies)
		apextract.apall.nsum=-1000
		apextract.apall.background='none'
		apextract.apall.t_nsum=50
		apextract.apall.t_step=15
		apextract.apall.t_nlost=100
		apextract.apall.t_niterate = 3 
		apextract.apall.t_naverage = 1
		apextract.apall.t_grow = 1.0
		# ask user to input a column where the SN is bright enough to start tracing from
		while True:
			colAnswer = raw_input("Please enter a column number where the SN is easily visible and can be traced: ")
			if int(colAnswer) > 0 or int(colAnswer) < 3100: # temporarily hard-coded upper limit for 2x4 binning and PG0900 grating
				break
			else:
				print "Invalid input: you must enter a column number between 0 and 3100."		
		apextract.apall.line=int(colAnswer)
	if faint=='4': # faint with extended galaxy emission (add local bkg counts; subtract local bkg)
		apextract.apall.nsum=-1000
		apextract.apall.background='fit' # to subtract local background for extended galaxy emission behind object
		apextract.apall.t_nsum=50
		apextract.apall.t_step=15
		apextract.apall.t_nlost=100
		apextract.apall.t_niterate = 3 
		apextract.apall.t_naverage = 1
		apextract.apall.t_grow = 1.0
		# ask user to input a column where the SN is bright enough to start tracing from
		while True:
			colAnswer = raw_input("Please enter a column number where the SN is easily visible and can be traced: ")
			if int(colAnswer) > 0 or int(colAnswer) < 3100: # temporarily hard-coded upper limit for 2x4 binning and PG0900 grating
				break
			else:
				print "Invalid input: you must enter a column number between 0 and 3100."		
		apextract.apall.line=int(colAnswer)
	first = twodimage[:3]
	if first == 'sci':
		apextract.apall.interactive='yes'
		print 'running apall interactively on science image: '+twodimage
	elif first=='std':
		apextract.apall.interactive='no'
		apextract.apall.background='fit' # since standards aren't background-subtracted; may not be necessary though
		print 'running apall non-interactively on standard star image: '+twodimage
	# since there are 2 extensions (0 and 1), need to specify extraction extension: 1
	twodimage = twodimage+'[1]'
	# If running task individually, pop up epar window to let users change parameters if they want.
	if customRun == True:
		while True: 
			eparAnswer = raw_input("Do you want to further edit the parameters? 0 for no, 1 for yes.")
			if eparAnswer == '0' or eparAnswer == '1':
				break
			else:
				print "Invalid input: you must enter either 0 (no) or 1 (yes)."		
		if eparAnswer == '1':		
			iraf.eparam(apextract.apall)
	# This runs apextract.apall
	apextract.apall(input=twodimage,output=spectrum,nfind=1,trace='yes',fittrace='yes',recenter='yes',resize='yes',edit='yes',extract='yes')


# This function sets the parameters for the pyraf task apextract.apsum and then runs it.  
def run_apsum(twodimage,refimage,spectrum,customRun=False):
	# this resets the parameters of apextract.apsum
	apextract.apsum.unlearn()
	apextract.apsum.interactive='no' # non-interactive
	apextract.apsum.review='no'
	apextract.apsum.background='none'
	apextract.apsum.format='multispec'
	# apextract.apsum.clean='yes'
	# apextract.apsum.weights='variance'
	apextract.apsum.nsum=50
	apextract.apsum.lsigma=2.0
	apextract.apsum.usigma=2.0
	# need to specify extraction extension (1, not 0) since there are 2
	twodimage = twodimage+'[1]'
	# If running task individually, pop up epar window to let users change parameters if they want.
	if customRun == True:
		while True: 
			eparAnswer = raw_input("Do you want to further edit the parameters? 0 for no, 1 for yes.")
			if eparAnswer == '0' or eparAnswer == '1':
				break
			else:
				print "Invalid input: you must enter either 0 (no) or 1 (yes)."		
		if eparAnswer == '1':
			iraf.eparam(apextract.apsum)
	# This runs apextract.apsum
	apextract.apsum(input=twodimage,output=spectrum,references=refimage)


# This function sets the parameters for the pyraf task apextract.apsum and then runs it
def run_apsumSci(inputimage,refimage,spectrum,faint,customRun=False):
	# this resets the parameters of apextract.apsum
	apextract.apsum.unlearn()
	apextract.apsum.interactive='no' # non-interactive
	apextract.apsum.review='no'
	apextract.apsum.background='fit' # we want the background subtracted for science
	apextract.apsum.extras = 'yes' # we want the additional spectral bands (sky (3), sigma (4) especially)
	apextract.apsum.format='multispec'
	apextract.apsum.clean='yes' # we want cosmic rays and bad pixels cleaned for science
	apextract.apsum.weights='variance' # we want a weighted extraction for science
	apextract.apsum.pfit='fit1d'
	apextract.apsum.nsum=50
	apextract.apsum.lsigma=2.0
	apextract.apsum.usigma=2.0
	# this makes additional parameter changes if spectrum is faint
	if faint == '1':
		apextract.apsum.nsum=-1000 # taken from apall for faint spectra above
	# this sets the gain parameter as specified in inputimage's 0-header
	hduin = pyfits.open(inputimage[:len(inputimage)-3]) # pyfits open fits file, not fits[1] "file"
	hdr = hduin[0].header
	thisgain = hdr.get('GAIN',2.146) # default for our usual SALT long-slit spectroscopy settings
	hduin.close()
	apextract.apsum.gain=thisgain
	# If running task individually, pop up epar window to let users change parameters if they want.
	if customRun == True:
		while True: 
			eparAnswer = raw_input("Do you want to further edit the parameters? 0 for no, 1 for yes.")
			if eparAnswer == '0' or eparAnswer == '1':
				break
			else:
				print "Invalid input: you must enter either 0 (no) or 1 (yes)."		
		if eparAnswer == '1':		
			iraf.eparam(apextract.apsum)
	# This runs apextract.apsum
	apextract.apsum(input=inputimage,output=spectrum,references=refimage)
	
	
# This function sets the parameters for and runs identify
def run_identify(arcsci,lamplines,customRun=False):
	# this resets the parameters of onedspec.identify
	onedspec.identify.unlearn()
	onedspec.identify.nsum=10
	onedspec.identify.match=-3.0
	onedspec.identify.maxfeatures=50
	onedspec.identify.zwidth=100.0
	onedspec.identify.ftype='emission'
	onedspec.identify.threshold=0.0
	onedspec.identify.function='spline3'
	onedspec.identify.order=1
	onedspec.identify.niterate=0
	onedspec.identify.low_reject=3.0
	onedspec.identify.high_reject=3.0
	onedspec.identify.grow=0.0
	onedspec.identify.autowrite='No'
	onedspec.identify.database='database'
	onedspec.identify.section='middle line'
	# If running task individually, pop up epar window to let users change parameters if they want.
	if customRun == True:
		while True: 
			eparAnswer = raw_input("Do you want to further edit the parameters? 0 for no, 1 for yes.")
			if eparAnswer == '0' or eparAnswer == '1':
				break
			else:
				print "Invalid input: you must enter either 0 (no) or 1 (yes)."		
		if eparAnswer == '1':			
			iraf.eparam(onedspec.identify)
	# this runs longslit.identify
	onedspec.identify(images=arcsci,coordlist=lamplines)
	

# This function sets the parameters for and runs the pyraf task reidentify.
def run_reidentify(input,ref,lamplines,customRun=False):
	# this resets the parameters of onedspec.identify
	onedspec.reidentify.unlearn()
	onedspec.reidentify.interactive='no'
	onedspec.reidentify.newaps='yes'
	onedspec.reidentify.override='no'
	onedspec.reidentify.refit='yes'
	onedspec.reidentify.trace='yes'
	onedspec.reidentify.step=10
	onedspec.reidentify.nsum=10
	onedspec.reidentify.shift=0.
	onedspec.reidentify.search=0.0
	onedspec.reidentify.nlost=0
	onedspec.reidentify.threshold=0.0
	onedspec.reidentify.addfeatures='no'
	onedspec.reidentify.match=-3.0
	onedspec.reidentify.maxfeatures=50
	onedspec.reidentify.verbose='yes'
	# If running task individually, pop up epar window to let users change parameters if they want.
	if customRun == True:
		while True: 
			eparAnswer = raw_input("Do you want to further edit the parameters? 0 for no, 1 for yes.")
			if eparAnswer == '0' or eparAnswer == '1':
				break
			else:
				print "Invalid input: you must enter either 0 (no) or 1 (yes)."		
		if eparAnswer == '1':		
			iraf.eparam(onedspec.reidentify)
	# this runs longslit.reidentify non-interactively
	onedspec.reidentify(reference=ref,images=input,coordlist=lamplines)


# This function sets the parameters for the pyraf task onedspec.dispcor and then runs it.
def run_dispcor(original,corrected,customRun=False):
	# this resets the parameters of onedspec.dispcor
	onedspec.dispcor.unlearn()
	onedspec.dispcor.verbose='yes'
	# If running task individually, pop up epar window to let users change parameters if they want.
	if customRun == True:
		while True: 
			eparAnswer = raw_input("Do you want to further edit the parameters? 0 for no, 1 for yes.")
			if eparAnswer == '0' or eparAnswer == '1':
				break
			else:
				print "Invalid input: you must enter either 0 (no) or 1 (yes)."		
		if eparAnswer == '1':		
			iraf.eparam(onedspec.dispcor)
	# This runs onedspec.dispcor (REFSPEC1 keyword addition was done by identifyarcs and reidentifyarcs)
	onedspec.dispcor(input=original,output=corrected)
	

# This function sets the parameters for the pyraf task onedspec.standard and then runs it.
# The user has already been alerted by fluxcal function above about the star_name.
def run_standard(standardimage,outputname,exp,air,starName,customRun=False):
	# this resets the parameters of onedspec.standard
	onedspec.standard.unlearn()
	onedspec.standard.caldir=params.standardsPath
	onedspec.standard.extinction=''
	# onedspec.standard.apertures='1' # only want to interactively choose bandpasses in aperture 1
	onedspec.standard.interact='NO' # non-interactive defining of bandpasses, need at least 15 bandpasses 
	onedspec.standard.airmass=air
	onedspec.standard.exptime=exp
	onedspec.standard.answer='NO'
	# magnitude of standard star (apparent/absolute) not given to us in SALT headers
	onedspec.standard.bandwidth = 50.0 # automatic definition of bandpasses for non-interactive flux calibration
	onedspec.standard.bandsep = 20.0
	# the star_name is the same as the root name of the dat file in the salt standards caldir
	onedspec.standard.star_name = starName
	# If running task individually, pop up epar window to let users change parameters if they want.
	if customRun == True:
		onedspec.standard.interact='YES' # give option to run interactively
		while True: 
			eparAnswer = raw_input("Do you want to further edit the parameters? 0 for no, 1 for yes.")
			if eparAnswer == '0' or eparAnswer == '1':
				break
			else:
				print "Invalid input: you must enter either 0 (no) or 1 (yes)."		
		if eparAnswer == '1':		
			iraf.eparam(onedspec.standard)
	# this runs onedspec.standard resulting in an std file named outputname 
	onedspec.standard(input=standardimage+'[1]',output=outputname)
	

# This function sets the parameters for the pyraf task onedspec.sensfunc and then runs it.
def run_sensfunc(stddata,sensname,customRun=False):
	# this resets the parameters of onedspec.sensfunc
	onedspec.sensfunc.unlearn()
	onedspec.sensfunc.apertures='1' # use only the first aperture (object data)
	onedspec.sensfunc.function='spline3'
	onedspec.sensfunc.order=2
	onedspec.sensfunc.extinction=''
	onedspec.sensfunc.ignoreaps='yes' # create only one sensfunc file
	onedspec.sensfunc.interactive='NO' # non-interactive fitting 
	onedspec.sensfunc.graphs='sri'
	onedspec.sensfunc.answer='NO' # user shouldn't have to input anything else, just check fit and quit
	# If running task individually, pop up epar window to let users change parameters if they want.
	if customRun == True:
		onedspec.sensfunc.interactive='YES' # give option to run interactively
		while True: 
			eparAnswer = raw_input("Do you want to further edit the parameters? 0 for no, 1 for yes.")
			if eparAnswer == '0' or eparAnswer == '1':
				break
			else:
				print "Invalid input: you must enter either 0 (no) or 1 (yes)."		
		if eparAnswer == '1':		
			iraf.eparam(onedspec.sensfunc)
	# runs onedspec.sensfunc resulting in a file named sensname
	onedspec.sensfunc(standards=stddata,sensitivity=sensname)
	

# This function sets the parameters for the pyraf task onedspec.calibrate and then runs it.
def run_calibrate(scienceimage,fluximage,sensfilename,customRun=False):
	# this resets the parameters of onedspec.calibrate
	onedspec.calibrate.unlearn()
	onedspec.calibrate.extinct='no'
	onedspec.calibrate.extinction=''
	onedspec.calibrate.sensitivity=sensfilename # the sensitivity function for flux calibration
	onedspec.calibrate.ignoreaps='yes' # look for and use sens*, not sens*0001, sens*0002, etc.
	# This stores some header information from the science image
	scihdr = pyfits.getheader(scienceimage,1)
	exptime = scihdr['EXPTIME']
	airmass = abs(scihdr['AIRMASS'])
	# this continues to set the necessary parameters of onedspec.standard
	onedspec.calibrate.airmass=airmass
	onedspec.calibrate.exptime=exptime
	
	# runs onedspec.calibrate resulting in a flux-calibrated spectrum with name fluximage (on [SCI] extension)
	scienceimage1 = scienceimage+'[1]'
	
	# If running task individually, pop up epar window to let users change parameters if they want.
	if customRun == True:
		while True: 
			eparAnswer = raw_input("Do you want to further edit the parameters? 0 for no, 1 for yes.")
			if eparAnswer == '0' or eparAnswer == '1':
				break
			else:
				print "Invalid input: you must enter either 0 (no) or 1 (yes)."		
		if eparAnswer == '1':		
			iraf.eparam(onedspec.calibrate)
	
	onedspec.calibrate(input=scienceimage1,output=fluximage)
	

# This function sets the parameters for the pyraf task onedspec.odcombine and then runs it.
def run_odcombine(inputlistseq,finalspectrumname,saltgain,customRun=False):
	# this resets the parameters of odcombine
	onedspec.odcombine.unlearn()
	onedspec.odcombine.group='all'
	onedspec.odcombine.combine='average'
	onedspec.odcombine.reject='avsigclip'
	onedspec.odcombine.apertures = 1
	onedspec.odcombine.outtype = 'real'
	onedspec.odcombine.smaskformat = 'bpmspectrum' # each input image has 'BPM' keyword linking it to its bpm
	onedspec.odcombine.smasktype = 'goodvalue'
	onedspec.odcombine.smaskvalue = 0.0
	onedspec.odcombine.blank = 0.0
	onedspec.odcombine.lsigma = 2.0
	onedspec.odcombine.hsigma = 2.0
	onedspec.odcombine.gain = saltgain
	onedspec.odcombine.masktype = 'goodvalue'
	onedspec.odcombine.maskvalue = 0.0
	# If running task individually, pop up epar window to let users change parameters if they want.
	if customRun == True:
		while True: 
			eparAnswer = raw_input("Do you want to further edit the parameters? 0 for no, 1 for yes.")
			if eparAnswer == '0' or eparAnswer == '1':
				break
			else:
				print "Invalid input: you must enter either 0 (no) or 1 (yes)."		
		if eparAnswer == '1':
			iraf.eparam(onedspec.odcombine)
	# running onedspec.odcombine resulting in the final wave- and flux-calibrated spectrum named finalspectrumname
	onedspec.odcombine(input=inputlistseq,output=finalspectrumname)
	

# This function creates a copy of a FITS file and preserves the 2-extension structure.
def run_imcopy(inputname,outputname,numExt=2,customRun=False):	
	if numExt == 1:
		imutil.imcopy(input=inputname+'[0]',output=outputname+'[append]')
	else:
		imutil.imcopy(input=inputname+'[0]',output=outputname+'[append]')
		imutil.imcopy(input=inputname+'[SCI]',output=outputname+'[append]')
		
		
# This function deletes a file.
def run_imdel(inputname,customRun=False):
	imutil.imdel(inputname)
	

# This function extracts bands from a spectrum stored in a FITS file.
def run_scopy(inputname,outputname,band,customRun=False):
	onedspec.scopy.unlearn()
	onedspec.scopy(input=inputname+'[1]',output=outputname,bands=band,format='multispec',clobber='yes',verbose='yes')

# This *generalized* function sets the parameters for the pyraf task imutil.imarith and then runs it.
def run_imarithGeneral(op1,op2,oper,outname,customRun=False):

	# This resets the parameters of imutil.imarith.
	imutil.imarith.unlearn()
	imutil.imarith.divzero=0.0
	# appends [1] to end of dividend image since that's where data lives
	# If running task individually, pop up epar window to let users change parameters if they want.
	if customRun == True:
		while True: 
			eparAnswer = raw_input("Do you want to further edit the parameters? 0 for no, 1 for yes.")
			if eparAnswer == '0' or eparAnswer == '1':
				break
			else:
				print "Invalid input: you must enter either 0 (no) or 1 (yes)."		
		if eparAnswer == '1':
			iraf.eparam(imutil.imarith)
	# This runs imutil.imarith
	imutil.imarith(operand1=op1,op=oper,operand2=op2,result=outname)	
	
	
	
	
			
>>>>>>> origin/master
