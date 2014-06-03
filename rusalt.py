#!/usr/bin/env python
'''
Created on Mar 16, 2014

@author: cmccully
'''
import sys,os,shutil
from glob import glob
from pyraf import iraf
import pyfits
from numpy import array, unique, median, std, append
import fnmatch
iraf.pysalt()
iraf.saltred()
iraf.saltspec()
iraf.twodspec()
iraf.longslit()

# System-specific paths for standards, linelists, etc.
standardsPath = '/usr/local/astro64/iraf/extern/pysalt/data/standards/spectroscopic/'
lineListPath = '/usr/local/astro64/iraf/extern/pysalt/data/linelists/'

#Define the stages: Pysalt, LAcosmicx, Flatcombine, flatten, 'mosaic', 2D Identify, Rectify, Extract (Sci and Arcs), 1D Identify, 
#if standard, calculate flux calibration, Flux calibrate, if standard calculate telluric, telluric correction
allstages = ['pysalt','makeflats','flatten','mosaic','combine2d','identify2d', 'rectify', 'background','lax','extract','identify1d','dispcor1d'
            'stdsensfunc','fluxscale','speccombine','stdtelluric','telluric']

def globr(files):
    matches = []
    for root, dirnames, filenames in os.walk('./'):
        for filename in fnmatch.filter(filenames, files):
            matches.append(os.path.join(root, filename))
    return matches

def tofits(filename, data, hdr=None,clobber=False):
    """simple pyfits wrapper to make saving fits files easier."""
    from pyfits import PrimaryHDU,HDUList
    hdu = PrimaryHDU(data)
    if not hdr is None: hdu.header = hdr
    hdulist = HDUList([hdu])
    hdulist.writeto(filename, clobber=clobber,output_verify='ignore')
    
# preserves AND modifies old FITS headers (unlike IRAF tasks) and maintains PySALT file structure
def preservefits(oldfilename,newfilename,auxfilename,keys=[],newheader=False,data_ext=1):
    shutil.copyfile(oldfilename,newfilename)
    if auxfilename != '':
	    hduaux = pyfits.open(auxfilename)
	    dataaux = hduaux[0].data.copy() # data from auxiliary output by IRAF task
	    hdraux = hduaux[0].header.copy() # header from aux output by IRAF task, good for apall/1d-identify, etc.
	    hduaux.close()
	    hdunew = pyfits.open(newfilename,mode='update')
	    hdunew[data_ext].data = dataaux # simply replace the old data
	    if newheader == True: # replace the header of data_ext with the aux header
		    hdraux.update('EXTNAME','SCI') # IRAF doesn't maintain these next two keywords
		    hdraux.update('EXTVER',data_ext)
		    hdunew[data_ext].header = hdraux
	    for k in keys:
		    (hdunew[0].header)[k] = hdraux[k] # e.g., avg EXPTIME
		    (hdunew[data_ext].header)[k] = hdraux[k]
	    hdunew.flush(output_verify='ignore')
	    hdunew.close()
	    os.remove(auxfilename)
		
def tods9(filename):
    d = ds9.ds9(start=True)
    d.set('file '+filename)
    d.set('zoom to fit')
    d.set('zscale')
    print filename+' has been briefly opened in ds9.'
    
def run(dostages='all',stdstar = False,interactive = False, files = None):
    #Make sure the stages parameters makes sense
    try:  
        if dostages =='all': 
            n0 = 0
            n = len(allstages)
        elif '-' in dostages:
            n0 = allstages.index(dostages.split('-')[0]) 
            n = allstages.index(dostages.split('-')[1]) 
        elif '+' in dostages:
            n0 = allstages.index(dostages[:-1]) 
            n = len(allstages)
        else:
            n0 = allstages.index(dostages) 
            n = allstages.index(dostages) 
    
    except: print "Please choose a valid stage."
    

    stages = allstages[n0:n+1]
    print('Doing the following stages:')
    print(stages)
   
    if ',' in dostages: stages = dostages.split(',')
    
    if 'pysalt' in stages: run_pysalt(files)
    
    if 'makeflats' in stages: run_makeflats(files)
    
    if 'flatten' in stages: run_flatten(files)
    
    if 'mosaic' in stages: run_mosaic()
    
    if 'identify2d' in stages: run_identify2d()
    
    if 'rectify' in stages: run_rectify()
    
    if 'background' in stages: run_background()
    
    if 'lax' in stages: run_lax()
    
    if 'combine2d' in stages: run_combine2d()
    
    if 'extract' in stages: run_extract()
    
    if 'identify1d' in stages: run_identify1d()
    
    if 'dispcor1d' in stages: run_dispcor1d()
    
    if stdstar and 'stdsensfunc' in stages: run_stdsensfunc()
    
    if 'fluxcal' in stages: run_fluxcal()
    
    if 'speccombine' in stages: run_speccombine()
    
    if stdstar and 'stdtelluric' in stages: run_stdtelluric()
    
    if 'telluric' in stages: run_telluric()
    
def run_pysalt(fs=None):
    #Run the pysalt pipeline on the raw data.
    if fs is None: fs = glob('P*.fits')
    if len(fs)==0:
        print "ERROR: No raw files to run PySALT pre-processing."
        return
    
    #Copy the raw files into a raw directory
    if not os.path.exists('raw'): os.mkdir('raw')
    if not os.path.exists('work'): os.mkdir('work')
    for f in fs: shutil.copy(f,'raw/')
    for f in fs: shutil.move(f, 'work/')
    os.chdir('work')
    #Run each of the pysalt pipeline steps deleting temporary files as we go
    
    #saltprepare
    iraf.unlearn(iraf.saltprepare)
    #Currently, there is not a bad pixel mask provided by SALT so we don't create one here.
    iraf.saltprepare(images = 'P*.fits',clobber=True,mode='h')
    
    for f in glob('P*.fits'): os.remove(f)
    #saltgain
    iraf.unlearn(iraf.saltgain)
    #Multiply by the gain so that everything is in electrons.
    iraf.saltgain(images='pP*.fits',gaindb = os.environ['PYSALT']+'/data/rss/RSSamps.dat', mult=True, usedb=True,mode='h')
    for f in glob('pP*.fits'): os.remove(f)
    
    #write a keyword in the header keyword gain = 1 in each amplifier
    fs = glob('gpP*.fits')
    for f in fs: 
        for i in range(1,7): pyfits.setval(f, 'GAIN', ext= i,value=1.0)

    #saltxtalk
    iraf.unlearn(iraf.saltxtalk)
    iraf.saltxtalk(images = 'gpP*.fits',clobber=True, usedb=True, xtalkfile = os.environ['PYSALT']+'/data/rss/RSSxtalk.dat',mode='h')
    for f in glob('gpP*.fits'): os.remove(f)
    
    #saltbias
    iraf.unlearn(iraf.saltbias)
    iraf.saltbias(images='xgpP*.fits',clobber=True,mode='h')
    for f in glob('xgpP*.fits'): os.remove(f)
    
    #Put all of the newly created files into the pysalt directory
    os.mkdir('pysalt')
    for f in glob('bxgP*.fits'): shutil.move(f,'pysalt')
    os.chdir('..')
    #Hold off on the the mosaic step for now. We want to do some processing on the individual chips
    
    return fs

def get_ims(fs,imtype):
    typekeys ={'sci':'OBJECT','arc':'ARC','flat':'FLAT'}
    ims = []
    grangles = []
    for f in fs:
        if pyfits.getval(f,'OBSTYPE')==typekeys[imtype]: 
            flats.append(f)
            grangles.append(pyfits.getval(f,'GR-ANGLE'))
    return array(flats),array(grangles)

def run_makeflats(fs=None):
    #Note the list of files need to not include any paths relative to the work directory
    #Maybe not the greatest convention, but we can update this later
    os.chdir('work')
    if fs is None:  
        fs = glob('pysalt/bxgp*.fits')
    if len(fs)==0:
        print "ERROR: No flat-fields to combine and normalize."
        #Maybe we need to change folders to fail gracefully, but I don't have this here for now.
        return
    #make a flats directory
    os.mkdir('flats')
    
    #Figure out which images are flats
    #Figure out which grating angles were used
    allflats, grangles = get_ims(fs,'flat')
    
    os.mkdir() 
    #For each grating angle
    for ga in unique(grangles):
        #grab the flats for this gr angle
        flats = allflats[grangles == ga]
        
        #For each chip
        for c in range(1,7):
            #run imcombine with average and crreject, weighted by exposure time
            flatlist = ''
            for f in flats: 
                flatlist += '%s[%i],'%(f,c)
                #Add the exptime keyword to extension
                pyfits.setval(f,'EXPTIME',ext = c, value = pyfits.getval(f,'EXPTIME'))
            
            #set the output combined file name
            combineoutname = 'flats/flt%0.2fcomc%i.fits'%(ga,c)
            if os.path.exists(combineoutname): os.remove(combineoutname)
            iraf.unlearn(iraf.imcombine); iraf.flpr()
            
            #don't forget to remove the last comma in the filelist
            iraf.imcombine(input = flatlist[:-1], output = combineoutname, combine='average',
                           reject='sigclip',lsigma = 3.0, hsigma = 3.0,
                           weight='exposure',expname='EXPTIME')
                      
            #Don't know if this is going to help, but it might.
            #We want to make an illumination correction file before running response:
            illumoutname = 'flats/flt%0.2fillc%i.fits'%(ga,c) 
            iraf.unlearn(iraf.illumination); iraf.flpr()
            iraf.illumination(images=combineoutname,illuminations=illumoutname,interactive=False, 
                              naverage=-40, order =11, low_reject=3.0,high_reject=3.0,niterate=5,mode='hl')
                         
            #Flag any pixels in the illumination correction< 0.1
            illumhdu = pyfits.open(illumoutname,mode='update')
            illumhdu[0].data[illumhdu[0].data <= 0.1] = 0.0
            illumhdu.flush()
                        
            #Get 40 pixels out of the middle of the image and median them to run response
            combinehdu = pyfits.open(combineoutname)
            ny = combinehdu[0].data.shape[0]
            #divide out the illumination correction before running response
            flat1d = median(combinehdu[0].data[ny/2 - 21: ny/2 +20,:].copy() / illumhdu[0].data[ny/2 - 21: ny/2 +20,:].copy() ,axis=0)
            #close the illumination file because we don't need it anymore
            illumhdu.close()
            
            flat1dfname = 'flats/flt%0.2f1dmc%i.fits'%(ga,c)
            tofits(flat1dfname,flat1d, hdr=combinehdu[0].header.copy())
            
            #run response
            resp1dfname = 'flats/flt%0.2f1drc%i.fits'%(ga,c)
            iraf.response(flat1dfname,flat1dfname, resp1dfname ,order = 31, 
                          interactive=False,  naverage=-5, low_reject=3.0,high_reject=3.0,
                          niterate=5, mode='hl'  )

            resp1dhdu = pyfits.open(resp1dfname)
            resp1d =resp1dhdu[0].data.copy()
            resp1dhdu.close()
            
            #After response divide out the response function
            #normailze the 1d resp to its median
            resp1d/= median(resp1d)
            
            #Chuck any outliers
            flatsig = std(resp1d - 1.0)
            resp1d[abs(resp1d - 1.0) > 5.0 * flatsig] = 1.0
            resp = flat1d/resp1d
            
            resp2dfname = 'flats/flt%0.2resc%i.fits'%(ga,c)
            resp2d = combinehdu[0].data.copy()/resp
            tofits(resp2dfname,resp2d,hdr = combinehdu[0].header.copy())
            combinehdu.close()
            
            #close the combined flat because we don't need it anymore
            combinehdu.close()
            
            pyfits.setval(resp2dfname,'DISPAXIS',value=1)
            
            #Reset any pixels in the flat field correction< 0.1
            #We could flag bad pixels here if we want, but not right now
            flathdu = pyfits.open(resp2dfname,mode='update')
            flathdu[0].data[resp2dfname[0].data <= 0.1] = 0.0
            flathdu.flush()
            flathdu.close()
    #Step back up to the top directory
    os.chdir('..')
    
def run_flatten(fs = None):
    if fs is None: fs = glob('pysalt/bxgpP*.fits')
    if len(fs)==0:
        print "ERROR: No images to flat-field."
        return
    if not os.path.exists('flts'): os.mkdir('flts')
    #Make sure there are science images or arcs and what grating angles were used
    scifs, scigas = get_ims(fs,'sci')
    arcfs, arcgas = get_ims(fs,'arc')

    ims = append(scifs,arcfs)
    gas = append(scigas, arcgas)
    #For each science and arc image
    for i,f in enumerate(ims):
        thishdu = pyfits.open(f)
        ga = gas[i]
        #For each chip
        for c in range(1,7):
            #open the corresponding response file
            resphdu = pyfits.open('flats/flt%0.2fresc%i.fits'%(ga,c))
            #divide out the illumination correction and the flatfield, make sure divzero = 0.0
            thishdu[c].data /= resphdu.copy()
            resphdu.close()
            
        #save the updated file
        if f in scifs: typestr = 'sci'
        else: typestr = 'arc'
        #get the image number
        #by salt naming convention, these should be the last 3 characters before the '.fits'
        imnum = f[-8:-5]
        outname = 'flts/'+typestr+'%0.2fflt%03i.fits'%(ga,imnum)
        thishdu.writeto(outname)
        thishdu.close()
        
def run_mosaic(fs=None):
    # If the file list is not given, grab the default files
    if fs is None: fs = glob('flts/*.fits') 
    #Abort if there are no files
    if len(fs)==0:
        print "ERROR: No flat-fielded images to mosaic."
        return
    
    #Get the images to work with
    scifs, scigas = get_ims(fs,'sci')
    arcfs, arcgas = get_ims(fs,'arc')
    ims = append(scifs,arcfs)
    gas = append(scigas, arcgas)
    
    if not os.path.exists('mos'):os.mkdir('mos')
    for i,f in enumerate(ims): 
            fname =  f.split('/')[1]
            typestr = fname[:3]
            #by our naming convention, these should be the last 3 characters before the '.fits'
            imnum = fname[-8:-5]
            outname = 'mos/'+typestr+'%0.2fmos%03i.fits'%(ga,imnum)
            iraf.unlearn(iraf.saltmosaic) ; iraf.flpr()# prepare to run saltmosaic
            iraf.saltmosaic(images=f,outimages=outname,outpref='',clobber=True,mode='h') 

def run_combine2d(fs=None):
    #Grab the mosaiced images
    if fs is None: fs = glob('*mos*.fits') 
    if len(fs)==0:
        print "There are no mosaiced images to combine."
        return
    #Find the science files and arcs (and their grating angles)
    (scifs,scigas),(arcfs,arcgas) = get_scis(fs),get_arcs(fs)
    for groupfs,groupgas in ((scifs,scigas),(arcfs,arcgas)): # do same thing for scis and arcs
        for i,f in enumerate(groupgas): 
            if len(groupgas) != len(unique(groupgas)): # this means there are multiple files for a given gr-angle
                repgas = set([a for a in groupgas if groupgas.count(a)>1]) # repeated gr-angles
                if f in arcfs: 
                    print "ERROR: multiple arcs found for gr-angles",repgas
                    return # user should inspect+delete extra arcs then re-run
                else:
                    for ga in repgas:
	                    suffix = ''
	                    for j in glob('sci%0.2fmos*.fits'%(ga)):
		                    suffix += j[-8:-5] # careful in the future: longer imgnum sequences now possible
	                    iraf.unlearn(iraf.imcombine); iraf.flpr()
	                    iraf.imcombine(input = 'sci%0.2fmos*.fits'%(ga)+'[1]',output='auxstk.fits',combine='sum',
                               	 	   reject='crreject',lsigma = 3.0, hsigma = 3.0,
                                       weight='exposure',expname='EXPTIME') 
                        # preserve file structure and summed EXPTIME keyword of stacked data
                        outname = 'sci%0.2fstk'%(ga)+suffix+'.fits'
                        preservefits(oldfilename=glob('sci%0.2fmos*.fits'%(ga))[0],newfilename=outname,
                        		     keys=['EXPTIME'],auxfilename='auxstk.fits')
				        		     
##### Need to find an efficient way to combine important header keys of individual images (e.g., date/time-obs in ext-0-header)
##### I think imcombine has an option for this that we could use for the ext-1-header
##### we could maybe just replace the old ext-1-header with the imcombine header
##### and use the keys arg in preservefits to change our most important keys in the 0-ext-header
##### Also, make sure combine='sum' & weight='exposure' ==> exptime_total = exptime1+exptime2 (sum, not average) -VP
	# combine the images using imcombine, sum, crreject, check to make sure the exptime keyword is updated - CM


def run_identify2d(fs=None):
    if fs is None: fs = glob('mos/arc*mos*.fits') 
    if len(fs)==0:
        print "ERROR: No mosaiced arcs for PySALT's (2D) specidentify."
        return
    (arcfs,arcgas) = get_arcs(fs)
    if not os.path.exists('sol'):os.mkdir('sol')
    for i,f in enumerate(arcfs):
	    # new way of finding imgnum based on naming convention for possibly stacked images
	    ga,imgnum = arcgas[i],f[11:f.index('.fits')]
	    # find lamp and corresponding linelist
	    lamp = pyfits.getval(f,'LAMPID')
	    if lamp == 'Th Ar':
		    print 'the lamp is '+lamp+' for '+f 
		    lamplines = lineListPath+'ThAr.salt' # global variable lineListPath defined in beginning
	    elif lamp == 'Xe':
		    print 'the lamp is '+lamp+' for '+f
		    lamplines = lineListPath+'Xe.txt'
	    elif lamp == 'Ne':
		    print 'the lamp is '+lamp+' for '+f
		    lamplines = lineListPath+'NeAr.salt' 
	    elif lamp == 'Cu Ar':
		    print 'the lamp is '+lamp+' for '+f
		    lamplines = lineListPath+'CuAr.txt'
	    elif lamp == 'Ar':
		    print 'the lamp is '+lamp+' for '+f
		    lamplines = lineListPath+'Argon_hires.salt'
	    elif lamp == 'Hg Ar':
		    print 'the lamp is '+lamp+' for '+f
		    lamplines = lineListPath+'HgAr.txt'
	    else:
		    print 'Could not find the proper linelist for '+lamp+' lamp.'
		    return
	    # run pysalt specidentify
	    idfile = 'sol/arc%0.2fsol'%(ga)+'.fits' # no need for imgnum complications
	    iraf.unlearn(iraf.specidentify)
	    iraf.specidentify(images=f,linelist=lamplines,outfile=idfile,guesstype='rss',automethod='Matchlines',
					      function='legendre',order=3,rstep=100,rstart='middlerow',mdiff=5,inter='yes',
					      startext=1,clobber='yes',verbose='yes')

def run_rectify(fs=None):
    if fs is None: fs = glob('sol/arc*sol*.fits') 
    if len(fs)==0:
        print "ERROR: No wavelength solutions for rectification."
        return
    # rectify each sci/arc and pop it open briefly in ds9
    (scifs,scigas),(arcfs,arcgas) = get_scis(glob('*mos*.fits')),get_arcs(glob('*mos*.fits'))
    if not os.path.exists('rec'):os.mkdir('rec')
    for groupfs,groupgas in ((scifs,scigas),(arcfs,arcgas)): # do same thing for scis and arcs
	    for i,f in enumerate(groupfs): 
		    if f in scifs: typestr = 'sci'
		    else: typestr = 'arc'
		    ga,imgnum = groupgas[i],f[11:f.index('.fits')]
		    outfile = 'rec/'+typestr+'%0.2frec'%(ga)+imgnum+'.fits'
		    iraf.unlearn(iraf.specrectify)
		    iraf.specrectify(images=f,outimages=outfile,solfile='arc%0.2fsol'%(ga)+'.fits',outpref='',caltype='line',
						     function='legendre',order=3,inttype='interp',clobber='yes',verbose='yes') 	   		  
		    tods9(outfile)

    # I think this step should just be combined with the identify2d() step, 
    # or at least the arc rectifications to verify good sol -VP

def split_by_chip(fs=None):
   if fs is None: fs = glob('rec/*rec*.fits') 
   if len(fs)==0:
       print "ERROR: No rectified images to split by chip."
       return
    # Grab the rectified science images (only they need to be split up for bkg+lax)
    (scifs,scigas) = get_scis(fs)
    # Pixel begin:end numbers (+/- epsilon included) for chip gaps based on different CCDSUM header key values (binning)
    chipGapPix = {'2 2':((1010,1095),(2085,2160)),'2 4':((1035,1121),(2115,2191)),'4 4':((500,551),(1035,1091))}
    # Split each image into 3 imgs such that the middle img has both chip gaps, and the other 2 imgs don't have chip gaps
    for i,f in enumerate(scifs):
        ga,imgnum = scigas[i],f[11:f.index('.fits')]
        (c1min,c1max),(c2min,c2max) = chipGapPix[pyfits.getval(f,'CCDSUM')][0],chipGapPix[pyfits.getval(f,'CCDSUM')][1]
        hdu = pyfits.open(f)
        data = hdu[1].data.copy()
        hdr1 = hdu[1].header.copy()
        for c in range(1,4):
            if c == 1: tofits(f[:-5]+'c1.fits',data[0:c1min],hdr=hdr1) # same rec directory, just adding c# before '.fits'
            elif c == 2: tofits(f[:-5]+'c2.fits',data[c1min:c2max+1],hdr=hdr1)    
            elif c == 3: tofits(f[:-5]+'c3.fits',data[c2max+1:],hdr=hdr1)               
    
def run_background(fs=None):
   if fs is None: fs = glob('rec/*rec*c*.fits') 
   if len(fs)==0:
       print "ERROR: No rectified chip-based images for 2D-background-subtraction."
       return
    # Get rectified science images and gr-angles
    (scifs,scigas) = get_scis(fs)
    if not os.path.exists('bkg'):os.mkdir('bkg')
    for i,f in enumerate(scifs):
        # Run automated 2D background subtraction on each individual chip image
        pyfits.setval(f,'DISPAXIS',value=1) # just in case since this is auto background
        hdu = pyfits.open(f)
        # the outfile name is very similar, just change folder prefix and 3-char stage substring
        outfile = 'bkg/'+f[4:12]+'bkg'+f[15:]
        iraf.background(input=f,output='auxbkg.fits',interactive='no',naverage='-100',function='legendre',
                        order=2,low_reject=1.0,high_reject=1.0,niterate=10,grow=0.0)
        hduaux = pyfits.open('auxbkg.fits')
        hdu[0].data = hduaux[0].data.copy()
        hduaux.close()
        hdu.writeto(outfile) # saving the updated file (data changed)
        os.remove('auxbkg.fits')
        ###### modify/shorten preservefits() to do the above instead if possible, or use data_ext keyword, saves lines
    
def run_lax(fs=None):
    fs = glob('bkg/*bkg*.fits')
    if len(fs)==0:
        print "ERROR: No background-subtracted files for LaCosmicX."
        return
    
    #Figure out which files are science files
    #For each science file
    #open the file in read only mode.
    #For each chip
    #run lacosmicx on a copy of the data
    #replace the data with cleaned lax data
    #save the updated file
    #cleanup
    
def remosaic_chips(fs=None):
    # for each rectified (non-split) image, make a copy
    # use numpy to remosaic three lacosmicx-processed chips according to pixel #
    # replace data in copied rectified image with re-mosaiced chips
    return

def run_extract(fs=None):
    #For each science image
    #run apall
    #run apsum on the corresponding arc image
    #run apsum on the corresponding sky image.
    return

def run_identify1d(): # unnecessary at the moment given that we are doing pysalt specidentify manually -VP
    #For each one d arc extraction
    #run iraf identify on the spectrum
    return

def run_dispcor1d(): # unnecessary at the moment given that we are doing pysalt specidentify manually -VP
    #For each 1d arc 
    #run iraf dispcor on the science and sky images
    return

def run_stdsensfunc():
    return

def run_fluxcal():
    return 

def run_speccombine():
    return 

def run_stdtelluric():
    return 

def run_telluric():
    return 

if __name__ == '__main__':
    if sys.version_info[1] >= 7 or sys.version_info[0] > 2:
        # Parse the input arguments. You have to use a different package based on the python version: What a pain!
        import argparse
        parser = argparse.ArgumentParser(description='Open a ds9 window and step through the image in a systematic way.')
        parser.add_argument('filenames', metavar='filenames', nargs='*',
                   help='Filenames of the images to be loaded in the ds9 window.')
        args = parser.parse_args()
        run()
        sys.exit("Thanks for using this pipeline!")
    else:
        import optparse
        parser = optparse.OptionParser()
        parser.add_option("-z", dest="zoom",
                          help="Zoom factor to view the image.", metavar="zoom", default=2.0)
        (options, args) = parser.parse_args()
        run()
        sys.exit("Thanks for using this pipeline!")
    
