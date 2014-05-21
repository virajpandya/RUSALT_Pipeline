#!/usr/bin/env python
'''
Created on Mar 16, 2014

@author: cmccully
'''
import sys,os,shutil
from glob import glob
from pyraf import iraf
import pyfits
from numpy import array, unique
iraf.pysalt()
iraf.saltred()
iraf.saltspec()
iraf.twodspec()
iraf.longslit()

#Define the stages: Pysalt, LAcosmicx, Flatcombine, flatten, 'mosaic', 2D Identify, Rectify, Extract (Sci and Arcs), 1D Identify, 
#if standard, calculate flux calibration, Flux calibrate, if standard calculate telluric, telluric correction
allstages = ['pysalt','makeflats','flatten','mosaic','combine2d','identify2d', 'rectify', 'background','lax','extract','identify1d','dispcor1d'
            'stdsensfunc','fluxscale','speccombine','stdtelluric','telluric']
    
def run(dostages='all',stdstar = False,interactive = False):
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
    
    if 'pysalt' in stages: run_pysalt()
    
    if 'makeflats' in stages: run_makeflats()
    
    if 'flatten' in stages: run_flatten()
    
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
    
def run_pysalt():
    #Run the pysalt pipeline on the raw data.
    fs = glob('P*.fits')
    if len(fs)==0:
        print "There are no files to run the PySALT pipeline on."
        return
    
    #Copy the raw files into a raw directory
    if not os.path.exists('raw'): os.mkdir('raw')
    if not os.path.exists('work'): os.mkdir('work')
    for f in fs: shutil.copy(f,'raw/')
    for f in fs: shutil.move(f, 'work/')
    os.chdir('work')
    #Run each of the pysalt pipeline steps
    
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
    
    os.chdir('..')
    #Hold off on the the mosaic step for now. We want to do some processing on the individual chips
    

def get_flats(fs):
    #define a list to hold all of the flat filename and the gr angles
    flats = []
    grangles = []
    for f in fs:
        if pyfits.getval(f,'OBSTYPE')=='FLAT': 
            flats.append(f)
            grangles.append(pyfits.getval(f,'GR-ANGLE'))
    return array(flats),array(grangles)

def run_makeflats():
    os.chdir('work')
    fs = glob('bxgp*.fits')
    if len(fs)==0:
        print "There are no files to run flat combine."
        return
    #Figure out which images are flats
    #Figure out which grating angles were used
    allflats, grangles = get_flats(fs)
    #For each grating angle
    for ga in unique(grangles):
        #grab the flats for this gr angle
        flats = allflats[grangles == ga]
        
        #For each chip
        for i in range(1,7):
            #run imcombine with average and crreject, weighted by exposure time
            flatlist = ''
            for f in flats: 
                flatlist += './%s[%i],'%(f,i)
                #Add the exptime keyword to extension
                pyfits.setval(f,'EXPTIME',ext = i, value = pyfits.getval(f,'EXPTIME'))
            
            #set the output combined file name
            outname = 'fltcombine%0.2fc%i.fits'%(ga,i)
            if os.path.exists(outname): os.remove(outname)
            iraf.unlearn(iraf.imcombine)
            #don't forget to remove the last comma in the filelist
            iraf.flpr()
            iraf.imcombine(input = flatlist[:-1], output = outname,
                           combine='average',reject='crreject',weight='exposure',expname='EXPTIME')
            
####FIX ME           
            #Don't know if this is going to help, but it might.
            #We want to make an illumination correction file using surf fit
            
            #divide out the illumination correction before running response
            
            #Then make a bad pixel mask for the combined flat field flagging everything below some threshold.
            #We need to make sure we get the normalization of the response correct as to not screw up the fringing
            #It looks like the response uses the mean of the image to set the normalization. That's all well and good 
            #if there are no dark sections (shadows) and the flats are illumination corrected.
            
######
            #run iraf's response on the combined image
            #set the dispaxis keyword so that response knows which direction to fit
            pyfits.setval(outname,'DISPAXIS',value=1)
            iraf.unlearn(iraf.response)
            respoutname = 'flt%0.2fc%i.fits'%(ga,i)
            if os.path.exists(respoutname): os.remove(respoutname)
            iraf.response(calibration=outname, normalization=outname, response=respoutname,
                          interactive = True, function = 'spline3', order =31, low_reject=3.0,
                          high_reject=3.0,niterate=10, mode='hl' )
            
###FIX ME
            #After we correct for the response, we should probably rescale to the median of the image which should be less sensitive
            #to dark patches etc.

#####
            #remove the fltcombine file
            os.remove(outname)
    
    #clean up the flat files
    for f in allflats:os.remove(f)
    os.chdir('..')
    
def run_flatten():
    #use the star to grab both lax files and regular.
    fs = glob('*bxgp*.fits')
    if len(fs)==0:
        print "There are no files to flatten."
        return
    #Make sure there are science images or arcs
    #Figure out which grating angles were used
    #For each grating angle
    #For each science and arc image at that grating angle
    #For each chip
    #divide out the illumination correction and the flatfield, make sure divzero = 0.0
    #save the updated file
    
def run_mosaic():
    #Grab the flattened images.
    #run saltmosaic
    return

def run_combine2d():
    #Grab the mosaiced images
    #Find the science files and arcs.
    #All of the grating angles
    #for each grating angle
    #If there are more than one science image at the grating angle
    # combine the images using imcombine, sum, crreject, check to make sure the exptime keyword is updated
    return

def run_identify2d():
    #For each grating angle
    #Find all of the arcs.
    #For each arc:
    #run pysalt specidentify.
    return

def run_rectify():
    #For each grating angle
    #For each arc and science image
    #run specrectify
    return

def run_background():
    #For each rectified science image
    #Run iraf background 
    #Save the sky image by taking the difference in the original image and the background subtracted image
    return

def run_lax():
    fs = glob('bxgp*.fits')
    if len(fs)==0:
        print "There are no files to run LaCosmicX on."
        return
    #Figure out which files are science files
    #For each science file
    #open the file in read only mode.
    #For each chip
    #run lacosmicx on a copy of the data
    #replace the data with cleaned lax data
    #save the updated file
    #cleanup
    

def run_extract():
    #For each science image
    #run apall
    #run apsum on the corresponding arc image
    #run apsum on the corresponding sky image.
    return

def run_identify1d():
    #For each one d arc extraction
    #run iraf identify on the spectrum
    return

def run_dispcor1d():
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
    