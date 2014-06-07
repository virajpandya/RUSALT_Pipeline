#!/usr/bin/env python
'''
Created on Mar 16, 2014

@author: cmccully
'''
import sys,os,shutil
from glob import glob
from pyraf import iraf
import pyfits
import numpy as np 
import fnmatch
iraf.pysalt()
iraf.saltred()
iraf.saltspec()
iraf.onedspec()
iraf.twodspec()
iraf.longslit()
iraf.apextract()
iraf.imutil()

# System-specific paths for standards, linelists, etc.
standardsPath = '/usr/local/astro64/iraf/extern/pysalt/data/standards/spectroscopic/'
lineListPath = '/usr/local/astro64/iraf/extern/pysalt/data/linelists/'
pipeStandardsPath = '/usr/local/astro64/rusaltsn/standards/'
skyLineSpectrum = '/Users/vgpandya/mastersky.fits' 

'''
There are two chip gaps, each spanning a range of columns (pixels). To avoid using data near the chip gap edges (primarily
due to skylines and the fact that we have data covering the chip gaps from overlapping spectra), we take the effective edges of
the chip gaps to be +/- ~20px beyond their definite edges. The begin:end pixel numbers of each chip gap will change if the binning
is different ('CCDSUM' header keyword), so the dict chipGapPix below allows one to access the pixel numbers based on 'CCDSUM'.
'''
chipGapPix = {'2 2':((1010,1095),(2085,2160)),'2 4':((1035,1121),(2115,2191)),'4 4':((500,551),(1035,1091))}

# Define a global variable (array) of possible standard star names 
possiblestds = glob(standardsPath+'m*.dat') # standardsPath is a global variable for PySALT standards directory
for n,s in enumerate(possiblestds): possiblestds[n] = (s.split('.'))[0][1:].lower() # remove m prefix & .fits suffix => std star name

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
    print 'DS9: '+filename+' has been briefly opened.'
    
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
        print "WARNING: No raw files to run PySALT pre-processing."
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
        if pyfits.getval(f,'OBSTYPE')==typekeys[imtype] or (pyfits.getval(f,'OBJECT')).lower() in possiblestds: 
            ims.append(f)
            grangles.append(pyfits.getval(f,'GR-ANGLE'))
    return np.array(ims),np.array(grangles)

def run_makeflats(fs=None):
    #Note the list of files need to not include any paths relative to the work directory
    #Maybe not the greatest convention, but we can update this later
    os.chdir('work')
    if fs is None:  
        fs = glob('pysalt/bxgp*.fits')
    if len(fs)==0:
        print "WARNING: No flat-fields to combine and normalize."
        #Maybe we need to change folders to fail gracefully, but I don't have this here for now.
        return
    #make a flats directory
    os.mkdir('flats')
    
    #Figure out which images are flats
    #Figure out which grating angles were used
    allflats, grangles = get_ims(fs,'flat')
    
    os.mkdir() 
    #For each grating angle
    for ga in np.unique(grangles):
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
            flat1d = np.median(combinehdu[0].data[ny/2 - 21: ny/2 +20,:].copy() / illumhdu[0].data[ny/2 - 21: ny/2 +20,:].copy() ,axis=0)
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
            resp1d/= np.median(resp1d)
            
            #Chuck any outliers
            flatsig = np.std(resp1d - 1.0)
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
        print "WARNING: No images to flat-field."
        return
    if not os.path.exists('flts'): os.mkdir('flts')
    #Make sure there are science images or arcs and what grating angles were used
    scifs, scigas = get_ims(fs,'sci')
    arcfs, arcgas = get_ims(fs,'arc')

    ims = np.append(scifs,arcfs)
    gas = np.append(scigas, arcgas)
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
        print "WARNING: No flat-fielded images to mosaic."
        return
    
    #Get the images to work with
    scifs, scigas = get_ims(fs,'sci')
    arcfs, arcgas = get_ims(fs,'arc')
    ims = np.append(scifs,arcfs)
    gas = np.append(scigas, arcgas)
    
    if not os.path.exists('mos'):os.mkdir('mos')
    for i,f in enumerate(ims): 
        ga = gas[i]
        fname =  f.split('/')[1]
        typestr = fname[:3]
        #by our naming convention, these should be the last 3 characters before the '.fits'
        imnum = fname[-8:-5]
        outname = 'mos/'+typestr+'%0.2fmos%03i.fits'%(ga,imnum)
        iraf.unlearn(iraf.saltmosaic) ; iraf.flpr()# prepare to run saltmosaic
        iraf.saltmosaic(images=f,outimages=outname,outpref='',clobber=True,mode='h') 


def run_identify2d(fs=None):
    if fs is None: fs = glob('mos/arc*mos*.fits') 
    if len(fs)==0:
        print "WARNING: No mosaiced arcs for PySALT's (2D) specidentify."
        return
    arcfs,arcgas = get_ims(fs,'arc')
    if not os.path.exists('id2'):os.mkdir('id2')
    
    lampfiles = {'Th Ar':'ThAr.salt','Xe':'Xe.txt', 'Ne':'NeAr.salt', 'Cu Ar':'CuAr.txt',
                 'Ar':'Argon_hires.salt', 'Hg Ar':'HgAr.txt' }
    for i,f in enumerate(arcfs):
        ga = arcgas[i]
        # find lamp and corresponding linelist
        lamp = pyfits.getval(f,'LAMPID')
        print 'the lamp is '+lamp+' for '+f 
        # linelistpath is a global variable defined in beginning, path to where the line lists are.
        try: lamplines = lineListPath+lampfiles[lamp]
        except KeyWARNING:
            print('Could not find the proper linelist for '+lamp+' lamp.')
            raise
        
        #img num should be right before the .fits
        imgnum = f[-8:-5]
        # run pysalt specidentify
        idfile = 'id2/arc%0.2fid2%03i'%(ga,imgnum)+'.db' 
        iraf.unlearn(iraf.specidentify); iraf.flpr()
        iraf.specidentify(images=f,linelist=lamplines,outfile=idfile,guesstype='rss',automethod='Matchlines',
        			      function='legendre',order=3,rstep=100,rstart='middlerow',mdiff=5,inter='yes',
        			      startext=1,clobber='yes',verbose='yes')

def run_rectify(fs=None):
    if fs is None: fs = glob('id2/arc*id2*.db') 
    if len(fs)==0:
        print "WARNING: No wavelength solutions for rectification."
        return
    # rectify each sci/arc and pop it open briefly in ds9
    scifs,scigas = get_ims(glob('mos/*mos*.fits'),'sci')
    arcfs,arcgas =get_ims(glob('mos/*mos*.fits'),'arc')
    
    ims = np.append(scifs,arcfs)
    gas = np.append(scigas, arcgas)
    
    if not os.path.exists('rec'):os.mkdir('rec')
    for i,f in enumerate(ims): 
            if f in scifs: typestr = 'sci'
            else: typestr = 'arc'
            ga,imgnum = gas[i],f[-8:-5]
            outfile = 'rec/'+typestr+'%0.2frec'%(ga)+imgnum+'.fits'
            iraf.unlearn(iraf.specrectify); iraf.flpr()
            iraf.specrectify(images=f,outimages=outfile,solfile='arc%0.2fsol'%(ga)+'.fits',outpref='',caltype='line',
						     function='legendre',order=3,inttype='interp',clobber='yes',verbose='yes') 	   		  

def run_unmosaic(fs=None):
    if fs is None: fs = glob('rec/*rec*.fits') 
    if len(fs)==0:
        print "WARNING: No rectified images to split by chip."
        return
    # Grab the rectified science images (only they need to be split up for bkg+lax)
    scifs,scigas = get_ims(fs,'sci')
    '''
    We're splitting by the 3 chips now instead of the 6 'amps' as in the pre-saltmosaic multi-FITS extension files from before for two 
    reasons. First, where the chip gaps begin and end is visible in the 2D spectra (counts == 0), but not for the amps since they are just
    based on the pre-saltmosaic images, e.g., the second amp begins directly after the first amp. More importantly, the rectification
    probably smeared out the column where one amp ends and the next begins (just look at how it curves the chip gaps).
    '''
    # Split each image into 3 imgs such that the middle img has both chip gaps, and the other 2 imgs don't have chip gaps
    for i,f in enumerate(scifs):
        ga,imgnum = scigas[i],f[11:f.index('.fits')]
        # Global variable chipGapPix returns the begin:end pixel numbers for the 2 chip gaps (depends on binning, i.e., CCDSUM).
        (c1min,c1max),(c2min,c2max) = chipGapPix[pyfits.getval(f,'CCDSUM')][0],chipGapPix[pyfits.getval(f,'CCDSUM')][1]
        # Copy 'GR-ANGLE' and 'OBSTYPE' from hdr0 to hdr1 for get_ims(); don't need to remove the keys since these split spectra
        # are only used in run_background() and run_lax(); we use the original *rec*.fits non-split spectra for apall.
        for k in ['GR-ANGLE','OBSTYPE']: pyfits.setval(f,ext=1,k,value=pyfits.getval(f,k))
        hdu = pyfits.open(f)
        data = hdu[1].data.copy()
        hdr1 = hdu[1].header.copy()
        for c in range(1,4):
            if c == 1: tofits(f[:-5]+'c1.fits',data[:,0:c1min],hdr=hdr1) # same rec directory, just adding c# before '.fits'
            elif c == 2: tofits(f[:-5]+'c2.fits',data[:,c1min:c2max+1],hdr=hdr1)    
            elif c == 3: tofits(f[:-5]+'c3.fits',data[:,c2max+1:],hdr=hdr1)               
    

def run_createbpm(fs=None):
    if fs is None: fs = glob('rec/*rec*c*.fits') 
    if len(fs)==0:
        print "WARNING: No rectified chip-based images for bad pixel mask creation."
        return
    ''' The BPM is used for lacosmicx '''
    # Get rectified science images and gr-angles (individual chips)
    scifs,scigas = get_ims(fs,'sci')
    if not os.path.exists('bpm'):os.mkdir('bpm')
    for i,f in enumerate(scifs):
        # the outfile name is very similar, just change folder prefix and 3-char stage substring
        outfile = 'bpm/'+f[4:12]+'bpm'+f[15:]
        iraf.unlearn(iraf.imexpr)
        iraf.imexpr(expr='a==0',output=outfile,a=f)
        # for record-keeping purposes
        pyfits.setval(f,'BPM',value=outfile)   
    

def run_background(fs=None):
    if fs is None: fs = glob('rec/*rec*c*.fits') 
    if len(fs)==0:
        print "WARNING: No rectified chip-based images for 2D-background-subtraction."
        return
    
    # Get rectified science images and gr-angles
    scifs,scigas = get_ims(fs,'sci')
    if not os.path.exists('bkg'):os.mkdir('bkg')
    for i,f in enumerate(scifs):
        # Run automated 2D background subtraction on each individual chip image
        pyfits.setval(f,'DISPAXIS',value=1) # just in case since this is automated 2D-bkg-sub
        hdu = pyfits.open(f)
        # the outfile name is very similar, just change folder prefix and 3-char stage substring
        outfile = 'bkg/'+f[4:12]+'bkg'+f[15:]
        iraf.unlearn(iraf.background)
        iraf.background(input=f,output='auxbkg.fits',interactive='no',naverage='-100',function='legendre',
                        order=2,low_reject=1.0,high_reject=1.0,niterate=10,grow=0.0)
        hduaux = pyfits.open('auxbkg.fits')
        hdu[0].data = hduaux[0].data.copy()
        hduaux.close()
        hdu.writeto(outfile) # saving the updated file (data changed)
        os.remove('auxbkg.fits')
        ###### modify/shorten preservefits() to do the above instead if possible, or use data_ext keyword, saves lines

    
def run_lax(fs=None):
    if fs is None: fs = glob('bkg/*bkg*.fits')
    if len(fs)==0:
        print "WARNING: No background-subtracted files for LaCosmicX."
        return
    
    # Get background-subtracted chip-based science images and gr-angles.
    scifs,scigas = get_ims(fs,'sci')
    if not os.path.exists('lax'):os.mkdir('lax')
    for i,f in enumerate(scifs):
        outname_img = 'lax/'+f[4:12]+'lax'+f[15:]
        outname_msk = 'lax/'+f[4:12]+'cpm'+f[15:] # cosmic (ray) pixel mask
        hdu = pyfits.open(f)
        datain = hdu[0].data.copy()
        # Get mode of image counts as a proxy for pre-sky-subtracted-level
        mode = float(imutil.imstat(images=f,fields='mode',format='no',Stdout=1)[0])
        # Get BPM data array; filename is similar but 3-char stage is bpm; remember to change folder prefix 'bkg/' to 'bpm/'
        namebpm = 'bpm/'+a.split('/')[1][0:8]+'bpm'+a.split('/')[1][11:]
        hdubpm = pyfits.open(namebpm)
        databpm = hdubpm[0].data.copy()
        hdubpm.close()
        # Run lacosmicx (assume that gain==1.0 because of iraf.pysalt.saltred.saltgain)
        ''' what might we do about RDNOISE? '''
        datalax = lacosmicx.run(inmat=datain,inmask=databpm,outmaskfile=outname_msk,
		                        sigclip=6.0,objlim=3.0,sigfrac=0.1,gain=1.0,pssl=mode,robust=True)
        # Update and save file
        hdu[0].data = datalax
        hdu.writeto(outname_img)
        hdu.close()
        
    
def run_remosaic(fs=None): 
    '''
    Our algorithm currently involves running apall on the non-2D-bkg-subtracted images
    so we can just use the non-run_unmosaic()-processed rectified images.
    But also, instead of using the lacosmicx output science images, we want to manually
    flag all cosmic ray pixels in the rectified image based on the lacosmicx cosmic ray pixel mask (CPM).
    So for now (unless we change the algorithm), we just want to mosaic the BPM and CPM images. - VP
    '''
    if fs is None: fs = glob('bpm/*bpm*.fits')
    if len(fs)==0:
        print "WARNING: No bad pixel masks to remosaic." # without BPM, wouldn't have made CPM anyway, so this check is sufficient
        return
    
    # Loop over the mosaiced rectified science images
    scifs,scigas = get_ims(glob('rec/*rec*.fits'),'sci') # THIS LOOP IS BAD, includes individual chips, maybe add len(name) check - VP
    for i,f in enumerate(scifs):
        ga,imgnum = scigas[i],f[11:f.index('.fits')]
        folder = {'bpm':'bpm/','cpm':'lax/'} # prefix for outfile names
        # Create BPM and CPM mosaics
        for masktype in ['bpm','cpm']:    
            outfile = folder[masktype]+f[4:12]+masktype+imgnum+'.fits' # gets rid of the 'c#' before '.fits'
            chips = []
            for c in range(1,4):
                maskfile = folder[masktype]+f[4:12]+masktype+imgnum+'c%i'%(c)+'.fits' # naming convention is very similar
                hdu = pyfits.open(maskfile)
                data = hdu[0].data.copy()
                chips.append(data)
                hdu.close()
            datamask = np.concatenate(chips,axis=1) # axis 1 is horizontal => the chips will be stacked horizontally
            tofits(outfile,datamask)
            # Add 'BPM' or 'CPM' to rectified image header for record-keeping purposes
            pyfits.setval(f,masktype.upper(),value=outfile)
        

def run_flagcosmics(fs=None):
    '''
    Instead of using the lacosmicx cleaned output image, we will flag cosmic ray pixels in the science data
    such that apall will reject those pixels when doing the variance-weighted extraction.
    '''
    if fs is None: fs = glob('lax/*cpm*.fits') # BAD glob, includes the unmosaiced images (c# before .fits) -VP
    if len(fs)==0:
        print "WARNING: No mosaiced cosmic ray pixel masks available."
        return
    
    # Loop over the mosaiced rectified science images
    scifs,scigas = get_ims(glob('rec/*rec*.fits'),'sci') # THIS LOOP IS BAD, includes unmosaiced imgs, maybe add len(name) check - VP
    for i,f in enumerate(scifs):
        ga,imgnum = scigas[i],f[11:f.index('.fits')]
        # Since CPM.shape == f.shape, just change each cosmic ray pixel in the data (f) to -100000 (way below apall threshold)
        cpmfile = 'lax/sci%0.2fcpm%03i.fits'%(ga,imgnum)
        hducpm = pyfits.open(cpmfile)
        datacpm = hducpm[0].data.copy()
        hducpm.close()
        hdu = pyfits.open(f,mode='update')
        data = hdu[1].data.copy()
        cosmicIndices = np.where(datacpm == 1)[0] # 1 ==> cosmic ray pixel
        data[cosmicIndices] = -1000000 # way below apall threshold, so will be rejected during extraction
        hdu[1].data = data 
        hdu.flush() # Just doing the flagging in the original image, not saving to new one
        hdu.close()

def run_extract(fs=None):
    '''
    We will only run apall with one set of default parameters so the user doesn't have to deal with a prompt.
    These default parameters were found by us to be an optimal generalized set.
    The user will be reminded in the terminal that they can change parameters (e.g., nsum) in the interactive window.
    It will be assumed that bright SNe can be extracted just as well using faint SNe apall parameters, but not vice versa.
    Since we are running apall on the non-2D-bkg-subtracted images, local bkg subtraction is necessary so it will be on by default.
    It will be up to the user to carefully define the background regions, especially for faint or extended objects.
    A pre-sky-subtracted-level (pssl) is not added back into the non-2D-bkg-sub image.
    However, it might be necessary based on previous reductions that 2d-bkg-sub images (+pssl) should be used for very faint spectra.
    '''
    if fs is None: fs = glob('rec/*rec*.fits') # BAD glob, includes the unmosaiced images (c# before .fits) -VP
    if len(fs)==0:
        print "WARNING: No rectified images available for extraction."
        return
    
    # For each science image, run apall
    scifs,scigas = get_ims(fs,'sci') # THIS LOOP IS BAD, includes unmosaiced imgs, maybe add len(name) check -VP
    if not os.path.exists('ext'):os.mkdir('ext')
    
    print "APALL: you can change parameters in the interactive apall window."
    print "       No continuum? Make nsum small (~-5) centered on an emission line."
    print "       See pipeline documentation for more help and apall tricks."
    for i,f in enumerate(scifs):
        ga,imgnum = scigas[i],f[11:f.index('.fits')]
        outfile = 'ext/sci%0.2fext%03i.fits'%(ga,imgnum)
        ##### IMPORTANT: make sure relevant keys in 0- or 1-ext-header to prevent pyfits crashes: GAIN (1?), RDNOISE?, DISPAXIS -VP
        ##### IMPORTANT: check if you can change 'line' inside of apall; that is important for non-continuum/extended/faint objects -VP
        #####            this is also where and why ds9 for extract() helps, just nice to see 2D. lots to remember/write otherwise.
        #####            plus, apall window is popping up anyway, unlike in rectify(), so clicks are inevitable
        ##### IMPORTANT:  make sure lsigma and usigma thresholds are enough to reject the cosmic ray pixels (value=1000000) from sum
        saltgain = pyfits.getval(f,'GAIN') 
        # If not a std star reduction, open rectified image in ds9 in case user has doubts
        ''' (apall window is popping up anyway, unlike in rectify, so there will be clicks regardless
            and it's better to be on the safe/convenient side than to have the user open another terminal to start ds9...)
            Could eventually make this optional (for the more interactive version) but since apall is opening anyway,
            and since the user has not yet seen any pipeline-processed data, this might be good before we jump into the 1D realm... '''
        # possiblestds is a global variable (array) of pysalt standard star names defined near the beginning of this code
        if (pyfits.getval(f,'OBJECT')).lower() not in possiblestds: tods9(f)
        iraf.unlearn(iraf.apall) # At long last, prepare to run apall
        iraf.apall(input=f+'[1]',output='auxext.fits',interactive='yes',review='no',line='INDEF',nsum=-1000,lower=-3.5,upper=3.5,
                   b_function='legendre',b_order=1,b_sample='-15:-10,10:15',b_naverage=-5,b_niterate=0,b_low_reject=3.0,
                   b_high_reject=3.0,b_grow=0.0,nfind=1,t_nsum=50,t_step=15,t_nlost=100,t_function='legendre',t_order=3,
                   t_naverage=1,t_niterate=3,t_low_reject=1.0,t_high_reject=1.0,t_grow=1.0,background='fit',weights='variance',
                   pfit='fit1d',clean='yes',readnoise=0,gain=saltgain,lsigma=2.0,usigma=2.0)
        hduaux = pyfits.open('auxext.fits')
        hdraux = hduaux[0].header.copy() # Need to replace old 2D 1-ext-header with this 1D one + transfer keywords as necessary
        hdraux['EXTNAME'] = 'SCI'
        hdraux['EXTVER'] = 1
        dataaux = hduaux[0].data.copy()
        hduaux.close()
        hdu = pyfits.open(f)
        hdu[1].header = hdraux
        hdu[1].data = dataaux
        hdu.writeto(outfile)
        hdu.close()
        os.remove('auxext.fits')
        
        ########### Maybe the apsum for arcs can be added into here somehow to shorten this entire function
        
    #run apsum on the corresponding arc image
    arcfs,arcgas = get_ims(fs,'arc')
    for i,f in enumerate(arcfs):
        ga,imgnum = arcgas[i],f[11:f.index('.fits')]
        outfile = 'ext/arc%0.2fext%03i.fits'%(ga,imgnum)
        pyfits.setval(f,'DISPAXIS',ext=1,value=1) # just in case
        ''' PROBLEM: can the search for refsci be made more efficient? (without choosing unmosaiced chip images) '''
        refsci = ''
        for s in scifs:
            # 19 is standard length for filenames if they don't have c# before .fits
            # and no imcombine was run (which would have changed the imgnum to imgnum1+imgnum2+...)
            if s[3:8] == '%0.2f'%(ga) and len(s) == 19:
                refsci = s # apsum will adopt trace of this corresponding science spectrum
        if refsci == '':
            print "WARNING: no reference science to run apsum on arcs with."
            return
        iraf.unlearn(iraf.apsum)
        iraf.apsum(input=f+'[1]',output='auxext_arc.fits',references=refsci,interactive='no',review='no',background='no'
                   nsum=50,lsigma=2.0,usigma=2.0)
        hduaux = pyfits.open('auxext_arc.fits')
        hdraux = hduaux[0].header.copy() # Need to replace old 1-ext-header with this, but transfer keywords as necessary
        hdraux['EXTNAME'] = 'SCI'
        hdraux['EXTVER'] = 1
        dataaux = hduaux[0].data.copy()
        hduaux.close()
        hdu = pyfits.open(f)
        hdu[1].header = hdraux
        hdu[1].data = dataaux
        hdu.writeto(outfile)
        hdu.close()
        os.remove('auxext_arc.fits')
        
        ##### IMPORTANT: can save lines by verifying that preservefits() does what we want
        ##### and using that to accomplish the aux.fits -> outfile.fits transfer (in other defs too)
            
    
    #run apsum on the corresponding sky image.
    ##### Not doing this right now since generally local bkg sub on non-2d-bkg-sub images yields a decent sky spectrum band.
    
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%% I will return to this tomorrow
def run_checksky(fs=None):
    '''
    Since we do pysalt specidentify manually, there's no reason to do identify1d+dispcor.
    Instead, do a quick cross-correlation of the skylines apall band against a good reference sky spectrum.
    Print a warning if the error (wavelength shift or 'lag') is more than 1AA.
    '''
    
    ##### Where do I get a sky line spectrum that includes the bright bluer lines like [O I] 5577AA, 6300AA, 6363AA?
    ##### mastersky.fits used currently is not good enough since its wavelength coverage begins at ~7250AA (to ~10000AA) -VP
    
    ##### Going with cross-correlation to compute the wavelength shift since it's faster than some kind of fitting.
    ##### Should the shift be applied to the dispersion parameters? I guess only if not > 1AA?
    
    ##### Here's my general code to cross-correlate and find the shift.
    ##### IMPORTANT: figure out what to use in place of mastersky.fits to cross-correlate blue bright lines.
    
     ''' on hold since i wanted to finish up a few of the subsequent tasks -VP '''
    
    if fs is None: fs = glob('ext/sci*ext*.fits')
    if len(fs)==0:
        print "WARNING: No extracted spectra to check sky line wavelengths."
        return

    scifs,scigas = get_ims(fs,'sci') # grabs the spectra and gr-angles
    for i,f in enumerate(scifs):
        # For each extracted science spectrum, grab the data in the 3rd band (by apall's defn: sky band)
        hdu = pyfits.open(f)
        datasky = (hdu[1].data.copy())[2][0]
        hdrsky = hdu[1].header.copy() # for dispersion parameters
        hdu.close()
        # Construct the wavelength array using the dispersion parameters in the extension 1 header
        crpixsky = hdrsky['CRPIX1']
	    cdsky = hdrsky['CD1_1']
	    crvalsky = hdrsky['CRVAL1']
	    wavesky = np.linspace(crvalsky,crvalsky+cdsky*(len(datasky)-crpixsky),len(datasky))
	    
	    # Similarly for skyLineSpectrum, which is a global variable for the path to the reference sky spectrum
	    hdu = pyfits.open(skyLineSpectrum)
	    
        # Need to resample skyLineSpectrum to the same # of data points as the science spectrum (for proper x-corr)
    
    return
    
    
def run_split1d(fs=None)
    '''
    We want to split the spectra into individual chip spectra.
    This will be done using the chipGapPix global dictionary (same as in run_unmosaic()) giving the 
    begin:end pixel numbers of the 2 chip gaps based on the binning ('CCDSUM' keyword).
    '''
    if fs is None: fs = glob('ext/sci*ext*.fits')
    if len(fs)==0:
        print "WARNING: No extracted spectra to split up."
        return

    # scifs,scigas = get_ims(fs,'sci') # not necessary; this way, can use this for stds too (which sometimes have blank 'OBSTYPE'=>crash)
    for i,f in enumerate(fs):
        # Global variable chipGapPix returns the begin:end pixel numbers for the 2 chip gaps (depends on binning, i.e., CCDSUM).
        # See description accompanying chipGapPix definition near the beginning of this script for more information.
        (c1min,c1max),(c2min,c2max) = chipGapPix[pyfits.getval(f,'CCDSUM')][0],chipGapPix[pyfits.getval(f,'CCDSUM')][1]
        hdu = pyfits.open(f)
        data = hdu[1].data.copy()
        # Henceforth, the headers of the split spectra shall consist of both hdr0 and hdr1 keywords.
        hdrcomb = hdu[0].header.copy()+hdu[1].header.copy() # No repetitive keys; left 'EXTNAME'='SCI' and 'EXTVER'=1 for now
        # Do the split such that the middle spectrum has both chip gaps, and the other two spectra have no chip gaps
        ##### We talked about removing the chip gaps entirely in this stage, but the chip gaps themselves
        ##### correspond to wavelengths, so we can't just remove them from the data array or else other data will be
        ##### mapped to the chip gap wavelengths.
        ##### What we could do in this step is to replace the chip gaps with interpolated data which is what we do 
        ##### in the working pipeline. But that method interpolates *over* the chip gap so it uses data in two different
        ##### chips to find the interpolated chip gap data. Not sure if we can do that since the whole point is that 
        ##### we are concerned about the data on different chips being offset, etc.
        ##### Basically: What's the best way to mask the chip gap data if both chip gaps are in the middle spectrum?
        
        # Begin:End pixel numbers for the three different chips
        # Remember: both chip gaps will be in the middle chip spectrum (there will be no chips in the 1st and 3rd chip spectra)
        chip = {'1':(0,c1min),'2':(c1min,c2max+1),'3':(c2max+1,len(data[0][0]))} # '3' end px: all bands have the same array length
        for c in range(1,4):
            data[0][0] = data[0][0][chip[str(c)][0]:chip[str(c)][1]] # band 1 (optimal extraction)
            data[1][0] = data[1][0][chip[str(c)][0]:chip[str(c)][1]] # band 2 (non-optimal extraction)
            data[2][0] = data[2][0][chip[str(c)][0]:chip[str(c)][1]] # band 3 (sky)
            data[3][0] = data[3][0][chip[str(c)][0]:chip[str(c)][1]] # band 4 (sigma)
            tofits(f[:-5]+'c%i.fits'%(c),data,hdr=hdrcomb) # save to same 'ext/' directory, just adding c# before '.fits'
        
        ##### CHECK: will the same dispersion parameters yield the correct wavelengths for the shorter array slices (chips) above?
        ##### Probably need to change CRVAL because the 2nd and 3rd chip data array pixels start at 0 => all chips, similar wavelengths


def run_stdsensfunc(fs=None):
    ''' 
    Produces 1D extracted std star spectrum and corresponding sensfunc. (So, runs onedspec.standard and onedspec.sensfunc
    and moves the non-flux-cal std star spectra and sensfuncs to the pipeline's standards directory (global variable pipeStandardsPath).
    By construction, this function will only run if it finds standards in the current directory. It is up to the user
    to make sure the pipeline is run in a *separate* std star reduction directory where the standards have their own arcs.
    For fast reductions, skip to run_prepstandards() which should pull the appropriate std star spectrum and sensfunc from the archive.
    '''
    if fs is None: fs = glob('ext/sci*ext*c*.fits')
    if len(fs)==0:
        print "WARNING: No extracted spectra to create sensfuncs from."
        return
    # Only run this function for std star reductions. Since it's assumed that the pipeline is run EITHER for science or 
    # std star data, but not a combination of both, it should suffice to check if just one object doesn't have a std star name.
    for f in fs: # hidden but reasonable assumption: if std star reduction, all objects have the same std star name
        if (pyfits.getval(f,'OBJECT')).lower() not in possiblestds: # possiblestds is a global variable of pysalt std star names
            print f+" is not a standard star object; skipping run_stdsensfunc()."
            return
    
    # Run iraf.onedspec.standard and iraf.onedspec.sensfunc for each split chip spectrum
    if not os.path.exists('flx'):os.mkdir('flx')
    stdfs,stdgas = get_ims(fs,'sci') # grabs the spectra and gr-angles (std == sci for std star reductions)
    for i,f in enumerate(stdfs):
        ga,binning,chipnum = '%0.2f'%(stdgas[i]),(pyfits.getval(f,'CCDSUM')).replace(' ',''),f[-6] # chipnum from naming convention
        outfile_std = 'flx/std'+ga+'c'+chipnum
        iraf.unlearn(iraf.onedspec.standard)
        iraf.onedspec.standard(input=f,output=outfile_std,caldir=standardsPath,extinction='',interact='NO',
                               airmass=pyfits.getval(f,'AIRMASS'),exptime=pyfits.getval(f,'EXPTIME'),
                               answer='NO',bandwidth=50.0,bandsep=20.0,star_name='m'+(pyfits.getval(f,'OBJECT')).lower())
        # Need to remove bandpasses (indexed by wavelengths) in the std file which fall in the telluric B or A band
        f = open(outfile_std,mode='r')
        lines = f.readlines()
        f.close()
        hdrline,datalines = lines[0],lines[1:] # don't want to include header line in data array
        f = open(outfile_std,mode='w') # deletes all previous content, will rewrite only non-telluric bandpass data lines
        f.write(hdrline)
        for l in datalines:
            w = float((l.split())[0]) # wavelength (first column)
            if (w>=7600.0 and w<=7630.0) or (w>=6860.0 and w<=6890.0) or (w>=7170.0 and w<=7350.0):
		        continue # don't write line because this is a bad (telluric) bandpass
	        else:
		        f.write(l) # good (non-telluric) bandpass
        f.close()
        # By construction, chip gap 1 is at the beginning of and chip gap 2 is at the end of '*c2.fits' (the middle ["chip 2"] spectrum)
        # So, if this file is *c2.fits, then we need to remove the bandpasses that fall in chip gaps.
        if chipnum == '2':
            # Need to compute the wavelengths of the chip gaps since bandpasses in std file are indexed by wavelength.
            hdu = pyfits.open(f); data = hdu[0].data.copy()[0][0]; hdu.close()
            crval,crpix,cd1 = pyfits.getval(f,'CRVAL1'),pyfits.getval(f,'CRPIX1'),pyfits.getval(f,'CD1_1')
            waves,pixels = np.linspace(crval,crval+cd1*(len(data)-crpix),len(data)),np.arange(len(data))
            # chipGapPix is a global variable (see beginning of code) which we can use to get the length of each chip gap
            # We can't use these pixel numbers since they are for the non-split (all-3-chip) spectra
            (c1min,c1max),(c2min,c2max) = chipGapPix[pyfits.getval(f,'CCDSUM')][0],chipGapPix[pyfits.getval(f,'CCDSUM')][1]
            len_c1,len_c2 = c1max-c1min+1,c2max-c2min+1 # tells us how many pixels the chip gaps span
            # Since the chip gaps are at the beginning and end, use their lengths to find their pix # in this c2*.fits spectrum slice
            c1_maxpix,c2_minpix = len_c1,len(data)-len_c2 # c1_minpix = 0, c2_maxpix = len(data) by construction
            wavesc1 = waves[np.logical_and(pixels>=0,pixels<=c1_maxpix)] # wavelengths of first chip gap
            wavesc2 = waves[np.logical_and(pixels>=c2_minpix,pixels<=len(data))] # wavelengths of second chip gap
            # Finally, remove chip gap bandpasses in the same way as was done for telluric bandpasses
            f = open(outfile_std,mode='r')
            lines = f.readlines()
            f.close()
            hdrline,datalines = lines[0],lines[1:] # don't want to include header line in data array
            f = open(outfile_std,mode='w') # deletes all previous content, will rewrite only non-chip-gap bandpass data lines
            f.write(hdrline)
            for l in datalines:
                w = float((l.split())[0]) # wavelength (first column)
                if (w>=np.min(wavesc1) and w<=np.max(wavesc1)) or (w>=np.min(wavesc2) and w<=np.max(wavesc2)):
		            continue # don't write line because this is a bad (chip gap) bandpass
	            else:
		            f.write(l) # good (non-chip-gap) bandpass
            f.close()
         
        # Now, run iraf.onedspec.sensfunc using the automatically created std bandpass file   
        outfile_sens = 'flx/sens'+ga+'c'+chipnum # Note sensfunc only wants root name, not .fits suffix
        iraf.unlearn(iraf.onedspec.sensfunc)
        iraf.onedspec.sensfunc(standards=outfile_std,sensitivity=outfile_sens,apertures='1',function='spline3',
                               order=2,extinction='',ignoreaps='yes',interactive='NO',answer='NO')
        # Transfer some relevant keywords for user convenience & record-keeping purposes to the sens FITS file
        for k in ['GR-ANGLE','CCDSUM']: pyfits.setval(outfile_sens,k,value=pyfits.getval(f,k))
        pyfits.setval(outfile_sens,'RUSTD',value=f)
        
    ''' Should the user manually move the single-chip sensfuncs and std*flx*.fits files to the pipeline standards archive?
    That's a lot of renaming to do given the naming convention defined in run_prepstandards(). 
    Should the pipeline automatically move the files when the user hasn't even manually inspected the combined std star spectrum? '''


def run_prepstandards(fs=None):
    '''
    For quick, automated reductions using standards with reliable wavelength solutions
    (i.e., standards which came with their own arcs), the pipeline should pull the relevant 1D
    standard star spectra and corresponding sensfuncs from the main pipeline standards directory 
    (defined by the global variable pipeStandardsPath).
    The four criteria for relevant standards are: (1) gr-angle (2) binning (3) most recent (4) chip number
    The naming convention for the std star spectra and corresponding sensfuncs was defined to 
    make it easy to find standards that satisy those criteria. Here's an example:
    '15.12b24std20140604c1.fits' and '15.12b24sen20140604c1.fits',
    where 'b24' is for binning 2x4 ('CCDSUM'). Note that both strings have the same length.
    '''
    if fs is None: fs = glob('ext/sci*ext*c*.fits') 
    if len(fs)==0:
        print "WARNING: No extracted & split spectra to prepare standards for."
        return
    
    # This function should not be run for std star reductions; if std star detected, return
    for f in fs:
        if (pyfits.getval(f,'OBJECT')).lower() in possiblestds: # possiblestds is a global variable of pysalt std star names
            print "Standard star detected ("+f+"); skipping run_prepstandards()."
            return
    
    # Loop over the chip-split extracted science spectra to find corresponding std star spectra and sensfuncs
    scifs,scigas = get_ims(fs,'sci') # grabs the spectra and gr-angles
    if not os.path.exists('flx'):os.mkdir('flx')
    beststds,bestsens = [],[]
    for i,f in enumerate(scifs):
        # First, find standards with the same gr-angle and binning
        binning = (pyfits.getval(f,'CCDSUM')).replace(' ','') # CCDSUM format is '2 4' with the extra space in middle
        ga = '%0.2f'%(scigas[i])
        prefix = ga+'b'+binning # prefix for archived std star filenames based on naming convention
        stds = glob(pipeStandardsPath+prefix+'std*.fits')
        if len(stds) == 0:
            print "WARNING: did not find a standard star to flux calibrate "+f
            return # if you can't fluxcal one spectrum, no point in doing it for other spectra (let user re-run for individual redux)
        # Now, find the most recent standard
        dates = [] 
        for s in stds:
            dates.append(int(s[11:19])) # According to the convenient naming convention
        ''' the most recent date will be the greatest integer if format is yyyymmdd, but since there will be 3 separate chip spectra
        at this date, need to choose the one with the same chip number as the current sci file f '''
        beststd = stds[np.where(dates==np.max(dates))[0][0]]
        # We just chose any one of the 3 separate chip spectra; now simply replace the chipnum with the one of sci file f
        beststd[-6] = f[-6] # 6th character from the end is the chip number
        # Now pick the sensfunc which has the same naming convention except 'std' becomes 'sen' in the filename
        bestsen = beststd[-24:-16]+'sen'+beststd[-13:]
        shutil.copyfile(beststd,'flx/'); shutil.copyfile(bestsen,'flx/')
        beststds.append(beststd); bestsens.append(bestsen)
        # Good idea to add the std and sensfunc filenames to science spectrum header for record-keeping and easy access in fluxcal()
        pyfits.setval(f,'RUSTD',ext=0,value=beststd[-24:]) # beststd contains path => reverse indices
        pyfits.setval(f,'RUSENS',ext=0,value=bestsen[-24:])
        
        ''' POSSIBLE ISSUE
        2014-06-06: I changed this to pull chip-split flux-calibrated std star spectra. Before, I had coded it to pull 
        the non-split extracted std star spectra and sensfuncs, and then call run_split1d() on them. The reason I liked that 
        is because if we ever change the chipGapPix *convention*, there could be a mismatch between the chip gap begin:end 
        pixel numbers in the current (new) science reduction and the archived already-chip-split spectra. It would be ideal 
        to pull the *non-chip-split* flux-calibrated std star spectra so run_split1d() can split them using chipGapPix 
        (which may have been updated since the std star spectra were archived), in the same manner as the science spectra.
        However, this pipeline doesn't put together the chip spectra to re-create the gr-angle spectra. It only combines them
        at the very end to create the full, combined spectrum. Again, this will only be an issue if the carefully chosen conventions 
        defined by me in chipGapPix are changed by someone in the future (could happen if SALT starts using new chips). 
        One solution is to turn the global variable chipGapPix into a function which runs as an early 2D stage, or if it is 
        defined as part of the run_createbpm() task. That task creates the BPM so it finds the bad pixel masks (where counts == 0).
        chipGapPix can then be created dynamically using the chip gap boundaries found by run_createbpm(). But this also requires
        that we figure out a reliable way to merge the chip-split std star spectra into the individual gr-angle spectra 
        (again, a step that is not normally done for science spectra due to the flux scaling + speccombine step). -VP
        '''


def run_fluxcal(fs=None):
    ''' We have the sensfunc for each chip for each science and std star spectrum. So just run calibrate. '''
    if fs is None: fs = glob('ext/sci*ext*c*.fits') # grabs individual chip science spectra
    if len(fs)==0:
        print "WARNING: No science chip spectra to flux calibrate."
        return
    
    # Loop over the science chip spectra, pull the relevant std star and sensfunc chip spectra as needed
    scifs,scigas = get_ims(fs,'sci') # grabs the spectra and gr-angles
    if not os.path.exists('flx'):os.mkdir('flx') # in case for some reason run_prepstandards() wasn't run
    for i,f in enumerate(scifs):
        ga,binning,chipnum = '%0.2f'%(scigas[i]),(pyfits.getval(f,'CCDSUM')).replace(' ',''),f[-6] # chipnum from naming convention
        outfile = f[-24:-13]+'flx'+f[-10:] # based on naming convention, just replace 3-char stage 'ext' with 'flx'
        # Grab the std star spectrum and sensfunc by assuming that only 1 of each is present for each gr-angle, see run_prepstandards()
        std,sens = glob('flx/'+ga+'b'+binning+'std*c'+chipnum+'.fits')[0],glob('flx/'+ga+'b'+binning+'sen*c'+chipnum+'.fits')[0]
        outstd = std[-24:-16]+'flx'+std[-13:] # See run_stdsensfunc() for naming convention; replacing 'std' with 'flx' preserves len(std)
        # Finally, calibrate the science and standard star spectra (note: sensitivity=sens without .fits)
        iraf.unlearn(iraf.calibrate)
        iraf.calibrate(input=f,output='auxflx_sci.fits',extinct='no',extinction='',sensitivity=sens[-24:-5],ignoreaps='yes',
                       exptime=pyfits.getval(f,'EXPTIME'),airmass=pyfits.getval(f,'AIRMASS'))
        iraf.unlearn(iraf.calibrate)
        iraf.calibrate(input=std,output='auxflx_std.fits',extinct='no',extinction='',sensitivity=sens[-24:-5],ignoreaps='yes',
                       exptime=pyfits.getval(std,'EXPTIME'),airmass=pyfits.getval(std,'AIRMASS'))
        # Transfer the data to a copy of the existing chip spectra (with the correct headers)
        transfers = [('auxflx_sci.fits',f,outfile),('auxflx_std.fits',std,outstd)]
        for aux,ref,new in transfers:
            hduaux = pyfits.open(aux)
            dataaux = hduaux[0].data.copy()
            hduaux.close()
            hdu = pyfits.open(ref)
            hdu[0].data = dataaux
            hdu.writeto(new)
            hdu.close()
            # Transfer calibrate's new keywords to the existing header: extinction, calibration, dispersion sampling flags, and flux units
            for k in ['EX-FLAG','CA-FLAG','DC-FLAG','BUNIT']: pyfits.setval(new,k,value=pyfits.getval(aux,k))
            os.remove(aux)


# %%%%%%%%%%%%%%%%%%%%%%%%%%%%% I will clean this up tomorrow or friday, after talking to Curtis
def run_pytelluric(fs=None):
    '''
    Essentially copied from our working pipeline's relatively standalone PyTelluric module.
    Go over this with Curtis piece by piece in great detail. Figure out how to incorporate distortion
    or improve the stretching of the std star telluric regions to fit the shape of the science telluric regions better.
    I'd like for us to turn this into a standalone code (like lacosmicx).
    This currently works some of the time (spikey corrections otherwise), and best if used with stds with their own arcs.
    Note: consider adding a step to correct for redshift due to Earth's velocity relative to object as in final.pro IDL routine.
    '''
    if fs is None: fs = glob('flx/*flx*c*.fits') # grabs individual chip science spectra
    if len(fs)==0:
        print "WARNING: No flux-calibrated science chip spectra to telluric-correct."
        return
    
    # Loop over the science chip spectra, pull the relevant std star and sensfunc chip spectra as needed
    scifs,scigas = get_ims(fs,'sci') # grabs the spectra and gr-angles
    if not os.path.exists('tel'):os.mkdir('tel')
    for i,f in enumerate(scifs):
        ga,binning,chipnum = '%0.2f'%(scigas[i]),(pyfits.getval(f,'CCDSUM')).replace(' ',''),f[-6] # chipnum from naming convention
        outfile = f[-24:-13]+'tel'+f[-10:] # based on naming convention, just replace 3-char stage 'flx' with 'tel'
        # Grab the std star spectrum by assuming that only 1 is present for each gr-angle&chip, see run_prepstandards()
        std = glob('flx/'+ga+'b'+binning+'flx*c'+chipnum+'.fits')[0]
        
        ''' begin copy paste, with very minor modifications (names, etc.) ... go through it in detail later '''
        hdusci = pyfits.open(f)
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

        telluricLines = [(beginA,endA,beginAinterp,endAinterp),(beginB,endB,beginBinterp,endBinterp)] # not doing a-band (called by me "C")

        # Begin telluric correction process.
        corsci = datsci.copy()
        for beginT,endT,beginTinterp,endTinterp in telluricLines:
            wavestdTOrig = wavestd[np.where(np.logical_and(wavestd>beginT,wavestd<endT))] # For interactive comparison purposes.
            datstdTOrig = datstd[np.where(np.logical_and(wavestd>beginT,wavestd<endT))] # For interactive comparison purposes.
            #print("Working in telluric region... "+str(beginT)+":"+str(endT)+" Angstroms.")
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
            #print("Standard star shifted by "+str(shiftT)+ " Angstroms.")
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
                #print("Correcting for len(datSTD) > len(datasci) by "+str(np.abs(len(datstdT)-len(datsciT)))+"px")
                for i in np.arange(np.abs(len(datstdT)-len(datsciT)))+1: # remove last element
	                wavestdT = wavestdT[0:-i]
	                datstdT = datstdT[0:-i]
            elif len(datstdT) < len(datsciT):
                #print("Correcting for len(datSTD) < len(datasci) by "+str(np.abs(len(datstdT)-len(datsciT)))+"px")
                for i in np.arange(np.abs(len(datstdT)-len(datsciT)))+1: # add element from shifted std star arrays
	                maxInd = (np.where(np.logical_and(wavestdS>beginT,wavestdS<endT)))[0].max()
	                wavestdT = np.append(wavestdT,wavestdS[maxInd+i]) # extend the smaller array to match the larger
	                datstdT = np.append(datstdT,datstd[maxInd+i])
            # Scale, shift, and stretch standard star sub-spectrum to match science sub-spectrum.
            fitfunc = lambda p,x: (p[0]+p[1]*x)*datstdT
            errfunc = lambda p,x,y: fitfunc(p,x) - y
            p0 = [1.0,0.0]
            p1,success = leastsq(errfunc,p0,args=(wavesciT,datsciT))
            #print("Standard star polynomial scaling parameters:")
            #print(p1)
            #datstdT = (p1[0] + p1[1]*wavestdT)*datstdT
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
            #print("Correction polynomial scaling parameters:")
            #print p1
            corsclT = (p1[0]+p1[1]*wavesciT+p1[2]*wavesciT*wavesciT)*corT
            # Replace telluric band in entire spectrum with the corrected and scaled data values.
            corsci[np.where(np.logical_and(wavesci>beginT,wavesci<endT))] = corsclT.astype(corsci.dtype)
        # De-regularize entire corrected science spectrum.
        sciFinal = corsci*np.std(datsciOrig)+np.mean(datsciOrig)
        # Copy telluric-corrected data to new file.
        hdu = pyfits.open(f)
        hdu[0].data = sciFinal
        hdu.writeto(outfile)
        hdu.close()
    
    
def run_sigmaclip():
    return
    
def run_fluxscale()
    return

def run_speccombine():
    return

def run_snapshots():
    ''' 
    Create an images folder and save snapshots of images/spectra after each reduction step.
    These will be used by the Rutgers Supernova Database program when creating the 
    webpage for each SN. And it's a quicker check on the outputs than manually opening each file.    
    This should not be run as part of the reduction process since it can take a while.
    '''
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
    
