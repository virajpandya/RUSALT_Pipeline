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
iraf.apextract()
iraf.imutil()

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
            ims.append(f)
            grangles.append(pyfits.getval(f,'GR-ANGLE'))
    return array(ims),array(grangles)

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
        print "ERROR: No mosaiced arcs for PySALT's (2D) specidentify."
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
        except KeyError:
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
        print "ERROR: No wavelength solutions for rectification."
        return
    # rectify each sci/arc and pop it open briefly in ds9
    scifs,scigas = get_ims(glob('mos/*mos*.fits'),'sci')
    arcfs,arcgas =get_ims(glob('mos/*mos*.fits'),'arc')
    
    ims = append(scifs,arcfs)
    gas = append(scigas, arcgas)
    
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
        print "ERROR: No rectified images to split by chip."
        return
    # Grab the rectified science images (only they need to be split up for bkg+lax)
    scifs,scigas = get_ims(fs,'sci')
    '''
    There are two chip gaps, each spanning a range of columns (pixels). To avoid using data near the chip gap edges (primarily
    due to skylines and the fact that we have data covering the chip gaps from overlapping spectra), we take the edges of the chip
    gaps to be +/- ~20px beyond their definite edges. The begin:end pixel numbers of each chip gap will change if the binning
    is different ('CCDSUM' keyword), so the dict chipGapPix below allows one to access the pixel numbers based on 'CCDSUM'.
    
    We're splitting by the 3 chips now instead of the 6 'amps' as in the pre-saltmosaic multi-FITS extension files before for two 
    reasons. First, where the chip gaps begin and end is visible in the 2D spectra (counts == 0), but not for the amps since they are just
    based on the pre-saltmosaic images, e.g., the second amp begins directly after the first amp. More importantly, the rectification
    probably smeared out where one amp ends and the next begins (just look at how it curves the chip gaps).
    '''
    chipGapPix = {'2 2':((1010,1095),(2085,2160)),'2 4':((1035,1121),(2115,2191)),'4 4':((500,551),(1035,1091))}
    # Split each image into 3 imgs such that the middle img has both chip gaps, and the other 2 imgs don't have chip gaps
    for i,f in enumerate(scifs):
        ga,imgnum = scigas[i],f[11:f.index('.fits')]
        (c1min,c1max),(c2min,c2max) = chipGapPix[pyfits.getval(f,'CCDSUM')][0],chipGapPix[pyfits.getval(f,'CCDSUM')][1]
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
        print "ERROR: No rectified chip-based images for bad pixel mask creation."
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
        # for easy access in run_lacosmicx() 
        pyfits.setval(f,'BPM',value=outfile) # will be propagated down throughout the stages ;)   
    

def run_background(fs=None):
    if fs is None: fs = glob('rec/*rec*c*.fits') 
    if len(fs)==0:
        print "ERROR: No rectified chip-based images for 2D-background-subtraction."
        return
    
    # Get rectified science images and gr-angles
    scifs,scigas = get_ims(fs,'sci')
    if not os.path.exists('bkg'):os.mkdir('bkg')
    for i,f in enumerate(scifs):
        # Run automated 2D background subtraction on each individual chip image
        pyfits.setval(f,'DISPAXIS',value=1) # just in case since this is auto background
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
    fs = glob('bkg/*bkg*.fits')
    if len(fs)==0:
        print "ERROR: No background-subtracted files for LaCosmicX."
        return
    
    # Get background-subtracted chip-based science images and gr-angles.
    scifs,scigas = get_ims(fs,'sci')
    if not os.path.exists('lax'):os.mkdir('lax')
    for i,f in enumerate(scifs):
        outname_img = 'lax/'+f[4:12]+'lax'+f[15:]
        outname_msk = 'lax/'+f[4:12]+'cpm'+f[15:] # cosmic (ray) pixel mask
        hdu = pyfits.open(f)
        datain = hdu[0].data.copy()
        ##### make sure gain == 1 (since mult='yes' in saltgain), but what about readnoise?
        hdr = hdu[0].header
        saltgain = hdr.get('GAIN',1.0)
        namebpm = hdr['BPM']
        # Get mode of image counts as a proxy for pre-sky-subtracted-level
        mode = float(imutil.imstat(images=f,fields='mode',format='no',Stdout=1)[0])
        # Get BPM data array
        hdubpm = pyfits.open(namebpm)
        databpm = hdubpm[0].data.copy()
        hdubpm.close()
        # Run lacosmicx
        datalax = lacosmicx.run(inmat=datain,inmask=databpm,outmaskfile=outname_msk,
		                        sigclip=6.0,objlim=3.0,sigfrac=0.1,gain=saltgain,pssl=mode,robust=True)
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
    fs = glob('bpm/*bpm*.fits')
    if len(fs)==0:
        print "ERROR: No bad pixel masks to remosaic." # without BPM, wouldn't have made CPM anyway, so this check is sufficient
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
            # Add 'BPM' or 'CPM' to rectified image header
            pyfits.setval(f,masktype.upper(),value=outfile) # will be propagated down throughout the stages ;)  
        

def run_flagcosmics(fs=None):
    '''
    Instead of using the lacosmicx cleaned output image, we will flag cosmic ray pixels in the science data
    such that apall will reject those pixels when doing the variance-weighted extraction.
    '''
    fs = glob('lax/*cpm*.fits') # BAD glob, includes the unmosaiced images (c# before .fits) -VP
    if len(fs)==0:
        print "ERROR: No mosaiced cosmic ray pixel masks available."
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
    fs = glob('rec/*rec*.fits') # BAD glob, includes the unmosaiced images (c# before .fits) -VP
    if len(fs)==0:
        print "ERROR: No rectified images available for extraction."
        return
    
    # For each science image, run apall
    scifs,scigas = get_ims(fs,'sci') # THIS LOOP IS BAD, includes unmosaiced imgs, maybe add len(name) check - VP
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
        #####            plus, apall window is popping up anyway, unlike in rectify()
        ##### IMPORANT:  make sure lsigma and usigma thresholds are enough to reject the cosmic ray pixels (value=1000000) from sum
        saltgain = pyfits.getval(f,'GAIN') 
        iraf.unlearn(iraf.apall)
        iraf.apall(input=f+'[1]',output='auxext.fits',interactive='yes',review='no',line='INDEF',nsum=-1000,lower=-3.5,upper=3.5,
                   b_function='legendre',b_order=1,b_sample='-15:-10,10:15',b_naverage=-5,b_niterate=0,b_low_reject=3.0,
                   b_high_reject=3.0,b_grow=0.0,nfind=1,t_nsum=50,t_step=15,t_nlost=100,t_function='legendre',t_order=3,
                   t_naverage=1,t_niterate=3,t_low_reject=1.0,t_high_reject=1.0,t_grow=1.0,background='fit',weights='variance',
                   pfit='fit1d',clean='yes',readnoise=0,gain=saltgain,lsigma=2.0,usigma=2.0)
        hduaux = pyfits.open('auxext.fits')
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
        os.remove('auxext.fits')
        
    #run apsum on the corresponding arc image
    arcfs,arcgas = get_ims(fs,'arc')
    for i,f in enumerate(arcfs):
        ga,imgnum = arcgas[i],f[11:f.index('.fits')]
        outfile = 'ext/arc%0.2fext%03i.fits'%(ga,imgnum)
        pyfits.setval(f,'DISPAXIS',ext=1,value=1) # just in case
        ''' PROBLEM: can the search for refsci be made more efficient? (without choosing unmosaiced chip images) '''
        refsci = ''
        for s in scifs:
            if s[3:8] == '%0.2f'%(ga) and len(s) == 19: # standard length for filenames if they don't have c# before .fits (and no imcombine)
                refsci = s # apsum will adopt trace of this corresponding science spectrum
        if refsci == '':
            print "ERROR: no reference science to run apsum on arcs with."
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
        
        ##### IMPORTANT: can save lines by verifying preservefits() and using that to accomplish the aux.fits -> outfile.fits transfer
            
    
    #run apsum on the corresponding sky image.
    ##### Not doing this right now since generally local bkg sub on non-2d-bkg-sub images yields a decent sky spectrum band.
    

def run_checksky(fs=None):
    '''
    Since we do pysalt specidentify manually, there's no reason to do identify1d+dispcor.
    Instead, do a quick cross-correlation of the skylines band against a good reference sky spectrum.
    Print a warning if the error (wavelength shift/lag) is more than 1AA.
    '''
    
    ##### What's the best way to estimate the "error"? Cross-correlation shift/lag or some kind of fitting?
    ##### Should the shift be applied to the dispersion parameters? I guess only if not > 1AA?
    
    
    return
    
def run_split1d(fs=None)
    '''
    We want to split the individual spectra into 
    '''

    ##### Turn chipGapPix in run_remosaic() into a global variable
    
    return

def run_stdsensfunc():
    return

def run_fluxcal():
    return 

def run_stdtelluric():
    return 

def run_telluric():
    return 
    
def run_fluxscale()
    return

def run_speccombine():
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
    
