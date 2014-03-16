'''
>>> Rutgers SALT Supernova Spectral Reduction Pipeline <<<

This module calls functions in the tasks.py module in the order 
needed to do a full reduction of the data in the working directory.

Please refer to the documentation for more information
about the full reduction process.

*** Modifications ***
Sept. 30, 2013: Created module. -Viraj Pandya

'''

import dicts # These are the pipeline's global dictionaries.
import params # Customizable parameters for the pipeline.
import pyrafCalls # These are the pipeline's PyRAF call functions.
import tasks # Core pipeline module.
import sortDicts # Pipeline's module for sorting files into global dictionaries.

def run():
	# sort through all the fits files and create the initial four dictionaries: flats, arcs, standards, sciences
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Sorting through all the images using: dictionaries()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	sortDicts.sort()
	# combine those flats that have the same GR-ANGLE header keyword value
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Combining flats with same GR-ANGLE using: combineflats()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	tasks.combineflats()
	# normalize all the combined flats
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Normalizing each combined flat using: normalizeflats()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	tasks.normalizeflats()
	# flat-field the science images
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Flat-fielding each science image with its associated normalized flat using: flattensciences()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	tasks.flattensciences()
	# flat-field the arc images
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Flat-fielding each arc image with its associated normalized flat using: flattenarcs()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	tasks.flattenarcs()
	# flat-field the standard star images
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Flat-fielding each standard star image with its associated normalized flat using: flattenstandards()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	tasks.flattenstandards()
	# create bad pixel masks for all flat-fielded sciences, arcs, standards
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Creating a bad pixel mask for each flat-fielded science, arc, and standard star image using maskimages()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	tasks.maskimages()
	# run lacosmicx on each flat-fielded science image
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Running LAcosmicx on each flat-fielded science image using LAxsciences()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	tasks.LAxsciences()
	# use pysalt.specidentify to produce wavelength solution files for each flat-fielded arc image
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Identifying the emission lines in each flat-fielded arc image using: specidentifyarcs()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	tasks.specidentifyarcs()
	# use wavelength solution files to wavelength-calibrate the 2D arc images
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Applying wavelength solution to each 2D flat-fielded arc image using: wavecalarc()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	tasks.wavecalarc()
	# use wavelength solution files to wavelength-calibrate the 2D science images
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Applying wavelength solution to each 2D flat-fielded science image using: wavecalsci()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	tasks.wavecalsci()
	# use wavelength solution files to wavelength-calibrate the 2D standard star images
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Applying wavelength solution to each 2D flat-fielded standard star image using: wavecalstd()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	tasks.wavecalstd()
	# use backgroundsciences() to background-subtract the 2D science images
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Subtracting background from each 2D wavelength-calibrated science image using: subtractbackground()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	tasks.subtractbackground()
	# extract the wavelength-corrected science spectra
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Extracted each wavelength-calibrated science spectrum using: extractsciences()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	tasks.extractsciences()
	# extract the wavelength-corrected standard star spectra
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Extracting each wavelength-calibrated standard star spectrum using: extractstandards()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	tasks.extractstandards()
	# extract the wavelength-corrected arc spectra (twice -- for science and for standard star images)
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Extracting each wavelength-calibrated arc spectrum twice (for science and for standard star) using: extractarcs()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	tasks.extractarcs()
	# run apsum on non-background subtracted science images to get "true" sky and sigma spectra
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Extracting the "true" sky and sigma spectra from non-background-subtracted science images using: extractSciSigma()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	tasks.extractSciSigma()
	# run identify on wavelength-corrected, extracted arcs to reduce errors further: just read in linelist, improve fit
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Identifying emission lines in the wavelength-calibrated arc spectra to reduce errors further using: identifyarcs'
	print 'Simply read in the linelist, and improve the fit by deleting outliers'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	tasks.identifyarcs()
	# run reidentify for remaining standard-star-associated wavelength-corrected, extracted arcs
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Automatically identifying the emission lines in the wavelength-corrected arcs associated with the standard stars using: reidentifyarcs'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	tasks.reidentifyarcs()
	# apply lower-error wavelength solution to wavelength-corrected, extracted science images
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Applying new lower-error wavelength solution to each science spectrum using: dispcorsci()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	tasks.dispcorsci()
	# apply lower-error wavelength solution to wavelength-corrected, extracted standard star images
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Applying new lower-error wavelength solution to each standard star spectrum using: dispcorstd()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	tasks.dispcorstd()
	# create bad wavelength masks for dispersion-corrected science (and standard star) spectra
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Creating bad wavelength masks for dispersion-corrected science and standard star spectra using: maskspectra()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	tasks.maskspectra()
	# flux-calibrate each extracted, dispersion-corrected science image
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Flux-calibrating the science and standard star spectra using: fluxcal()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	tasks.fluxcal()
	# Sigma-clip science and standard star spectra
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Sigma-clipping science and standard star spectra using: fluxcal()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	tasks.sigmaClipSpectra()
	# scale flux-calibrated science spectra so that flux values are similar based on overlapping wavelengths
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Scaling flux-calibrated science spectra using: scaleFluxScienceSpectra()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	tasks.scaleFluxScienceSpectra()
	# scale flux-calibrated standard star spectra so that flux values are similar based on overlapping wavelengths
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Scaling flux-calibrated standard star spectra using: scaleStdStarSpectra()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	tasks.scaleStdStarSpectra()
	# extract the sky and sigma spectra from the flux-calibrated science spectra
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Copying the sky and sigma spectra to separate spectral files using: copySciSkySig()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	tasks.copySciSkySig()
	# scale dispersion-corrected science spectra for combinespectra()
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Scaling dispersion-corrected science spectra using: scaleDispScienceSpectra()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	tasks.scaleDispScienceSpectra()
	# combine dispersion-corrected (and flux-calibrated) science spectra
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Combining dispersion-corrected and perhaps (scaled) flux-calibrated science spectra using: combinespectra()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	tasks.combinespectra()
	# clean spectra of leftover cosmic rays and background oversubtraction
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Cleaning spectra of leftover residuals using: cleanCombinedSpectra()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	tasks.cleanCombinedSpectra()
	# telluric-correct spectra
	print '>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	print 'Correcting telluric features in spectra: telluricCorrect()'
	print '<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	tasks.telluricCorrect()
