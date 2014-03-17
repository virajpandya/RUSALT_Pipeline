<<<<<<< HEAD
'''
>>> Rutgers SALT Supernova Spectral Reduction Pipeline <<<

This module initializes the global dictionaries which will be 
used across all other pipeline modules to store and sort the 
data filenames which may be accessed by the pipeline. 

A dictionary's elements will consist of FITS filenames while 
its indices will consist of the grating angles corresponding
to this files. The grating angles are taken from the file's header.

Please refer to the documentation for more information
on what each dictionary's elements represents and how the pipeline
accesses, modifies, and works with the global dictionaries.

*** Modifications ***
Sept. 26, 2013: Created module. -Viraj Pandya

'''

flats = {} # Unprocessed 2-D flat-field images.
arcs = {} # Unprocessed 2-D arc lamp images.
standards = {} # Unprocessed 2-D standard star images.
sciences = {} # Unprocessed 2-D science images.
combflats = {} # Combined 2-D flat-field images.
normflats = {} # Normalized 2-D flat-field images.
flatsciences = {} # Flat-fielded 2-D science images.
flatarcs = {} # Flat-fielded 2-D arc images.
flatstandards = {} # Flat-fielded 2-D standard star images.
bpmsciences = {} # Bad pixel masks for 2-D science images.
bpmarcs = {} # Bad pixel masks for 2-D arc images.
bpmstandards = {} # Bad pixel masks for 2-D standard star images.
laxsciences = {} # Lacosmicx-corrected 2-D science images.
bpmlaxsciences = {} # Bad pixel masks for Lacosmicx-corrected 2-D science images.
wavesols = {} # Preliminary wavelength solutions.
wavearcs = {} # Wavelength-corrected (rectified) 2-D arc images.
wavesciences = {} # Wavelength-corrected (rectified) 2-D science images.
wavestandards = {} # Wavelength-corrected (rectified) 2-D standard star images.
backgroundsciences = {} # Background-subtracted 2-D science images.
extractedsciences = {} # Extracted 1-D science spectra.
extractedstandards = {} # Extracted 1-D standard star spectra.
extractedarcs_std = {} # Extracted 1-D arc spectra corresponding to standard star spectra.
extractedarcs_sci = {} # Extracted 1-D arc spectra corresponding to science spectra.
apsumsciences = {} # Extracted 1-D science spectra (non-background-subtracted).
dispsciences = {} # Dispersion-corrected 1-D science spectra (with secondary wavelength solutions).
dispstandards = {} # Dispersion-corrected 1-D standard star spectra (with secondary wavelength solutions).
bwmsciences = {} # Bad wavelength masks for 1-D science spectra.
bwmstandards = {} # Bad wavelength masks for 1-D standard star spectra.
stdfiles = {} # Flux calibration standard star bandpass-flux files.
sensfiles = {} # Flux calibration sensitivity function FITS files.
fluxsciences = {} # Flux-calibrated 1-D science spectra.
fluxstandards = {} # Flux-calibrated 1-D standard star spectra.
tellsciences = {} # Telluric-corrected 1-D science spectra.
tellstandards = {} # Telluric-corrected 1-D standard star spectra.
scaledfluxsciences = {} # Flux-scaled 1-D science spectra.
scaledstandards = {} # Flux-scaled 1-D standard star spectra.
skysciences = {} # Sky spectra corresponding to 1-D science spectra.
sigmasciences = {} # Sigma spectra corresponding to 1-D science spectra.
scaleddispsciences = {} # Count-scaled, dispersion-corrected (non-flux-calibrated) 1-D science spectra.
=======
'''
>>> Rutgers SALT Supernova Spectral Reduction Pipeline <<<

This module initializes the global dictionaries which will be 
used across all other pipeline modules to store and sort the 
data filenames which may be accessed by the pipeline. 

A dictionary's elements will consist of FITS filenames while 
its indices will consist of the grating angles corresponding
to this files. The grating angles are taken from the file's header.

Please refer to the documentation for more information
on what each dictionary's elements represents and how the pipeline
accesses, modifies, and works with the global dictionaries.

*** Modifications ***
Sept. 26, 2013: Created module. -Viraj Pandya

'''

flats = {} # Unprocessed 2-D flat-field images.
arcs = {} # Unprocessed 2-D arc lamp images.
standards = {} # Unprocessed 2-D standard star images.
sciences = {} # Unprocessed 2-D science images.
combflats = {} # Combined 2-D flat-field images.
normflats = {} # Normalized 2-D flat-field images.
flatsciences = {} # Flat-fielded 2-D science images.
flatarcs = {} # Flat-fielded 2-D arc images.
flatstandards = {} # Flat-fielded 2-D standard star images.
bpmsciences = {} # Bad pixel masks for 2-D science images.
bpmarcs = {} # Bad pixel masks for 2-D arc images.
bpmstandards = {} # Bad pixel masks for 2-D standard star images.
laxsciences = {} # Lacosmicx-corrected 2-D science images.
bpmlaxsciences = {} # Bad pixel masks for Lacosmicx-corrected 2-D science images.
wavesols = {} # Preliminary wavelength solutions.
wavearcs = {} # Wavelength-corrected (rectified) 2-D arc images.
wavesciences = {} # Wavelength-corrected (rectified) 2-D science images.
wavestandards = {} # Wavelength-corrected (rectified) 2-D standard star images.
backgroundsciences = {} # Background-subtracted 2-D science images.
extractedsciences = {} # Extracted 1-D science spectra.
extractedstandards = {} # Extracted 1-D standard star spectra.
extractedarcs_std = {} # Extracted 1-D arc spectra corresponding to standard star spectra.
extractedarcs_sci = {} # Extracted 1-D arc spectra corresponding to science spectra.
apsumsciences = {} # Extracted 1-D science spectra (non-background-subtracted).
dispsciences = {} # Dispersion-corrected 1-D science spectra (with secondary wavelength solutions).
dispstandards = {} # Dispersion-corrected 1-D standard star spectra (with secondary wavelength solutions).
bwmsciences = {} # Bad wavelength masks for 1-D science spectra.
bwmstandards = {} # Bad wavelength masks for 1-D standard star spectra.
stdfiles = {} # Flux calibration standard star bandpass-flux files.
sensfiles = {} # Flux calibration sensitivity function FITS files.
fluxsciences = {} # Flux-calibrated 1-D science spectra.
fluxstandards = {} # Flux-calibrated 1-D standard star spectra.
tellsciences = {} # Telluric-corrected 1-D science spectra.
tellstandards = {} # Telluric-corrected 1-D standard star spectra.
scaledfluxsciences = {} # Flux-scaled 1-D science spectra.
scaledstandards = {} # Flux-scaled 1-D standard star spectra.
skysciences = {} # Sky spectra corresponding to 1-D science spectra.
sigmasciences = {} # Sigma spectra corresponding to 1-D science spectra.
scaleddispsciences = {} # Count-scaled, dispersion-corrected (non-flux-calibrated) 1-D science spectra.
>>>>>>> origin/master
combinedspectra = {} # Final combined spectra (indexed by 'dsp', 'flx', etc. instead of grating angles).