
# >>> Rutgers SALT Supernova Spectral Reduction Pipeline <<<
#
# This CUSTOMIZABLE module allows the user to specify parameters
# which will be used by the pipeline or PyRAF tasks. These
# parameters can include, for example, system-specific file paths,
# or specific values for a PyRAF task's parameters, and they will
# replace optimal parameters set automatically by the pipeline.
#
# The user is allowed to edit this file. However, please refer
# to the documentation to preserve the structure of this file 
# while adding to it because the pipeline will not be able to
# parse parameters unless it can follow the file structure. It is
# also not recommended to change PyRAF parameters related to input
# and output file names.
# 
# *** Modifications ***
# Sept. 30, 2013: Created module. -Viraj Pandya
#
# This global variable will prevent repetitive calling of sortDicts.sort()
firstIndivReduce = True

# Line lists for arc line identification:
lineListPath = '/home/viraj/orion/python/pysalt/data/linelists/'
Ar = lineListPath+'Argon_hires.salt'
ThAr = lineListPath+'ThAr.salt'
Xe = lineListPath+'Xe.txt'
Ne = lineListPath+'NeAr.salt' 
CuAr = lineListPath+'CuAr.txt'
HgAr = lineListPath+'HgAr.txt'

# Standard star flux calibration files:
standardsPath = '/home/viraj/orion/python/pysalt/data/standards/spectroscopic/'

# The following are pixels defining the two chip gaps in each spectrum.
# There are two px limits (+/- error): for CCDSUM='2 4' (2x4 binning) or '4 4' (4x4).
# The pipeline will use either one of these based on header keywords.
# The structure is a tuple of 2 tuples: ((chip1min,chip1max),(chip2min,chip2max))
chipGapPix22 = ((1010,1095),(2085,2160)) # added 2014-05-15 for Luke Hovey's SNR - VP
chipGapPix24 = ((1035,1121),(2115,2191))
chipGapPix44 = ((500,551),(1035,1091))

# 'SLOW' and 'FAINT' RSS multiplicative gain = simple average of values for each of the 6 amps
# 2014-07-18: should be 1.00 because the PySALT task saltgain already does the multiplication.
gain = 1.00
gain_dict = {'faint_slow':1.00}




