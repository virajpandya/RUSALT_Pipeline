RUSALT_Pipeline
===============

modularized RU SALT pipeline


initial upload: 2014-03-16

- dicts.py: initializes global dictionaries used by other modules to refer to images
- driver.py: run this to begin the pipeline
- fullReduce: runs tasks in tasks.py in the correct order for a full reduction
- indivReduce: lets one run individual tasks
- params.py: customizable parameters like linelist / standards path
- pipeHistory: functions that add/update keys in FITS headers after each task
- pyrafCalls: second-to-main function that sets parameters for and runs PyRAF tasks
- sortDicts: used to sort images in directory into global dictionaries (dicts.py)
- tasks.py: CORE of the pipeline, one function per reduction step, calls functions from pyrafCalls.py to actually run PyRAF tasks. The order of the functions mimic the order/steps of an actual full reduction as seen in fullReduce.py
