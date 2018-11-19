TUWmodel
--------
Python implementation of the rainfall-runoff TUWmodel developed by Parajka et al. (2007). The model is based on the original Fortran code, updated to Fotran90 and optimized for use in the Python library. This python package is still work in progress, but can be used already (2018-11-19). 

Installation
------------
This Fortran code has been compiled for MacOS 10.14. A windows compiled file will soon be added. The python package f2py is used for compilation of the Fortran90 code. 

For installation of the python-package run::

  pip install https://github.com/raoulcollenteur/tuwmodel/zipball/master

References
----------
Parajka, J., R. Merz, G. Bloeschl (2007) Uncertainty and multiple objective calibration in regional water balance modelling: case study in 320 Austrian catchments, Hydrological Processes, 21, 435-446.
