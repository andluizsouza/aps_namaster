Description

The present pipeline computes the angular power spectrum (APS) of a masked scalar field, that is usually a galaxy catalog with limited survey area, using NaMaster library (https://github.com/LSSTDESC/NaMaster) developed by David Alonso, Javier Sanchez and Anže Slosar for the LSST Dark Energy Science Collaboration.

In general description, NaMaster is a C library, Python module and standalone program to compute full-sky angular cross-power spectra of masked, spin-0 and spin-2 fields with an arbitrary number of known contaminants using a pseudo-Cl (aka MASTER) approach. The code also implements E/B-mode purification and is available in both full-sky and flat-sky modes. Details about this approach  can be found in: https://arxiv.org/abs/1809.09603. 

The present pipeline runs with the Python module of NaMaster, the pymaster (https://namaster.readthedocs.io/en/latest/index.html). It has already been optimized to compute APS of many galaxy catalogs simultaneously using parallelization method of mpi4py. One direct application of measures pipeline is to calculate covariance matrix on a sample mock catalogs for cosmological analysis.
