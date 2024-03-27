# Preporessing-HMI-Magnetograms
Python scripts to pre-process HMI-SHARP magnetograms with all three magnetic field components (Bx, By and Bz) 
when a list of SHARP regions is supplied (e.g. SHARPS.txt).


The Python script prepare_data.py 
downloads three components of HMI SHARP patches
using hmi.sharp_ces_720s data series (http://jsoc.stanford.edu/ajax/lookdata.html?ds=hmi.sharp_cea_720s), 
the chosen time for each patch is when it's closest to the central meridian. Thus, we can approximate 
Bp: B_phi -> Bx, Bt: B_theta -> -By, Br -> Bz


The Python script AR_cutout.py is to be used after downloading the data
This script processes magnetograms to 
1. extract the major structures from a given magnetic field distribution on a Cartesian plane
2. It also calculates the tilt angle created between the positive and negative polarities
 of a given magnetogram, and rotates it such that the tilt will become zero
3. Add extra padding with zero field strength, such that the centroid of the magnetogram 
 coincides with the origin of the new coordinates
4. These preprocessing are aimed for performing spectral decomposition of the magnetograms 
 (e.g., using Zernike polynomials)

The Python script util.py has necessary tools to run prepare_data.py and AR_cutout.py
