# Preporessing-HMI-Magnetograms
Python scripts to pre-process HMI-SHARP magnetograms with all three magnetic field components (Bx, By and Bz). See the 


The Python script prepare_data.py 
downloads three components of HMI SHARP patches
using hmi.sharp_ces_720s data series (http://jsoc.stanford.edu/ajax/lookdata.html?ds=hmi.sharp_cea_720s), 
the chosen time for each patch is when it's closest to the central meridian. Thus, we can approximate 
Bp: B_phi -> Bx, Bt: B_theta -> -By, Br -> Bz

The Python script 
