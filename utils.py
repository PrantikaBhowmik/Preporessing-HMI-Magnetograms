"""
    Utilities.

    PB -- 27 March 2024
    
    necessary tools to run prepare_data.py and AR_cutout.py
"""
import os
from natsort import natsorted

#------------------------------------------------------------------------------
def compile_f90(codepath='./fortran/'):
    """
    Compile the fortran tracing code.
    """
    oldd = os.getcwd()
    os.chdir(codepath)
    #os.system('make clean')
    os.system('make')
    os.chdir(oldd)

# A function that returns the length of the value:
def myFunc(e):
  return len(e)


#------------------------------------------------------------------------------
def list_snaps(path, prefix='bz'):
    """
    List all netcdf snapshot files in given directory.
    """
    
    snaps = []
    for file in os.listdir(path):
        if (file.startswith(prefix)):
            if (len(file)==len(prefix)+7 or len(file)==len(prefix)+8 or len(file)==len(prefix)+9):
                snaps.append(file)
    snaps = natsorted(snaps)

    return snaps

#------------------------------------------------------------------------------
def make_im_dir(path, name='im'):
    """
    Make a subdirectory in a run directory.
    """
    if not os.path.exists(path+name):
        os.makedirs(path+name)

