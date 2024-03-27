"""
 PB -- 27 March 2024
 
 Python script to download three components of HMI SHARP patches
 using hmi.sharp_ces_720s data series, 
 the chosen time for each patch is when it's closest to the central merideian
 thus we can approximate, Bp: B_phi -> Bx, Bt: B_theta -> -By, Br -> Bz
 # read: http://hmi.stanford.edu/doc/magnetic/VectorTransformation.pdf
 # read: http://hmi.stanford.edu/doc/magnetic/VectorCoordinates.pdf
 # read: http://hmi.stanford.edu/hminuggets/?p=1757
 # read: http://jsoc.stanford.edu/doc/data/hmi/sharp/sharp.htmv (CEA VERSION - Unique Segments)
"""
from __future__ import division, print_function
import os.path
import matplotlib.pyplot as plt
from astropy.io import fits
import drms
import numpy as np
import utils
import os
import sys

path = '/%%/%%/%%/SHARPs/' # directory path where the data will be downloaded
series = 'hmi.sharp_cea_720s'
segments = ['Bp', 'Bt', 'Br', 'magnetogram'] 
kwlist = ['T_REC', 'LON_FWT', 'OBS_VR', 'CROTA2',
          'CRPIX1', 'CRPIX2', 'CDELT1', 'CDELT2', 'CRVAL1', 'CRVAL2']
          


def read_fits_data(fname):
    """Reads FITS data and fixes/ignores any non-standard FITS keywords."""
    hdulist = fits.open(fname)
    hdulist.verify('silentfix+warn')
    return hdulist[1].data

sharps = np.loadtxt(path+'SHARPs.txt') # .txt file with the numbers odf the HMI-SHARPs regions
sharpsnum = sharps[:,0]
sharpsn = len(sharpsnum)

print('----------------------------')
print(path)
print('Number of SHARPs to study %i' % len(sharpsnum))
print('----------------------------') 

datar =  np.zeros((sharpsn,5)) # details regarding X-Y dimension of the Sharp patch in Mm unit
 
for i in range(0, len(sharpsnum)):
    sharpnum = int(sharpsnum[i])
    utils.make_im_dir(path, name='SHARP_'+str(sharpnum)+'/')
    path1 = path +'SHARP_'+str(sharpnum)+'/'  # separate folder with each SHARP region
    print('folder name: ', path1)
     
    # Create DRMS client and query metadata (use your own registered email address)     
    c = drms.Client(email='%%%%%@gmail.com', verbose=True)
    k = c.query('%s[%d]' % (series, sharpnum), key=kwlist, rec_index=True)
    
    # Find the record that is clostest to the central meridian, 
    # by using the minimum of the patch's absolute longitude:
    rec_cm = k.LON_FWT.abs().idxmin()
    k_cm = k.loc[rec_cm]
    t_cm = drms.to_datetime(k.T_REC[rec_cm])
    print(rec_cm, '@', k.LON_FWT[rec_cm], 'deg')
    print('Timestamp:', t_cm)

    # first check if the files were already downloaded
    t_cm_str = t_cm.strftime('%Y%m%d_%H%M%S_TAI')
    fname_mask = '{series}.{sharpnum}.{tstr}.{segment}.fits'
    fnames = {
        s: fname_mask.format(
            series=series, sharpnum=sharpnum, tstr=t_cm_str, segment=s)
        for s in segments}
    
    # now submit the export request and download all files that were not downloaded in any previous run
    download_segments = []
    for k, v in fnames.items():
        if not os.path.exists(v):
            download_segments.append(k)

    if download_segments:
        exp_query = '%s{%s}' % (rec_cm, ','.join(download_segments))
        r = c.export(exp_query)
        dl = r.download('.')
    
    # now finally read the data from the downloaded FITS files
    mag = read_fits_data(fnames['magnetogram'])
    bz = read_fits_data(fnames['Br'])
    bx = read_fits_data(fnames['Bp'])
    by = -(read_fits_data(fnames['Bt']))

    print('shape: ', np.shape(mag), np.shape(bz), np.shape(bx), np.shape(by))
    print('Maximum values: ', np.max(mag), np.max(bz), np.max(bx), np.max(by))
    
    # save in .txt files in respective folders
    np.savetxt(path1+'bx.txt', bx)
    np.savetxt(path1+'by.txt', by)
    np.savetxt(path1+'bz.txt', bz)

    nx, ny = mag.shape
    Rs = 6.9634e+2 # solar radius in Mm
    dtheta = k_cm.CDELT1*(np.pi/180.0) # approximate value of delta theta or delta phi
    x = np.linspace(0,nx+1,nx)*dtheta
    y = np.linspace(0,ny+1,ny)*dtheta
    print('shape of x and y : ', np.shape(x), np.shape(y))
    xmin = 0
    xmax = nx*dtheta*Rs
    ymin = 0
    ymax = ny*dtheta*Rs
    
    extent = (xmin, xmax,
              ymin, ymax)
              
    print(extent)
    
    datar[i,0] = int(sharpnum)
    datar[i,1] = int(nx)
    datar[i,2] = int(ny)
    datar[i,3] = nx*dtheta*Rs
    datar[i,4] = ny*dtheta*Rs
    print('the multiplicative constant: ', k_cm.CDELT1, dtheta*Rs)


# plotting data          
    plt.rc('mathtext', default='regular')
    plt.rc('image', origin='lower', interpolation='nearest', cmap='gray')

    fig = plt.figure(figsize=(12, 6.75))
    gs = plt.GridSpec(2, 5, width_ratios=[1, 0.1, 0.1, 1, 0.1])
    ax = [fig.add_subplot(gsi, anchor='E')
          for gsi in [gs[0, 0], gs[0, 3], gs[1, 0], gs[1, 3]]]
    cax = [fig.add_subplot(gsi, anchor='W', aspect=15)
          for gsi in [gs[0, 1], gs[0, 4], gs[1, 1], gs[1, 4]]]

    ax_mag, cax_mag = ax[0], cax[0]
    im_mag = ax_mag.imshow(mag/1e3, extent = extent, vmin=-1, vmax=1)
    cb_mag = plt.colorbar(im_mag, cax_mag, label='$B_{\mathrm{los}}$ [kG]')

    ax_bz, cax_bz = ax[1], cax[1]
    im_bz = ax_bz.imshow(bz/1e3, extent = extent, vmin=-1, vmax=1)
    cb_bz = plt.colorbar(im_bz, cax_bz, label='$B_{z}$ [kG]')

    ax_bx, cax_bx = ax[2], cax[2]
    im_bx = ax_bx.imshow(bx/1e3, extent = extent, vmin=-1, vmax=1)
    cb_bx = plt.colorbar(im_bx, cax_bx, label='$B_{x}$ [kG]')

    ax_by, cax_by = ax[3], cax[3]
    im_by = ax_by.imshow(by/1e3, extent = extent, vmin=-1, vmax=1)
    cb_by= plt.colorbar(im_by, cax_by, label='$B_{y}$ [kG]')

    ax_mag.set_xticklabels([])
    ax_bz.set_xticklabels([])
    ax_bz.set_yticklabels([])
    ax_by.set_yticklabels([])
    ax_bx.set_xlabel('X (in Mm)')
    ax_by.set_xlabel('X (in Mm)')
    ax_mag.set_ylabel('Y (in Mm)')
    ax_bx.set_ylabel('Y (in Mm)')


    fig.subplots_adjust(
        left=0.07, bottom=0.10, right=0.94, top=0.95, wspace=0.05, hspace=0.10)
    plt.draw()

    fig.savefig(path1+'sharp_mag_f'+str(sharpnum)+'.pdf', dpi=200)
    fig.savefig(path1+'sharp_mag_f'+str(sharpnum)+'.png', dpi=200)

#plt.show()

    print('----------------------------')
    print('Done: SHARP %i' % (sharpnum))
    print('----------------------------') 

np.savetxt(path+'dimension.txt', datar)
