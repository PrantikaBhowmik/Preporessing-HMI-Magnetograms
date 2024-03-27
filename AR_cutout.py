"""
 PB -- 27 March 2024
 
 Primarily to be used on HMI-SHARP magnetogram patches generated by the 
 associated prepare_data.py, which downloads three magnetic field components (Bx, By, Bz) 
 of HMI SHARP patches using hmi.sharp_ces_720s data series
 
 This script processes magnetograms to 
 1. extract the major structures from a given magnetic field distribution on a Cartesian plane
 2. It also calculates the tilt angle created between the positive and negative polarities
 of a given magnetogram, and rotate it such that the tilt will become zero
 3. Add extra padding with zero field strength, such that the centroid of the magnetogram 
 coincides with the origin of the new coordinates
 4. These preprocessing are aimed at performing spectral decomposition of the magnetograms 
 (e.g., using Zernike polynomials)
 
"""
import matplotlib.pyplot as plt
import numpy as np
import math
import utils
import os
import sys
from scipy import ndimage
from scipy.interpolate import interp2d, griddata

path = '/%%/%%/%%/SHARPs/' # directory path where the data will be downloaded

sharps = np.loadtxt(path+'SHARPs.txt')
sharpsnum = sharps[:,0]
sharpsn = len(sharpsnum)

print('----------------------------')
print(path)
print('Number of SHARPs to study %i' % len(sharpsnum))
print('----------------------------') 

# function calculating the centre of mass (total, positive nad negative)
def f(b):
     return ndimage.center_of_mass(np.abs(b))
 
def padding(b,btm,top,lft,rgt):
# padding around the boundary when not enough data is available 
     return np.pad(b, pad_width=((btm, top), (lft, rgt)), mode = 'linear_ramp', end_values = ((0,0), (0,0)))
 
def padding_c(b,pad_wx,pad_wy):
    return np.pad(b, pad_width=((pad_wy, pad_wy), (pad_wx, pad_wx)), mode='constant', constant_values=0)
    # pad_wy : top/bottom extensions, pad_wx: left/right extensions

#########################################
# parameters to play with
bthresh = 40 # in Gauss, a threshold to discard very low intensity field
fc = 1.0 # decides the multiple of dx and dy on the both side of the CM in the final data patch
fcl = 0.9 # controls the colour saturation in the colourbars
s = 7.0 # sigma: level of gaussian smoothing
r = 0.5 # what portion of the structure with maximum area should be considered
dxcor = 0.0003 # scaling to convert grid index to X-Y coordinate values
dycor = dxcor
#########################################


dtls = np.zeros((14,1))
tm = 0
 
for i in range(0, len(sharpsnum)):
    sharpnum = int(sharpsnum[i])
    path1 = path +'SHARP_'+str(sharpnum)+'/'  
    print('folder name: ', path1)

    utils.make_im_dir(path1, name='analysis/')
    path2 = path1 + 'analysis/'  
    print('folder name: ', path2)
    
    bx = np.loadtxt(path1+'bx.txt')
    by = np.loadtxt(path1+'by.txt')
    bz = np.loadtxt(path1+'bz.txt')
    xn = int(np.size(bz,0))
    yn = int(np.size(bz,1))
    
    bx1 = bx
    by1 = by
    bz1 = bz
    print('----------------------------')
    print(' ')
    # preparing for the zernike decomposition after tilt correction based on bz distribution only
    print('PREPARE data for Zernike decomposition')

    a2 = bz1
    a2 = np.transpose(a2)
    ax2 = np.transpose(bx1)
    ay2 = np.transpose(by1)
    # add padding around the data
    wdt = max(int(xn), int(yn))
    pad_wx = int(wdt*0.5) # same number of arrays with zeros are addeded in left and right
    pad_wy = int(wdt*0.5) # same number of arrays with zeros are addeded in top and bottom
    an = padding_c(a2,pad_wx,pad_wy)
    axn = padding_c(ax2,pad_wx,pad_wy)
    ayn = padding_c(ay2,pad_wx,pad_wy)

    
    a1 = ndimage.gaussian_filter(an, sigma=s, order=0)
    xn1 = np.size(an,0)
    yn1 = np.size(an,1)
    a1_masked = a1*0.0
    mask1 = a1*0.0
    print('xn, yn and shape of a2: ', xn, yn, np.shape(a2))
    print('xn1, yn1 and shape of the padded an: ', xn1, yn1, np.shape(an))
    
    ##### creating the mask #########
    mask = np.abs(a1) > bthresh
    label_im, nb_labels = ndimage.label(mask)
    print('number of structures', nb_labels)
    sizes = ndimage.sum(mask, label_im, range(nb_labels + 1))
    mask_size = sizes < r*np.max(sizes) # criterion of selection

    remove_pixel = mask_size[label_im]
    label_im[remove_pixel] = 0
    np.savetxt(path2+'maximum_mask_bz.txt', label_im)
    
    k = 0
    xlst = []
    ylst = []
    for i in range(0,int(xn1)):
        for j in range(0,int(yn1)):
            if(label_im[i,j]>0):
                mask1[i,j] = 1.0
                a1_masked[i,j]=an[i,j]
                xlst.append(i)
                ylst.append(j)
                k = k+1
    
    # decide the X-Y extent
    xlst = np.array(xlst)
    ylst = np.array(ylst)
    xmin = np.min(xlst)
    xmax = np.max(xlst)
    ymin = np.min(ylst)
    ymax = np.max(ylst)
    print('x limit: ', xmin, xmax, 'y limit: ', ymin, ymax)
    
    # determine the centroid (centre of mass calculation)
    [xc, yc] = f(an[xmin:xmax,ymin:ymax])
    xc = xc + xmin
    yc = yc + ymin
    ap2 = an*0.0
    an2 = an*0.0
    for i in range(0,int(xn1)):
        for j in range(0,int(yn1)):
            if(an[i,j] > 0.0):
                ap2[i,j] = an[i,j]
            else:
                an2[i,j] = an[i,j]
    
    # determine the centroid for positive field (centre of mass calculation)
    [xpc, ypc] = f(ap2[xmin:xmax,ymin:ymax])
    xpc = xpc + xmin
    ypc = ypc + ymin
    [xnc, ync] = f(an2[xmin:xmax,ymin:ymax])
    xnc = xnc + xmin
    ync = ync + ymin
    
    print('amplitude at CM: ', xc, yc, an[int(xc),int(yc)])
    print('amplitude at positive CM: ', xpc, ypc, an[int(xpc),int(ypc)])
    print('amplitude at negative CM: ', xnc, ync, an[int(xnc),int(ync)])

    # rotate an coordinates for tilt corrections
    # asign coordinate values to the grids
    
    xcor = np.zeros((xn1, yn1))
    ycor = np.zeros((xn1, yn1))
    
    for i in range(0, int(xn1)):
        for j in range(0, int(yn1)):
            xcor[i,j] = (i-int(xc))*dxcor
            ycor[i,j] = (j-int(yc))*dycor
        
    print('cell width along X and Y coordinates: ', dxcor, dycor)
    xcor1 = np.zeros((xn1, yn1))
    ycor1 = np.zeros((xn1, yn1))
    # decide the sign of theta, i.e., direction of rotation
    theta = math.atan(np.abs(ypc-ync)/np.abs(xpc-xnc))
    theta1 = 0.0
    if (xpc<xnc and ypc>ync):
        theta1 = theta
        print('rotation: anti-clockwise, -theta')
    if (xpc<xnc and ypc<ync):
        theta1 = -theta
        print('rotation: clockwise, theta')
    if (xpc>xnc and ypc>ync):
        theta1 = (math.pi - theta)
        print('rotation: anti-clockwise, pi - theta')
    if (xpc>xnc and ypc<ync):
        theta1 = (math.pi + theta)
        print('rotation: anti-clockwise, pi + theta')

    print('tilt angle: ', np.rad2deg(theta1))
    axn1 = axn*0.0
    ayn1 = ayn*0.0
    
    # calculating the components in new rotated coordinates
    for i in range(0, xn1):
        for j in range(0, yn1):
            xcor1[i,j] = xcor[i,j]*math.cos(theta1) - ycor[i,j]*math.sin(theta1)
            ycor1[i,j] = xcor[i,j]*math.sin(theta1) + ycor[i,j]*math.cos(theta1)
            axn1[i,j] = axn[i,j]*math.cos(theta1) + ayn[i,j]*math.sin(theta1)  
            ayn1[i,j] = -axn[i,j]*math.sin(theta1) + ayn[i,j]*math.cos(theta1)
               
    print('test shapes: ', np.shape(xcor1), np.shape(axn1))  
            
    
    aznf = an
    axnf = axn1
    aynf = ayn1
    
    # now decides the range of data to be saved 
    # distance between centroid and X-Y extent
    dx = np.max([np.abs(xmin-int(xc)), np.abs(xmax-int(xc))])
    dy = np.max([np.abs(ymin-int(yc)), np.abs(ymax-int(yc))])
    print('extents from CM: ', dx, dy)
    # take equal extents along X and Y based on dx and dy 
    xminf = int(int(xc) - int(dx*fc))
    xmaxf = int(int(xc) + int(dx*fc))
    yminf = int(int(yc) - int(dy*fc))
    ymaxf = int(int(yc) + int(dy*fc))
    print('new x limit: ', xminf, xmaxf, 'new y limit: ', yminf, ymaxf)

    dtls[0] = int(sharpnum) # AR number
    dtls[1] = tm # time for saving data
    dtls[2] = int(xmaxf-xminf) # number of data points along X
    dtls[3] = int(ymaxf-yminf) # number of data points along Y
    dtls[4] = int(xc) # X-location of CM of abs(Bz)
    dtls[5] = int(yc) # Y-location of CM of abs(Bz)
    dtls[6] = int(xpc) # X-location of CM of positive(Bz)
    dtls[7] = int(ypc) # Y-location of CM of positive(Bz)
    dtls[8] = int(xnc) # X-location of CM of negative(Bz)
    dtls[9] = int(ync) # Y-location of CM of negative(Bz)
    dtls[10] = np.rad2deg(theta) # tilt angle
    dtls[11] = np.rad2deg(theta1) # the actual angle after rotation
    dtls[12] = s # level of gausian smoothing, sigma
    dtls[13] = dxcor # the contant to be multiplied to give ccordinate values to the grids

    # saving the data within this new limits
    np.savetxt(path2+'x_mcor.txt',xcor1[xminf:xmaxf,yminf:ymaxf])
    np.savetxt(path2+'y_mcor.txt',ycor1[xminf:xmaxf,yminf:ymaxf])
    np.savetxt(path2+'Bzm_f.txt',aznf[xminf:xmaxf,yminf:ymaxf])
    np.savetxt(path2+'Bxm_f.txt',axnf[xminf:xmaxf,yminf:ymaxf])
    np.savetxt(path2+'Bym_f.txt',aynf[xminf:xmaxf,yminf:ymaxf])
    np.savetxt(path2+'Rdetails_max.txt',dtls)

    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2, 3, figsize=(16, 8))
    levels = np.linspace(-np.max(np.abs(a2))*fcl, np.max(np.abs(a2))*fcl, 20)
    pl1 = ax1.contourf(np.transpose(a2), cmap='bwr', levels = levels)
    plt.colorbar(pl1, ax = ax1)
    ax1.title.set_text('original')

    pl2 = ax2.contourf(xcor, ycor, an, cmap='bwr', levels = levels)
    plt.colorbar(pl2, ax = ax2)
    ax2.scatter((xpc-xc)*dxcor,(ypc-yc)*dycor, c = 'm', marker = "+")
    ax2.scatter((xnc-xc)*dxcor,(ync-yc)*dycor, c = 'c', marker = "_")
    ax2.scatter((xc-xc)*dxcor,(yc-yc)*dycor, c = 'k', marker = "*")
    ax2.title.set_text('Bz distribution with CMs')

    pl3 = ax3.contourf(np.transpose(mask1), cmap='gray')
    plt.colorbar(pl3, ax = ax3)
    ax3.plot(xmin,ymin,'yo')
    ax3.plot(xmin,ymax,'go')
    ax3.plot(xmax,ymin,'mo')
    ax3.plot(xmax,ymax,'co')
    ax3.title.set_text('filtered cutout, time')

    pl4 = ax4.contourf(np.transpose(xcor1[xminf:xmaxf,yminf:ymaxf]), np.transpose(ycor1[xminf:xmaxf,yminf:ymaxf]), np.transpose(axnf[xminf:xmaxf,yminf:ymaxf]), cmap='bwr', levels = levels)
    plt.colorbar(pl4, ax = ax4)
    ax4.title.set_text('Bx: tilt correction, time')

    pl5 = ax5.contourf(np.transpose(xcor1[xminf:xmaxf,yminf:ymaxf]), np.transpose(ycor1[xminf:xmaxf,yminf:ymaxf]), np.transpose(aynf[xminf:xmaxf,yminf:ymaxf]), cmap='bwr', levels = levels)
    plt.colorbar(pl5, ax = ax5)
    ax5.title.set_text('By: tilt correction')

    pl6 = ax6.contourf(np.transpose(xcor1[xminf:xmaxf,yminf:ymaxf]), np.transpose(ycor1[xminf:xmaxf,yminf:ymaxf]), np.transpose(aznf[xminf:xmaxf,yminf:ymaxf]), cmap='bwr', levels = levels)
    plt.colorbar(pl6, ax = ax6)
    ax6.title.set_text('Bz: tilt correction')

    plt.savefig(path2+'processed_SHARP_'+str(int(sharpnum))+'.png', bbox_inches='tight')
    #plt.show()
    plt.close()
    print('----------------------------')
    print(' ')

    print('----------------------------')
    print('Done: SHARP %i' % (sharpnum))
    print('----------------------------')
