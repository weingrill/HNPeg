#!/usr/bin/python
# -*- coding: utf-8 -*-
'''
Created on May 13, 2016

@author: Joerg Weingrill <jweingrill@aip.de>
'''

import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
       
def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    http://stackoverflow.com/questions/2413522/weighted-standard-deviation-in-numpy
    """
    import math
    
    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2, weights=weights)  # Fast and numerically precise
    return (average, math.sqrt(variance))

def storespectra():
    from glob import glob
    import shutil

    filenames = glob('/work2/jwe/Projects/HNPeg/data/science*.txt')
    filenames = sorted(filenames)
    print len(filenames)
    nspectra = len(filenames)
    nwavelengths = 30000
    
    wavelengths = np.linspace(4000., 6800., nwavelengths)
    a = np.zeros((nspectra,nwavelengths))
    
    reference = np.genfromtxt('/work2/jwe/Projects/HNPeg/results/reference.txt')
    
    for i,filename in enumerate(filenames):
        data = np.genfromtxt(filename)
        wv = data[:,0]
        fl = data[:,1]
        er = data[:,2]
        #plt.plot(wv, fl, 'b') 
        #plt.xlim(4000., 8000.)
        #plt.ylim(0., 1.)
        #plt.title(filename)
        
        # clip fluxes to remove wrong pixels
        cfl = np.clip(fl, 0.0, 1.0)
        int_flux = interp1d(wv, cfl, kind='slinear')
        fluxes = int_flux(wavelengths)
        #plt.plot(wavelengths, fluxes,'g')
        #plt.show()
        
        print filename, np.std(fluxes - reference[:,1])
        if np.std(fluxes - reference[:,1])>0.03:
            shutil.move(filename,'/work2/jwe/Projects/HNPeg/data/bad/')
            #plt.plot(wavelengths, fluxes,'g')
            #plt.plot(wavelengths, fluxes - reference[:,1],'r')
            #plt.xlim(4000., 6800.)
            #plt.title(filename)
            #plt.show()
        
        a[i,:] = fluxes
    
    mflux = np.mean(a, axis=0)
    stdflux = np.std(a, axis=0)
    import pickle
    with open('/work2/jwe/Projects/HNPeg/results/spectra.pickle','wb') as picklefile:
        pickle.dump(a, picklefile)
    plt.plot(wavelengths, mflux, 'b') 
    plt.plot(wavelengths, stdflux, 'r') 
    np.savetxt('/work2/jwe/Projects/HNPeg/results/analysis.txt', np.column_stack((wavelengths, mflux, stdflux)))
    #np.savetxt('/work2/jwe/Projects/HNPeg/results/filenames.txt', filenames)
    plt.show()

def saveimage(filename, myarray):
    """
    First ensure your numpy array, myarray, is normalised with the max value at 1.0.
    Apply the colormap directly to myarray.
    Rescale to the 0-255 range.
    Convert to integers, using np.uint8().
    Use Image.fromarray().
    """
    from PIL import Image  # @UnresolvedImport
    
    mavg = np.nanmean(myarray)
    mstd = 2.0*np.nanstd(myarray)
    
    myarray -= mavg
    i = np.where(np.isnan(myarray))
    myarray[i] = 0.0

    myarray = np.clip(myarray, -mstd, mstd)
    
    myarray *= -1.0
    myarray -= np.nanmin(myarray)
    myarray *= np.trunc(256.0/np.nanmax(myarray))
    
    simg = np.rint(myarray)
    simg = simg.astype('uint8')
    
    im = Image.fromarray(plt.cm.jet(simg, bytes=True))  # @UndefinedVariable
    im.save(filename)


def scalespectrum():
    import pickle
    
    nwavelengths = 30000
    wavelengths = np.linspace(4000., 6800., nwavelengths)
        
    with open('/work2/jwe/Projects/HNPeg/results/spectra.pickle','rb') as picklefile:
        a = pickle.load(picklefile)
        
    reference = np.genfromtxt('/work2/jwe/Projects/HNPeg/results/reference.txt')
    #plt.plot(a[0,:])
    #plt.plot(reference[:,2],'g')
    #plt.plot(a[:,12345],'r')
    #plt.plot(a[:,12346],'b')
    #plt.plot(a[:,24414],'k')
    m = np.mean(a, axis=0)
    print a.shape
    ms = np.resize(m, a.shape)
    a = a - ms
    saveimage('/work2/jwe/Projects/HNPeg/results/spectra.png', a)
    plt.imshow(a, aspect='auto')
    xticks = np.arange(4000., 6800., 10.)
    xlocs = np.arange(0., nwavelengths, 10.)
    plt.xticks(xlocs, xticks)
    plt.xlim(6550., 6574.)
    plt.show()
    
    #plt.plot([0,1],[0,1],'r--')
    
    
       
if __name__ == '__main__':
    #storespectra()
    scalespectrum()