#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan 10 15:57:15 2019

@author: jishnu
"""
# =============================================================================
# This is part one of the reduction script.
# this script is placed inside a folder which sits in the directory containing
# all the FITS files
#
# Science files, Bias files and FLat files.
#
#
#> This script runs through the directory; collects and segregates all files 
# into sci, bias and flat.
#> Computes masterbias and masterflats
#
#> Calibrates all the files and writes to a new sub-directory with a .proc.fits
# extension
# =============================================================================

#%%
from astropy.io import fits
from astropy.visualization import astropy_mpl_style
from astropy.stats import sigma_clipped_stats
import astroalign as aa
import matplotlib.pyplot as plt
import numpy as np
import glob
import os

#%%
# =============================================================================
# creating empty lists/containers to store values
# =============================================================================
files,fitsf,obj,ut,ra,dec,airmass,obsdate = [],[],[],[],[],[],[],[]
exp,gain,rdnoise,filtr = [],[],[],[]
hjd,bjd = [],[] #for ut to HJD conversion
bfiles,ffiles = [],[]
bias,flat = [],[]
masterbias,masterflat = 0,0
y,x = 0,0
pattern = ''
#%%
class fread():

    def fitsfiles(): #find all image files and add to a list
        for fitsName in glob.glob('*.fits'):
            files.append(fitsName)

    def biasfiles(): #find all bias files and add to a list
        for Bfile in glob.glob('*bia*.fits') or glob.glob('*Bia*.fits'):
            bfiles.append(Bfile)
            
    def flatfiles(): #find all flat files and add to a list
        for Ffile in glob.glob('*sflat_V*.fits') or glob.glob('*Flat*.fits'):
            ffiles.append(Ffile)


flread = '''fread.fitsfiles()
fread.biasfiles()
fread.flatfiles()
'''

#%%
def calibration():
    
    global masterbias, masterflat,y,x
    
    y,x = get_img_dims(bfiles[0])[1],get_img_dims(bfiles[0])[2] 
    # x,y dimensions of the images
    print('Found %d bias files'%len(bfiles)) #number of bias files
    print('Computing masterbias')
    # A 3D array to hold all bias images with images as a X*Y array
    biasImages = np.zeros((y, x, len(bfiles)))
    for i in range(len(bfiles)): # load all images into the array above.
        biasImages[:,:,i] = get_img(bfiles[i])[0]
    # Take median of the array along the Z axis (x=0,y=1,Z=2)
    masterbias = np.median(biasImages, axis=2) #updates global 'masterbias'
    del biasImages #del the array & save memory
    
    
    print('Found %d flat files'%len(ffiles))
    print('Computing masterflat')
    # A 3D array to hold all flat images with images as a X*Y array
    flatImages = np.zeros((y, x, len(ffiles)))
    for i in range(len(ffiles)): # load all images into the array above.
        flatImages[:,:,i] = get_img(ffiles[i])[0]  
    # Take median of the array along the Z axis (x=0,y=1,Z=2)
    masterflat = np.median(flatImages, axis=2) #updates global 'masterflat'
    masterflat = masterflat - masterbias
    masterflat = masterflat/np.mean(masterflat)
    
    del flatImages #del the array & save memory
    
    print('Calculating bias and flat corrections')

#%%
def align(self): #magic!!
    try:
        self = np.log10(self)
        p, (pos_self,pos_refimage) = aa.find_transform(self,refimage)
        return p #transformation parameters
    except RuntimeWarning:
        self += abs(np.amin(self)+10)
        self = np.log10(self)
        p, (pos_self,pos_refimage) = aa.find_transform(self,refimage)
        return p #transformation parameters
    

#%%
def get_img(self): #smart way #returns 2D array from fits file data
    hdul = fits.open(self)
    image = np.asarray(hdul[0].data[0])
    #image = np.array(image).astype(np.float16)
    #reads data section of fits into an array as float
    hdul.close() #closing the file. GOTTA save RAM!
    return image

#%%
''' not very useful now
    this is here for sentimental value'''
def get_img_slow(self): #returns 2D array from fits file data
    hdul = fits.open(self)
    image = np.asarray(hdul[0].data).astype(float) #reads data section of 
                                                   #fits into an array as float
    #image = np.array(image) #converting image to an array
    #usually fits data comes as a 3D array, of which we need only one layer.
    y,x = image.shape[-2],image.shape[-1] #getting y and x axes
    image = image.reshape(y,x) #converting 3D array to 2D
    hdul.close() #closing the file. GOTTA save RAM!
    return image

#%%
def get_img_dims(self): #returns an array from fits file data
    hdul = fits.open(self)
    image = np.asarray(hdul[0].data)
    y,x = image.shape[-2],image.shape[-1] #getting y and x axes
    hdul.close() #closing the file. GOTTA save RAM!
    return image,y,x

#%%
os.chdir('../') #going up one directory in folder tree
file_list = '''files,bfiles,ffiles = [],[],[]'''
exec(file_list) #declare all lists
exec(flread)
files = list(set(files)-set(bfiles))
files = list(set(files)-set(ffiles))
files.sort()
refimage = get_img_slow(files[0]) #file to be used as reference frame
refimage += 70
refimage = np.log10(refimage)


#%%
def processed():
    
    folder = os.getcwd() #get the current working directory
    '''
    Now the program should be in the directory where the FITS files are stored
    '''
    directory = 'reduced' #a directory name to store processed files
    
    if not os.path.exists(directory):
        print('create a new directory %s'%directory) #create a new directory
        os.makedirs(directory)              #to write all the reduced files.
    else:
        print('''!! %s already exists..!!'''%directory)
        
    procFolder = (folder+'/reduced/') #path to the directory 'reduced'
    
    #file_list = '''files,bfiles,ffiles = [],[],[]'''
    #exec(file_list) #declare all lists
    #exec(flread) #read through the directory using glob and find, sci, 
                 #bias and flat files
    calibration()
    
    for i in range(len(files)): #loop through the files and create reduced file
        rawHDU = fits.open(files[i]) #open a Sci file
        rawDATA = get_img(files[i]) #Image part of the file
        rawHEADER = rawHDU[0].header #header of the image
        #print(files[i]) 
        
        name = str(files[i]).strip('.fits') #remove fits extension
        try:
            #p = align(rawDATA) #getting transformation
            #rawDATA = aa.apply_transform(p,rawDATA,refimage) #apply transform
            procDATA = (rawDATA-masterbias)/masterflat #one step, sweet and simple
            procHDU = fits.PrimaryHDU(procDATA) #new file, processed
            procHDU.header = rawHEADER #replacing header
        except:
            print('oops file %s didnt meet alignment criteria'%files[i])
            print('removing file %s from disk'%files[i])
            os.remove(files[i])
            continue
        #y1,y2,x1,x2 = 100,y-100,100,x-100 #excluding a 100 pixel border
        #rawDATA = rawDATA[y1:y2, x1:x2].astype(float) #taking a section

        try:
            procHDU.header.remove('BZERO')
            procHDU.header.remove('BSCALE')
        except:
            print('BZERO & BSCALE keywords are not present')
        #writing the reduced file to disk
        print('writing %s to disk'%(name+'.proc.fits')) #seeing is believing
        procHDU.writeto(procFolder+name+'.proc.fits', overwrite=True)
    
    status = input('>>> What is 0/1? in python?\n>>> ')
    if status == '0.0':
        print('good boi!!')
        print('Done!!, time to run the photometry script!')
    else:
        print('''you shouldn't be even allowed to touch a computer!''')
        print('but who am I to say anything?...')
        print('''D...Done..!''')

#%%
def imgshow(self): #display the image
    plt.figure(figsize=(8,8)) #size
    plt.style.use(astropy_mpl_style) #styling and gimmicks of astropy
    
    mean, median, sigma = sigma_clipped_stats(self) #use this for contrast 
                                                    #stretching
    plt.imshow(self,cmap='gray',vmin=median-6*sigma, vmax=median+6*sigma)
    plt.colorbar() #a bar showing the scale of the image
    plt.show()
    plt.close() #clean up after yourself!!
#%%
processed()
