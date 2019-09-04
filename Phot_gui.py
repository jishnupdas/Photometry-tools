#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  3 16:31:54 2019

@author: jishnu
"""

#%%

import os
import glob
import photutils
import time as zzz
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy import ndimage
from astropy.io import fits
from astropy.time import Time
from astropy.visualization import astropy_mpl_style
from astropy import time, coordinates as coord, units as u


#%%
class Phot:
    '''Object to hold the FITS file and manipulate image/header to perform
    differential photometry on a single source along with a comparison and
    check star
    '''

    def __init__(self,file,image=None,header=None,clck=None,src_coords=None,
                 aperture=None,annulus=None,cutout=None,obj_num=None,
                 var=None,check=None,comp=None):
        '''opens the image and inits variables

        Parameters:
        ----------
        image,header   : the image and header data from the fits file
        clck           : pixel coordinates received from user click input
        src_coords     : coodinates of all the sources to be analyzed (pixels)
        var,chk,comp   : cutouts of variable,comparison and check star

        '''

        self.image      = fits.getdata(file)
        self.header     = fits.getheader(file)
        self.error      = np.sqrt(self.image)
        self.src_coords = src_coords
        self.clck       = []
        self.var        = var
        self.check      = check
        self.comp       = comp
        self.aperture   = aperture
        self.annulus    = annulus
        self.cutout     = 70
        self.obj_num    = 3


    def imgshow(self,image): #display the image
        '''displays the image and does some histogram stretching'''

        plt.imshow(image,vmin=np.median(image)-1*np.std(image),
                   vmax=np.median(image)+5*np.std(image),cmap='gray')
        plt.gca().invert_yaxis()


    def single_click(self,event): ##to read the event of a click
        '''records click events from the plot'''
        if event.button == 3: ##for right click
            self.clck.append((event.xdata,event.ydata))
            #append array with value of coordinates of click
            #print("clicked on (x,y) : (%.2f,%.2f)"%(event.xdata,event.ydata))
        if len(self.clck) >= self.obj_num:
            plt.close()

    def on_key(self,event):
        'close plot when a key is pressed'
        plt.close()


    def img_click(self):
        '''plot the image & connect the click event to the plot'''

        fig, ax = plt.subplots(figsize=(10,10))
        plt.style.use(astropy_mpl_style) #styling template of astropy
        plt.title('Right click on 3 sources')
        ax.imshow(self.image,vmin=np.median(self.image)-1*np.std(self.image),
           vmax=np.median(self.image)+5*np.std(self.image),cmap='gray')
        plt.gca().invert_yaxis()
        cid = fig.canvas.mpl_connect('button_press_event', self.single_click)
        plt.show()


    def get_center(self,coords):
        '''get the peak of the given cutout'''

        x,y     = int(coords[0]),int(coords[1])

        section = ndimage.gaussian_filter(
                                self.image[y-self.cutout:y+self.cutout,
                                           x-self.cutout:x+self.cutout],
                                           sigma=2)
        #gets a cutout section of the image and applies a gaussian filter
        #to it, and then proceeds to find the maxima, the filter helps to
        #remove random spikes and hence solution converges to maxima

        c = ndimage.maximum_position(section)
        X,Y = x,y
        Y  += c[0]-self.cutout
        X  += c[1]-self.cutout
#        x,y = int(X),int(Y)
        CRD = (int(X),int(Y))

        return CRD,(X,Y),c


    def get_section(self,x,y):
        section = self.image[y-self.cutout:y+self.cutout,
                             x-self.cutout:x+self.cutout]
        return section

    def plot_src(self,coords):

        positions      = coords
        radius         = self.find_optimal_aperture()
        self.aperture  = photutils.CircularAperture(positions, r=radius)
        self.annulus   = photutils.CircularAnnulus(positions, r_in=radius+8,
                                                              r_out=radius+16)
        fig = plt.figure(figsize=(10,10))
        plt.style.use(astropy_mpl_style)
        cid = fig.canvas.mpl_connect('key_press_event', self.on_key)
        self.imgshow(self.image)
        self.aperture.plot(color='white', lw=0.5)
        self.annulus.plot(color='red', lw=0.5)
        plt.xlabel('Press any key to close.')

        plt.show()

    def find_optimal_aperture(self):
        '''function to fin the optimal aperture'''
        radii    = list(np.linspace(1,15,50))
        flx      = []
        position = self.src_coords[0]

        for r in radii:
            aperture = photutils.CircularAperture(position,r)
            tbl      = photutils.aperture_photometry(self.image,
                                                     aperture)
            flx.append(tbl['aperture_sum'][0])

        slope = [dflx/dr for dflx,dr in zip(flx,radii)]

        'plots the curve of growth and identifies the inflection point'
        'by taking the numerical differentiation of the cog'
        'the differntial wiil max out at FWHM'
        'optimal aperture then becomes 2 x FWHM radius'
#        plt.plot(radii,flx,'-o')

        opt_rad = radii[slope.index(max(slope))]*2

        return opt_rad


    def apphot(self):
        '''perform aperture photometry on the selected sources'''

        positions      = self.src_coords
        radius         = self.find_optimal_aperture()
        self.aperture  = photutils.CircularAperture(positions, r=radius)
        self.annulus   = photutils.CircularAnnulus(positions, r_in=radius+6,
                                                         r_out=radius+8)

        tbl       = photutils.aperture_photometry(self.image,
                                                  [self.aperture,self.annulus],
                                                  error=self.error)

        bg_mean   = tbl['aperture_sum_1'] / self.annulus.area
        bg_flx    = bg_mean * self.aperture.area
        final_flx = tbl['aperture_sum_0'] - bg_flx

        tbl['res_flx'] = final_flx
        tbl['flx_err'] = tbl['aperture_sum_err_0']+tbl['aperture_sum_err_1']
        tbl['SNR']     = tbl['res_flx']/bg_flx

        for col in tbl.colnames:
            tbl[col].info.format = '%.6g'
            # for consistent table output

        self.var  = (tbl['res_flx'][0],tbl['flx_err'][0],tbl['SNR'][0])
        self.comp = (tbl['res_flx'][1],tbl['flx_err'][1],tbl['SNR'][1])
        self.check= (tbl['res_flx'][2],tbl['flx_err'][2],tbl['SNR'][2])

        return tbl

    def get_info(self):
        '''get information from header'''
        ut=(self.header['UT'])
#        ra=(self.header['RA'])
#        dec=(self.header['DEC'])
        am=(self.header['AIRMASS'])
        t = Time(self.header['DATE-OBS'], scale='utc')
        jd=(t.jd)

        return ut,am,jd

#%%
os.chdir('../reduced/')

files = glob.glob('*proc.fits')
#%%
'testing code'
#phot   = Phot(files[0])
#
#phot.img_click()
#
#phot.plot_src(phot.clck)
#
#phot.src_coords = [phot.get_center(c)[0] for c in phot.clck]
#
#phot.plot_src(phot.src_coords)
#
#print(phot.clck)
#print(phot.src_coords)
#
#phot_table = phot.apphot()
#
#print(phot_table)
#%%
'some more testing'
#phot   = Phot(files[0])
#phot.src_coords = [(964, 1090), (408, 726), (835, 328)]
#
#phot_table = phot.apphot()

#%%

print('File,UT,AIRMASS,JD,v_flx,err_v,c1,err_c1,c2,err_c2,SNR,flx,err')

for f in files:
    phot       = Phot(f)
    phot.img_click()
    phot.src_coords = [phot.get_center(c)[0] for c in phot.clck]
    phot.plot_src(phot.src_coords)
    plt.close()
    phot_table = phot.apphot()
    ut,am,jd   = phot.get_info()
    v,c1,c2    = phot.var[0],phot.comp[0],phot.check[0]
    ev,ec1,ec2 = phot.var[1],phot.comp[1],phot.check[1]
    snr        = phot.var[2] #SNR of the source

    flx,err    = c1-v,ev+ec1
    print('%s,%s,%.2f,%.6f,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%.2f,%.3f,%.3f'%(
            f,str(ut),am,jd,v,ev,c1,ec1,c2,ec2,snr,flx,err))
