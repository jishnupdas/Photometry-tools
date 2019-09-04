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
#from astroquery.vizier import Vizier
#from astropy.coordinates import SkyCoord
#from astropy.stats import sigma_clipped_stats
from astropy.visualization import astropy_mpl_style
from astropy import time, coordinates as coord, units as u

#%%
os.chdir('/home/jishnu/Documents/data_red/200643/09/reduced/')

files = glob.glob('*proc.fits')

#%%
class Phot:
    '''Object to hold the FITS file and manipulate image/header to perform
    differential photometry on a single source along with a comparison and
    check star
    '''

    def __init__(self,file,image=None,header=None,clck=None,src_coords=None,
                 cutout=None,obj_num=None,var=None,check=None,comp=None):
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
        self.cutout     = 100
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
            print("clicked on (x,y) : (%.2f,%.2f)"%(event.xdata,event.ydata))
        if len(self.clck) >= self.obj_num:
            plt.close()


    def img_click(self):
        '''plot the image & connect the click event to the plot'''

        fig, ax = plt.subplots(figsize=(10,10))
        plt.style.use(astropy_mpl_style) #styling template of astropy
        ax.imshow(self.image,vmin=np.median(self.image)-1*np.std(self.image),
           vmax=np.median(self.image)+5*np.std(self.image),cmap='gray')
        plt.gca().invert_yaxis()
        cid = fig.canvas.mpl_connect('button_press_event', self.single_click)
        plt.show()


    def get_center(self,coord):
        '''get the peak of the given cutout'''

        x,y     = int(coord[0]),int(coord[1])

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
        fig = plt.figure(figsize=(10,10))
        plt.style.use(astropy_mpl_style)
        ax = fig.gca()
        self.imgshow(self.image)
        circles = [plt.Circle(
                (coords[i][0], coords[i][1]),
                radius = 30,edgecolor='r',
                facecolor='None') for i in range(len(coords))
                ]
        for c in circles:
            ax.add_artist(c)

        plt.show()

    def find_optimal_aperture(self):
        '''function to fin the optimal aperture'''
        radii    = list(np.linspace(1,19,50))
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

        positions = self.src_coords
        radius    = self.find_optimal_aperture()
        aperture  = photutils.CircularAperture(positions, r=radius)
        annulus   = photutils.CircularAnnulus(positions, r_in=radius+6,
                                                         r_out=radius+8)
        tbl       = photutils.aperture_photometry(self.image,
                                                  [aperture,annulus],
                                                  error=self.error)

        bg_mean   = tbl['aperture_sum_1'] / annulus.area
        bg_flx    = bg_mean * aperture.area
        final_flx = tbl['aperture_sum_0'] - bg_flx

        tbl['res_flx'] = final_flx
        tbl['flx_err'] = tbl['aperture_sum_err_0']+tbl['aperture_sum_err_1']

        for col in tbl.colnames:
            tbl[col].info.format = '%.6g'
            # for consistent table output

        self.var  = (tbl['res_flx'][0],tbl['flx_err'][0])
        self.comp = (tbl['res_flx'][1],tbl['flx_err'][1])
        self.check= (tbl['res_flx'][2],tbl['flx_err'][2])

        return tbl

    def get_info(self):

        ut=(self.header['UT'])
#        ra=(self.header['RA'])
#        dec=(self.header['DEC'])
        am=(self.header['AIRMASS'])
        t = Time(self.header['DATE-OBS'], scale='utc')
        jd=(t.jd)

        return ut,am,jd
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
phot   = Phot(files[0])
phot.src_coords = [(964, 1090), (408, 726), (835, 328)]

phot_table = phot.apphot()


#%%
#positions = [(964, 1090), (408, 726), (835, 328)]
#aperture = photutils.CircularAperture(positions, r=7.)
#data = phot.image
#error = 0.1 * data
#
#phot_table = photutils.aperture_photometry(data, aperture, error=error)
#for col in phot_table.colnames:
#    phot_table[col].info.format = '%.4g'  # for consistent table output
#print(phot_table)
##%%
#def get_img(self):
#    '''returns the image as an np array from fits file'''
#    with fits.open(self) as hdul:
#        image = np.asarray(hdul[0].data)
#    return image
#
#def single_click(event): ##to read the event of a click
#    '''records click events from the plot'''
#    global xclickarr, yclickarr,n
#    if event.button == 3: ##for right click
#        clickarr.append((event.xdata,event.ydata))
#        #append array with value of coordinates of click
#        print("clicked (x,y) : (%.2f,%.2f"%(event.xdata,event.ydata))
#    if len(clickarr) >= n:
#        plt.close()
#
##%%
#def imgshow(self): #display the image
#    plt.figure(figsize=(8,8)) #size
#    plt.style.use(astropy_mpl_style) #styling template of astropy
#    plt.imshow(self,vmin=np.median(self)-1*np.std(self),
#               vmax=np.median(self)+10*np.std(self),cmap='gray')
#    plt.gca().invert_yaxis()
#    #use this for contrast stretching
#
##    plt.colorbar() #a bar showing the scale of the image
##    plt.show()
#
##%%
#def img_click(self): #display the image
#    fig, ax = plt.subplots(figsize=(10,10))
#    plt.style.use(astropy_mpl_style) #styling template of astropy
#    ax.imshow(self,vmin=np.median(self)-1*np.std(self),
#               vmax=np.median(self)+5*np.std(self),cmap='gray')
#    #use this for contrast stretching
#    cid = fig.canvas.mpl_connect('button_press_event', single_click)
#    plt.show()
#
#def dist(p1,p2):
#    x1,y1 = int(p1[0]),int(p1[1])
#    x2,y2 = int(p1[0]),int(p1[1])
#    dist = np.sqrt(x1-x2)**2 +(y1-y2)**2
#    return dist
##%%
#def get_center(image,coord,cutout):
#    x,y  = int(coord[0]),int(coord[1])
#
#    c = ndimage.maximum_position(image[y-cutout:y+cutout,x-cutout:x+cutout])
#    Y,X = y,x
#    Y += cutout-c[0]
#    X += cutout-c[1]
#    y,x = int(Y),int(X)
#    CRD = (x,y)
#
#    return CRD,(X,Y),c
#
##%%
#def get_c(image,coord,cut):
#    x,y  = int(coord[0]),int(coord[1])
#    #print(y,x)
#    c = np.where(image[y-cut:y+cut,
#                       x-cut:x+cut] == np.amax(image[y-cut:y+cut,x-cut:x+cut]))
#
#    Y,X = y,x
#
#    Y += 50-c[0]
#    X += 50-c[1]
#
#    return c,(X,Y)
#
#
#
#def plot_src(self,coords):
#    fig = plt.figure(figsize=(10,10))
#    plt.style.use(astropy_mpl_style)
#    ax = fig.gca()
#    #plt.imshow(data, vmin=median-3*sigma, vmax=median+3*sigma,cmap='gray')
#    plt.imshow(self,vmin=np.median(self)-1*np.std(self),
#               vmax=np.median(self)+5*np.std(self),cmap='gray')
#    circles = [plt.Circle(
#            (coords[i][0], coords[i][1]),
#            radius = 30,edgecolor='r',
#            facecolor='None') for i in range(len(coords))
#            ]
#    for c in circles:
#        ax.add_artist(c)
#
#    plt.show()
#
##%%
#def test_plot(image,coord,correc):
#    imgshow(image)
#    plt.plot(coord[0],coord[1],'ro')
#    plt.plot(correc[0],correc[1],'go')
##%%
#clickarr,n = [],1
#
#img_click(get_img(files[0]))
#
#c_cord = []
#
#for c in clickarr:
#    image = get_img(files[0])
#    CRD,cent,cord = get_center(image,c,50)
#    c_cord.append(cord)
#
#
#print('clicked array',clickarr)
#print('corrected coordinates',c_cord)
#
#plot_src(get_img(files[0]),clickarr)
#plot_src(get_img(files[0]),c_cord)
