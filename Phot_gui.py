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
#from astropy.time import Time
#from astroquery.vizier import Vizier
#from astropy.coordinates import SkyCoord
#from astropy.stats import sigma_clipped_stats
from astropy.visualization import astropy_mpl_style
#from astropy import time, coordinates as coord, units as u

#%%
os.chdir('/home/jishnu/Documents/data_red/200643/09/reduced/')

files = glob.glob('*proc.fits')

#%%
#class Phot:
#    '''Object to hold the FITS file and manipulate image/header'''
#
#    def __init__(self,image=None,header=None,sources=None,)

#%%
def get_img(self):
    '''returns the image as an np array from fits file'''
    with fits.open(self) as hdul:
        image = np.asarray(hdul[0].data)
    return image

def single_click(event): ##to read the event of a click
    '''records click events from the plot'''
    global xclickarr, yclickarr,n
    if event.button == 3: ##for right click
        clickarr.append((event.xdata,event.ydata))
        #append array with value of coordinates of click
        print("clicked (x,y) : (%.2f,%.2f"%(event.xdata,event.ydata))
    if len(clickarr) >= n:
        plt.close()

#%%
def imgshow(self): #display the image
    plt.figure(figsize=(8,8)) #size
    plt.style.use(astropy_mpl_style) #styling template of astropy
    plt.imshow(self,vmin=np.median(self)-1*np.std(self),
               vmax=np.median(self)+10*np.std(self),cmap='gray')
    plt.gca().invert_yaxis()
    #use this for contrast stretching

#    plt.colorbar() #a bar showing the scale of the image
#    plt.show()

#%%
def img_click(self): #display the image
    fig, ax = plt.subplots(figsize=(10,10))
    plt.style.use(astropy_mpl_style) #styling template of astropy
    ax.imshow(self,vmin=np.median(self)-1*np.std(self),
               vmax=np.median(self)+5*np.std(self),cmap='gray')
    #use this for contrast stretching
    cid = fig.canvas.mpl_connect('button_press_event', single_click)
    plt.show()

def dist(p1,p2):
    x1,y1 = int(p1[0]),int(p1[1])
    x2,y2 = int(p1[0]),int(p1[1])
    dist = np.sqrt(x1-x2)**2 +(y1-y2)**2
    return dist
#%%
def get_center(image,coord,cutout):
    x,y  = int(coord[0]),int(coord[1])

    c = ndimage.maximum_position(image[y-cutout:y+cutout,x-cutout:x+cutout])
    Y,X = y,x
    Y += cutout-c[0]
    X += cutout-c[1]
    y,x = int(Y),int(X)
    CRD = (x,y)

    return CRD,(X,Y),c

#%%
def get_c(image,coord,cut):
    x,y  = int(coord[0]),int(coord[1])
    #print(y,x)
    c = np.where(image[y-cut:y+cut,
                       x-cut:x+cut] == np.amax(image[y-cut:y+cut,x-cut:x+cut]))

    Y,X = y,x

    Y += 50-c[0]
    X += 50-c[1]

    return c,(X,Y)


def get_cutout(image,x,y,size):
    cutout = image[y-size:y+size,x-size:x+size]
    return cutout

def plot_src(self,coords):
    fig = plt.figure(figsize=(10,10))
    plt.style.use(astropy_mpl_style)
    ax = fig.gca()
    #plt.imshow(data, vmin=median-3*sigma, vmax=median+3*sigma,cmap='gray')
    plt.imshow(self,vmin=np.median(self)-1*np.std(self),
               vmax=np.median(self)+5*np.std(self),cmap='gray')
    circles = [plt.Circle(
            (coords[i][0], coords[i][1]),
            radius = 30,edgecolor='r',
            facecolor='None') for i in range(len(coords))
            ]
    for c in circles:
        ax.add_artist(c)

    plt.show()

#%%
def test_plot(image,coord,correc):
    imgshow(image)
    plt.plot(coord[0],coord[1],'ro')
    plt.plot(correc[0],correc[1],'go')
#%%
clickarr,n = [],1

img_click(get_img(files[0]))

c_cord = []

for c in clickarr:
    image = get_img(files[0])
    CRD,cent,cord = get_center(image,c,50)
    c_cord.append(cord)


print('clicked array',clickarr)
print('corrected coordinates',c_cord)

plot_src(get_img(files[0]),clickarr)
plot_src(get_img(files[0]),c_cord)
