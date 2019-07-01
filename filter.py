#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May  6 13:50:57 2019

@author: jishnu
"""
#%%
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from astropy import time, coordinates as coord, units as u
from astropy.io import fits
from astropy.time import Time
import glob
import os

#%%
# =============================================================================
# This script runs through files containing detected sources from Sextractor and Psfex
# source files have the following columns, refer back to photometry script 
# and sextractor to add/remove other columns
# Index(['X_IMAGE', 'Y_IMAGE', 'ALPHAWIN_J2000', 'DELTAWIN_J2000', 'SNR_WIN',
#        'MAG_POINTSOURCE', 'MAGERR_POINTSOURCE'],
#       dtype='object')
#
# The output is a dataframe and Lightcurve csv written to disk
# =============================================================================
#%%
'''
The coordinates for variable and check stars
'''
V_x,V_y = 1230,1750
C1x,C2x = 1670,288
C2x,C2y = 960,2054

#%%
#root_folder = '../../'
#os.chdir(root_folder)
os.chdir('/home/jishnu/JCBT/TYC3315/14Nov2/reduced/sources')
os.chdir('..')
ffiles = glob.glob('*proc.fits')
ffiles.sort()

#%%
def read(filename): #read a csv file and returns a dataframe
    df = pd.read_csv(filename,sep='\t')
    df = df.drop([df.columns[0]],axis=1)
    return df

#%%
def search(x,y,df,t):
    '''
    #here x,y are pixel coords of target star, t is the size of the box
    used for search
    '''
    flt = df[ #setting condtion so to consider sources inside a box around(x,y)
            (df['X_IMAGE'] < (x+t)) &
            (df['X_IMAGE'] > (x-t)) &
            (df['Y_IMAGE'] < (x+t)) &
            (df['Y_IMAGE'] > (x-t))
            ]
    print(flt)
    return flt['MAG_POINTSOURCE'],flt['MAGERR_POINTSOURCE'] #returns pandas series object 

#%%
def cone_search(x,y,df,dist): 
    '''
    #another implematation of previous function
    takes in (x,y) coordinates, a dataframe and search radius
    returns the magnitude and error if any source is 
    found around the given coordinates
    '''
    #uses distance formula to check how close the point is, to given threshold
    fdf = df[ 
            (np.sqrt((df['X_IMAGE']-x)**2 +(df['Y_IMAGE']-y)**2) < dist)
            ]
    print(fdf)
    
    return fdf['MAG_POINTSOURCE'],fdf['MAGERR_POINTSOURCE']

#%%
def coord_search(ra,dec,df,dist):
    '''
    Uses RA,DEC instead of pixel coordinates to search
    Use this function only if there is accurate WCS information and the 
    fits file has been platesolved
    '''
    fdf = df[ #uses distance formula to check how close the point is, to given threshold
            (np.sqrt((df['ALPHAWIN_J2000']-ra)**2 +
                     (df['DELTAWIN_J2000']-dec)**2) < dist)
            ]
    print(fdf)
    if len(fdf['MAG_POINTSOURCE'])==0:
        return 0,0 #returns 0 if no sources found
    else:
        return fdf['MAG_POINTSOURCE'],fdf['MAGERR_POINTSOURCE']
#%%
def sorting(x,y,files,t):
    '''
    Takes (x,y), files and distance threshold and returns the magnitude
    and error of the source located at (x,y)
    '''
    mag,err = [],[]
    for f in files: #honestly this could have been done in one function..
        print(f+'\n'*2)
        mags = cone_search(x,y,read(f),t) #pick your poison here
        mag.append(mags[0])
        err.append(mags[1])
        print('\n'*2+'-'*50)
    return mag,err #returns the magnitude and error columns, as pandas series


#%%
def filter(lst,err): #serves same purpose as the above function, now you might ask why there are two tools for the same job, to which i have no answer. They are nice, first one i came up with, in my head; and second one, well Internet people! https://stackoverflow.com/questions/12141150/from-list-of-integers-get-number-closest-to-a-given-value
    '''
    Input is pandas series of mag and error, if there are more than one elements
    go into the function, else just return the value from series
    Takes the magnitude and associated error columns, and removes duplicates
    (sometimes there are more than 1 source in the given (x,y) and d search region)
    The function tries to find the correct source by taking mean of the list
    and looking to reject other sources based on that
    
    Should not be used if the sources are clumped up and has the same mag
    
    Needs better filtering, until then; stuck with this dirty way
    '''
    opt,Err = [0],[]
    for i in range(len(lst)):
        try:
            if len(lst[i])==1:
                opt.append(lst[i].iloc[0])
                Err.append(err[i].iloc[0])
            else: #fancy pants way to do the earlier function
                opt.append(min([lst[i].iloc[n] for n in range(len(lst[i]))],
                               key=lambda x:abs(x-np.mean(opt))))
                Err.append(min([err[i].iloc[n] for n in range(len(err[i]))],
                               key=lambda x:abs(x-np.mean(opt))))
        except:
            opt.append(0)
            Err.append(0)
            
    return opt[1:],Err #returns 2 lists, mag and err. Do what you want with them


#%%
def Filter(lst,err):
    opt = [0] # creating list with one element to avoid np raising error 
              # while taking mean
    Err = []
    for i in range(len(lst)):
        if len(lst[i])==1:
            opt.append(lst[i].iloc[0])
            Err.append(err[i].iloc[0])
        else:
            if (abs(np.mean(opt)-lst[i].iloc[0])) > (abs(np.mean(opt)-lst[i].iloc[1])):
                opt.append(lst[i].iloc[1])
                Err.append(err[i].iloc[1])
            else:
                opt.append(lst[i].iloc[0])
                Err.append(err[i].iloc[0])
    
    return opt[1:],Err


#%%
def get_info(imageName):
    loc = coord.EarthLocation.of_site('vbo')
    hdul = fits.open(imageName)
    ut=(hdul[0].header['UT'])
    ra=(hdul[0].header['RA'])
    dec=(hdul[0].header['DEC'])
    am=(hdul[0].header['AIRMASS'])
    t = Time(hdul[0].header['DATE-OBS'], scale='utc')
    jd=(t.jd) #we only use JD.
    
    obs_target = coord.SkyCoord(hdul[0].header['RA'],hdul[0].header['DEC'],
                                unit = (u.hourangle,u.deg),frame = 'icrs'
                                )
    times = time.Time(t.jd, format='mjd',scale='utc', location=loc
                      )
    ltt_bary = times.light_travel_time(obs_target)
    time_barycentre = times.tdb + ltt_bary
    bjd=(time_barycentre.value)
    
    hdul.close()
    return ut,ra,dec,am,jd,bjd
    #USAGE
    #ut,ra,dec,am,jd,bjd = get_info(imagename)
#%%
df = pd.DataFrame()
ut,ra,dec,am,jd,bjd = [],[],[],[],[],[]
#%%
def main():
    for f in range(len(ffiles)):
        UT,RA,DEC,AM,JD,BJD = get_info(ffiles[f])
        ut.append(UT)
        ra.append(RA)
        dec.append(DEC)
        am.append(AM)
        jd.append(JD)
        bjd.append(BJD)
    
    folder = 'sources/'
    os.chdir(folder)

    Files = [i+'_source_cat.csv' for i in ffiles] #can just read the files off disk but where's the fun in that.
    m,e = sorting(V_x,V_y,Files,100) #loopy logic of this function! yay me
    Vmag,Verr = filter(m,e)
    
    m,e = sorting(C1x,C2x,Files,100)
    C1mag,C1err = filter(m,e)
    
    m,e = sorting(C2x,C2y,Files,100)
    C2mag,C2err = filter(m,e)

#    m,e = sorting(50.381181,47.433372,Files,0.01) #loopy logic of this function! yay me
#    Vmag,Verr = filter(m,e)
    
    df = pd.DataFrame()
    df = df.assign(
            **{'UT':ut,'RA':ra,'DEC':dec,'airmass':am,
               'JD':jd,'BJD':bjd,'Vmag':Vmag,'v_err':Verr,
               'C1':C1mag,'C1_err':C1err,'C2':C2mag,'C2_err':C2err
                    })
    df['Vmag_d'] = df['Vmag']-df['C2']
    df.to_csv('LC.csv')
    return df

#%%
plt.figure(figsize=(8,6))
plt.plot(bjd,df.Vmag,'o',label='v')
plt.plot(bjd,df.C1,'o',label='C1')
plt.plot(bjd,df.C2,'o',label='C2')
plt.legend()

#%%
plt.plot(bjd,df.Vmag-df.C1,'o',label='v-c1')
plt.plot(bjd,df.Vmag-df.C2,'o',label='v-c2')
plt.plot(bjd,df.C1-df.C2,'o',label='c1-c2')
plt.legend()