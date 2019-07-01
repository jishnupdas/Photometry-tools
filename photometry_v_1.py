#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec  8 14:52:20 2018

@author: jishnu

Modified on Fri Mar 1 15:32:00 2019
"""
#%%
import numpy as np
import numpy.ma as ma
import os
import time as zzz
from astropy.io import ascii
from astropy.coordinates import SkyCoord
from astropy import time, coordinates as coord, units as u
from astropy.wcs import WCS
from astropy.stats import sigma_clipped_stats, sigma_clip
from astropy.io import fits
from astropy.time import Time
from astropy.visualization import astropy_mpl_style
import matplotlib.pyplot as plt
import pandas as pd
import glob
import subprocess

#%%
# =============================================================================
# creating empty lists/containers to store values
# =============================================================================
obj,ut,ra,dec,airmass,obsdate = [],[],[],[],[],[]
exp,filtr = [],[]
jd,bjd = [],[] #for ut to HJD conversion
bias,flat = [],[]
masterbias,masterflat = 0,0
var,c1,c2,c3 = 0,0,0,0
instMagV,instMagVerr = [],[]
magV,magC1 = [],[]
magVerr,source_index = [],[]
magC1err,C1_index = [],[]
instMagC1,instMagC1err = [],[]
y,x = 0,0 #image dimensions
zero_psfmean, zero_psfmed, zero_psfstd = 0,0,0
'''Update the gain value below'''
g = 0.745 #gain read from the header
photoDistThresh = 50 #Pixel scale:0.29 arcsec/pixel
data = pd.DataFrame()
vbo = coord.EarthLocation.of_site('vbo')
#os.chdir('/home/jishnu/JCBT/11nov2018/Photometry')
#%%

Objects = ['Object','Comparison']
RA  = ['Object RA','Comp RA']
DEC = ['Obj RA','Comp Dec']


#%%
def object(idx):

    coordinates = coord.SkyCoord([RA[idx]],[DEC[idx]],
                            unit = (u.hourangle,u.deg),frame = 'icrs'
                            )
    return coordinates


#%%
Source = object(0)
C1 = object(1)
Source_name = 'TYC3315'
Source_name.strip(',')
Source_name.replace(' ','_')

#%%
os.chdir('../reduced')
f = os.getcwd()
try:
    os.system('mkdir sources')
    os.system('mkdir plots')
    os.system('mkdir frames')
    os.system('mkdir radialprofiles')
    os.system('mkdir PSF')
except:
    pass

os.chdir('../Photometry/config')
os.system('cp * '+f)
os.chdir(f)

#%%
files = glob.glob('*proc.fits') #getting all the fits files in the directory
files.sort()

#%%

'''code block to keep ut1 time updated'''
from astroplan import download_IERS_A
# ### Using UT1
# - To keep accurate time, the changes in earth's rotation period have to be
# - taken into account.
# - AstroPy does this using a convention called UT1, that is tied to the
# - rotation of earth with respect to the position of distant quasars.
# - IERS - International Earth Rotation and Reference Systems Service keeps
# - continuous tabs on the orientation of the earth and updates the data in
# - the IERS bulletin.
# Update the bulletin:

from astropy.utils import iers
iers.IERS_A_URL = 'http://maia.usno.navy.mil/ser7/finals2000A.all' #default
#iers.IERS_A_URL = 'http://toshi.nofs.navy.mil/ser7/finals2000A.all'
#iers.IERS_A_URL = 'https://datacenter.iers.org/data/latestVersion/19_EOP_C01.1846-1899_V2013_0119.txt'
#default server throws up, HTTPError :403, so we are using alternate mirror
download_IERS_A()
#astroplan.get_IERS_A_or_workaround()
#%%
def get_img(self): #smart way #returns 2D array from fits file data
    hdul = fits.open(self)
    image = np.asarray(hdul[0].data)
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
def get_img_dims(self): #returns an array and x,y shappe of image
    hdul = fits.open(self)
    image = np.asarray(hdul[0].data)
    y,x = image.shape[-2],image.shape[-1] #getting y and x axes
    hdul.close() #closing the file. GOTTA save RAM!
    return image,y,x


#%%
def progress(i):
        print('''
    	File:%s ---- %.2f%s '''%(i,
                                      (files.index(i)/len(files)*100),
                                      '%'))
#%% 
astrometry = '../Photometry/'+'autoastrometry.py'

def WCS_astrometry():
    print('Running astrometry...')
    for i in range(len(files)):
        try:
            subprocess.run(
                    ['python',astrometry,files[i],
                     '-px','0.2898','-c','tmc'],
                     check=True)
            progress(i)
        except subprocess.CalledProcessError as err:
            print('''Could not run autoastrometry with error %s.
                  Check if file exists.'''%err)
            
#%%
WCS_astrometry()

#%%
def imgshow(self): #display the image
    plt.figure(figsize=(8,8)) #size
    plt.style.use(astropy_mpl_style) #styling and gimmicks of astropy
    
    mean, median, sigma = sigma_clipped_stats(self) #use this for contrast 
                                                    #stretching
    plt.imshow(self,cmap='gray',vmin=median-6*sigma, vmax=median+6*sigma)
    plt.colorbar() #a bar showing the scale of the image
    #plt.show()
    plt.close() #clean up after yourself!!

#%%
def astrometry(hdul): #single use
    global raImage,decImage,W,Image
    #Image = get_img_slow(hdul)
    Image = hdul[0].data
    
    '''
    obs_target = coord.SkyCoord(hdul[0].header['RA'],hdul[0].header['DEC'],
                                unit = (u.hourangle,u.deg),frame = 'icrs'
                                )
    '''
    #[ra, dec] = np.array(obs_target.ra.degree),np.array(obs_target.dec.degree)
    #Get the RA and Dec of the center of the image
    [raImage, decImage] =W.all_pix2world(
                                        Image.shape[0]/2,
                                        Image.shape[1]/2,1
                                        )

    return raImage,decImage,[ra,dec]

#%%
def Query(raImage, decImage,i,W):    #single use
    #Set the box size to search for catalog stars
    boxsize = 25 # arcminutes
    #Magnitude cut-offs of sources to be cross-matched against
    maxmag = 16
    
    from astroquery.vizier import Vizier
    #Vizier.VIZIER_SERVER = 'vizier.ast.cam.ac.uk'
    catNum = 'II/349'#This is the catalog number of PS1 in Vizier
    print('''\nQuerying Vizier %s around RA %.4f,
          Dec %.4f with a radius of %.4f arcmin'''
          %(catNum, raImage, decImage, boxsize))

    try:
        #You can set the filters for the individual columns 
        #(magnitude range, number of detections) inside the Vizier query
        v = Vizier(columns=['*'],
                   column_filters={"gmag":"<%.2f"%maxmag,
                                   "Nd":">6", "e_gmag":"<1.086/3"},
                                   row_limit=-1)
        Q = v.query_region(
                    SkyCoord(ra = raImage, dec = decImage,
                             unit = (u.deg, u.deg)),
                             radius = str(boxsize)+'m',
                             catalog=catNum, cache=False
                        )
        #query vizier around (ra, dec) with a radius of boxsize
        #ascii.write(Q[0], i+'cat_sources.csv', format='csv',overwrite=True)
        #print(Q[0])
        
    except:
        print('cannnot reach the Vizier database. check your internet!')
        
    global ps1_imCoords
    #Convert the world coordinates of these stars to image coordinates
    ps1_imCoords = W.all_world2pix(Q[0]['RAJ2000'], Q[0]['DEJ2000'], 1)
    
    #Another round of filtering where we reject sources close to the edges
    good_cat_stars = Q[0][np.where((ps1_imCoords[0] > 100) & 
                      (ps1_imCoords[0] < 1900) & 
                      (ps1_imCoords[1] > 100) & 
                      (ps1_imCoords[1] < 3900))]
    print(good_cat_stars)
    #taking only good catalogue stars
    ps1_imCoords = W.all_world2pix(good_cat_stars['RAJ2000'],
                                   good_cat_stars['DEJ2000'], 1)
    return Q[0],ps1_imCoords,good_cat_stars
    #plot_sources(Image,ps1_imCoords)


#%%

reference = files[1]
refimage = get_img(reference) #file to be used as reference frame
refimage += 100
refimage = np.log10(refimage)
hdul = fits.open(reference)
raImage,decImage = 0,0
Image = 0
ps1_imCoords = []
W = WCS(hdul[0].header)
#%%

astrometry(hdul)
Q,ps1_imCoords,ref_imCoords = Query(raImage, decImage,0,W)
ref_imCoords.write('sources/'+'reference_sources_cat.csv',overwrite=True)
#%%
''' 
    this gives the coords of sources in the 
    reference catalogue,
    this will be used while cross matching sources, and 
    computing ZP for the frames
'''
ps1CatCoords = SkyCoord(
        ra=ref_imCoords['RAJ2000'], 
        dec=ref_imCoords['DEJ2000'], 
        frame='icrs', unit='degree'
        )

#%%
def plot_sources(data,ps1_imCoords):
    '''plot the sources on the image'''
    mean, median, sigma = sigma_clipped_stats(data)
    fig = plt.figure(figsize=(10,10))
    plt.style.use(astropy_mpl_style)
    ax = fig.gca()
    plt.imshow(data, vmin=median-3*sigma, vmax=median+3*sigma,cmap='gray')
    circles = [plt.Circle(
            (ps1_imCoords[0][i], ps1_imCoords[1][i]),
            radius = 30,edgecolor='r',
            facecolor='None') for i in range(len(ps1_imCoords[0]))
            ]
    for c in circles:
        ax.add_artist(c)
        
    #plt.show()
    plt.close()


#%%
def source_detection(imageName):
    
    c = os.getcwd()
    
    os.chdir('../reduced/')
    f = os.getcwd()
    os.chdir('../Photometry/config')
    os.system('cp photomCat.sex '+f)
    os.system('cp photomCat.param '+f)
    os.chdir(f)
    
    
    configFile = 'photomCat.sex'
    paramName  = 'photomCat.param'
    
    open(imageName+'.cat','w+')
    catalogName = imageName+'.cat'
    
    try:
        subprocess.run(
                ['sextractor', '-c', configFile, 
                 imageName,
                 '-CATALOG_NAME', catalogName, 
                 '-PARAMETERS_NAME', paramName],
                 check=True
                 )
    except subprocess.CalledProcessError as err:
        print('Could not run sextractor with exit error %s'%err)
    
    os.chdir(c)  

#%%
def get_table_from_ldac(filename, frame=1):
    """
    Load an astropy table from a fits_ldac by frame (Since the ldac format has 
    column info for odd tables, giving it twce as many tables as a regular 
    fits BinTableHDU,match the frame of a table to its corresponding frame in
    the ldac file).
    
    Parameters
    ----------
    filename: str
        Name of the file to open
    frame: int
        Number of the frame in a regular fits file
    """
    from astropy.table import Table
    if frame>0:
        frame = frame*2
    tbl = Table.read(filename, hdu=frame)
    return tbl
#%%
def LDAC_cat(catalogName,filename):
    #This is a python wrapper for reading LDAC files produced by SExtractor
    sourceTable = get_table_from_ldac(catalogName)
    #Let's look at the contents of the table
    #print(sourceTable.colnames)
    #print(sourceTable)
    '''
    ascii.write(
            sourceTable[0], 
            filename+'_extracted_sources.csv', 
            format='csv',
            overwrite=True)
    '''
    #sourceTable.write(filename+'_extracted_sources.csv', format='csv')
    return sourceTable

#%%
def clean_sources(imageName):
    #filter on the sources to select the ones satisfying our criteria
    catalogName = imageName+'.psf.cat'
    sourceTable = get_table_from_ldac(catalogName)
    cleanSources = sourceTable[
            (sourceTable['FLAGS']==0)          &
            (sourceTable['FWHM_WORLD'] < 2   ) & 
            (sourceTable['XWIN_IMAGE'] < 1800) & 
            (sourceTable['XWIN_IMAGE'] > 200 ) &
            (sourceTable['YWIN_IMAGE'] < 3800) &
            (sourceTable['YWIN_IMAGE'] > 200 ) &
            (sourceTable['SNR_WIN']    > 150 )
            ]

    return cleanSources


#%%
def detected_sources_plot(imageName):
    data = get_img(imageName)
    cleanSources = clean_sources(imageName)
    fig = plt.figure(figsize=(10,10))
    plt.style.use(astropy_mpl_style)
    mean, median, sigma = sigma_clipped_stats(data)
    ax = fig.gca()
    plt.imshow(data, vmin=median-6*sigma, vmax=median+8*sigma,cmap='gray')
    #plotting circles on top of all detected sources
    circles = [plt.Circle(
            (source['XWIN_IMAGE'], 
             source['YWIN_IMAGE']), 
             radius = 30, edgecolor='r',
             facecolor='None') for source in cleanSources]
    for c in circles:
        ax.add_artist(c)
    
    plt.savefig('frames/'+imageName+'field.png',dpi=150)
    #plt.show()
    plt.close()
    

#%%    
def psf_photometry(imageName):
    
    psfConfigFile = 'psfex_conf.psfex'
    catalogName = imageName+'.cat'
    
    try:
        subprocess.run(
                ['psfex',
                 '-c', psfConfigFile,
                 catalogName], 
                 check=True
                 )
    except subprocess.CalledProcessError as err:
        print('Could not run psfex with exit error %s'%err)

#%%
def Curve_of_growth(imageName):    
    psfModelHDU = fits.open('moffat_'+imageName+'.fits')[0]
    psfModelData = psfModelHDU.data
    
    mean, median, std = sigma_clipped_stats(psfModelData)
    plt.figure(figsize=(6,6))
    plt.title(imageName)
    plt.style.use(astropy_mpl_style)
    plt.imshow(psfModelData, vmin=0, vmax=median+20*std)
    plt.savefig('PSF/'+imageName+'_psf.png',dpi=150)
    #plt.show()
    plt.close()
    

    psfImageCenter = [(psfModelData.shape[0]-1)/2, (psfModelData.shape[1]-1)/2]
    y, x = np.indices(psfModelData.shape)
    r = np.sqrt((x-psfImageCenter[0])**2 + (y-psfImageCenter[1])**2)
    r = r.astype(np.int)
    
    tbin = np.bincount(r.ravel(), psfModelData.ravel())
    nr = np.bincount(r.ravel())
    radialprofile = tbin/nr
    
    plt.figure(figsize=(6,6))
    plt.title(imageName)
    plt.style.use(astropy_mpl_style)    
    plt.plot(range(len(radialprofile)), radialprofile, 'k.', markersize=15)
    plt.xlabel('Radial distance (pixels)', fontsize=15)
    plt.ylabel('PSF Amplitude', fontsize=15)
    plt.xlim(-1,20)
    plt.savefig('radialprofiles/'+imageName+'_COG.png',dpi=150)
    #plt.show()
    plt.close()


#%%
def psf_fitting(imageName):
    c = os.getcwd()
    
    os.chdir('../reduced/')
    f = os.getcwd()
    os.chdir('../Photometry/config')
    os.system('cp photomCat.sex '+f)
    os.system('cp photomPSF.param '+f)
    os.chdir(f)
    
    configFile = 'photomCat.sex'
    psfName = imageName + '.psf'
    open(imageName+'psf.cat','w+')
    psfcatalogName = imageName+'.psf.cat'
    psfparamName = 'photomPSF.param' 
    
    '''
    This is a new set of parameters to be obtained from SExtractor,
    including PSF-fit magnitudes
    '''
    
    try:
        '''
        We are supplying SExtactor with the PSF model with the
        PSF_NAME option
        '''
        subprocess.run(
                ['sextractor',
                 '-c', configFile,
                 imageName,
                 '-CATALOG_NAME', psfcatalogName,
                 '-PSF_NAME', psfName,
                 '-PARAMETERS_NAME',
                 psfparamName],
                 check=True
                 )
    
    except subprocess.CalledProcessError as err:
        print('Could not run sextractor with exit error %s'%err)
    
    os.chdir(c)

#%%
def psf_catalogue(imageName):
    PSFSources = clean_sources(imageName)
    cat = pd.DataFrame()
    cat = cat.assign(#adding values to dataframe
            **{'X_IMAGE':PSFSources['X_IMAGE'],
               'Y_IMAGE':PSFSources['Y_IMAGE'],
               'ALPHAWIN_J2000':PSFSources['ALPHAWIN_J2000'],
               'DELTAWIN_J2000':PSFSources['DELTAWIN_J2000'],
               'SNR_WIN':PSFSources['SNR_WIN'],
               'MAG_POINTSOURCE':PSFSources['MAG_POINTSOURCE'],
               'MAGERR_POINTSOURCE':PSFSources['MAGERR_POINTSOURCE'],
               }
            )
    print('writing %s sources to file'%imageName)
    #cat.dropna()
    cat.to_csv('sources/'+imageName+'_source_cat.csv',sep='\t')

#%%
'''
    WORK to do HERE!!!
    this can be modularized, and broken. Current implementation has
    repetetive blocks for the target and standard star.
    If this can be made into a function, it can take the image name and
    targets(V,C1,C2...) as arguments, at the same time record() function
    can be modified to take arguements.

'''

def Corrected_Magnitude(imageName,zero_psfmed, zero_psfstd):
    cleanPSFSources = clean_sources(imageName)
    sourceCatCoords = SkyCoord(
            ra=cleanPSFSources['ALPHAWIN_J2000'],
            dec=cleanPSFSources['DELTAWIN_J2000'],
            frame='icrs', unit='degree'
            )
    idx_Source, idx_cleanpsf_Source, d2d, d3d = sourceCatCoords.search_around_sky(
            Source, photoDistThresh*u.arcsec)
    #print('Found the source at index %d'%idx_cleanpsf_Source[0])
    
    try:
        psfinstmag = cleanPSFSources[idx_cleanpsf_Source]['MAG_POINTSOURCE'][0]
        psfinstmagerr = cleanPSFSources[idx_cleanpsf_Source]['MAGERR_POINTSOURCE'][0]
        
        psfmag = zero_psfmed + psfinstmag
        psfmagerr = np.sqrt(psfinstmagerr**2 + zero_psfstd**2)
        
        instMagV.append(psfinstmag)
        instMagVerr.append(psfinstmagerr)
        magV.append(psfmag)
        magVerr.append(psfmagerr)            
        source_index.append(idx_cleanpsf_Source)
                
        print('''
              ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
              | PSF-fit magnitude of SOURCE is %.2f +/- %.2f |
              | at index %d                                    |
              |_______________________________________________|
              '''%(psfmag, psfmagerr,idx_cleanpsf_Source[0]),
              idx_cleanpsf_Source)
    except:
        instMagV.append(0)
        instMagVerr.append(0)
        magV.append(0)
        magVerr.append(0)
        source_index.append(-1)
        print('Trouble finding the source!')
        
#%%    
#   This bit is to find psf mag and W mag of comparison
    idx_Source, idx_cleanpsf_Source, d2d, d3d = sourceCatCoords.search_around_sky(
        C1, photoDistThresh*u.arcsec)
    #print('Found the source at index %d'%idx_cleanpsf_Source[0])
    
    try:
        psfinstmag = cleanPSFSources[idx_cleanpsf_Source]['MAG_POINTSOURCE'][0]
        psfinstmagerr = cleanPSFSources[idx_cleanpsf_Source]['MAGERR_POINTSOURCE'][0]
        
        psfmag = zero_psfmed + psfinstmag
        psfmagerr = np.sqrt(psfinstmagerr**2 + zero_psfstd**2)
        
        instMagC1.append(psfinstmag)
        instMagC1err.append(psfinstmagerr)
        magC1.append(psfmag)
        magC1err.append(psfmagerr)            
        C1_index.append(idx_cleanpsf_Source)
        
    except:
        instMagC1.append(0)
        instMagC1err.append(0)
        magC1.append(0)
        magC1err.append(0)
        C1_index.append(-1)
        print('Trouble finding the C1!')

#%%
def psf_cross_matching(imageName):
    cleanPSFSources = clean_sources(imageName)
    
    psfCatCoords = SkyCoord(
            ra=cleanPSFSources['ALPHAWIN_J2000'], 
            dec=cleanPSFSources['DELTAWIN_J2000'], 
            frame='icrs', unit='degree'
            )
    '''
     Now cross match sources
     cross-match distance threshold was set in the beginning, 
     update to your liking 
     Currently set to a generous value, reduce for crowded fields
     or precise pointing and good seeing. (unrealistic at VBO)
    '''
    
    idx_psfimage, idx_psfps1, d2d, d3d = ps1CatCoords.search_around_sky(
            psfCatCoords, photoDistThresh*u.arcsec
            )
    
    print('Found %d good cross-matches'%len(idx_psfimage))
    
    psfoffsets = ma.array(
            ref_imCoords['gmag'][idx_psfps1] - 
            cleanPSFSources['MAG_POINTSOURCE'][idx_psfimage])
    #Compute sigma clipped statistics
    global zero_psfmean, zero_psfmed, zero_psfstd
    zero_psfmean, zero_psfmed, zero_psfstd = sigma_clipped_stats(psfoffsets)
    
    Cross_match_plot(cleanPSFSources,imageName,idx_psfps1,idx_psfimage)
    
    Corrected_Magnitude(imageName,zero_psfmed, zero_psfstd)
    #return zero_psfmean, zero_psfmed, zero_psfstd

#%%
def Cross_match_plot(cleanPSFSources,imageName,idx_psfps1,idx_psfimage):
    plt.figure(figsize=(8,8))
    plt.title(imageName)
    plt.style.use(astropy_mpl_style)
    plt.plot(
            ref_imCoords['gmag'][idx_psfps1], 
            cleanPSFSources['MAG_POINTSOURCE'][idx_psfimage],
            'r.', markersize=14,
            )
    plt.text(-0.03,0.82,
         '''
          PSF Mean ZP    : %.2f
          PSF Median ZP : %.2f
          PSF STD ZP      : %.2f
         '''%(zero_psfmean, zero_psfmed, zero_psfstd),
         transform=plt.gca().transAxes)
    #plt.ylim(-16, -7.5)
    plt.xlabel('PS1 magnitude', fontsize=15)
    plt.xlim(11,16)
    plt.ylabel('Instr PSF-fit mag', fontsize=15)
    plt.ylim(-15,-10)
    plt.savefig('plots/'+imageName+'crossmatch.png',dpi=150)
    #plt.show()
    plt.close()

#%%
def record(Source_name):
    os.chdir('../')
    global data
    data = data.assign(#adding values to dataframe
            **{'UT':ut,'RA':ra,'Dec':dec,'Filter':filtr,
               'File':files,'Idx':source_index,'Exposure':exp,
               'Airmass':airmass,'Date':obsdate,'JD':jd,
               'mag V':magV, 'magV Err':magVerr,
               'instMagV':instMagV,'instMagVErr':instMagVerr
               }
            )
    data = data.assign(
            **{'mag C1':magC1, 'magC1 Err':magC1err,
               'instMagC1':instMagC1,'instMagErr':instMagC1err
               }
            )
    data = data.assign(#performing extinction correction
            **{'EC instmagV':data['instMagV']-(data['Airmass']*0.2),
               'EC instmagC1':data['instMagC1']-(data['Airmass']*0.2),
               }
            )
    data = data.assign(
            **{'BJD':bjd,
               'Diff mag V':data['EC instmagC1']-data['EC instmagV'],
               'Diff mag Err':instMagVerr
               }
            )
    #print(data)
    data.dropna()
    data.to_csv('LC_'+str(Source_name)+'_'+obsdate[0]+'_.csv',sep=',')
    
#%%
'''extinction m(X)  =  m0   +  k*X
   corrected m0  =  m(X)  -  k*X 
   
   passband         k
 ----------------------
      U            0.6
      B            0.4
      V            0.2
      R            0.1
      I            0.08
      
-------------------------------------------------------
formula; mag = K-2.5log(flux), K is zeropoint magnitude
-------------------------------------------------------

'''
#%%
'''
 Things to do
 1.Copy WCS into all files #DONE!!!
     (http://www.astropy.org/astropy-tutorials/FITS-header.html)
     a.New implementation, using Dan Pereley's astrometry code
 
 2.write the tableinto file. -                          DONE!
 3.loop through all the images. -                       DONE!
 4.psf fitting using psfex                              DONE!
 5.cross matching, ZP deviation for each file           DONE!
 6.read mag from all the files, apply zp correction     DONE!
 7.write to lightcurve file as BJD,Mag & mag error      DONE!
 8.perform differential photometry on the frames        Yet to do.
'''

#%%
def get_info(imageName):
    hdul = fits.open(imageName)
    ut.append(hdul[0].header['UT'])
    obj.append(hdul[0].header['OBJECT'])
    ra.append(hdul[0].header['RA'])
    dec.append(hdul[0].header['DEC'])
    exp.append(hdul[0].header['EXPTIME'])
    airmass.append(hdul[0].header['AIRMASS'])
    filtr.append(hdul[0].header['FILTER'])
    obsdate.append(hdul[0].header['DATE-OBS'])
    t = Time(hdul[0].header['DATE-OBS'], format='fits', scale='utc')
    jd.append(t.jd) #we only use JD.
    
    obs_target = coord.SkyCoord(hdul[0].header['RA'],hdul[0].header['DEC'],
                                unit = (u.hourangle,u.deg),frame = 'icrs'
                                )
    times = time.Time(t.jd, format='mjd',scale='utc', location=vbo
                      )
    ltt_bary = times.light_travel_time(obs_target)
    time_barycentre = times.tdb + ltt_bary
    bjd.append(time_barycentre.value)
    
    hdul.close()
    
#%%
def main(): 
    os.chdir('../reduced/') #cd into the directory where reduced files reside.
    #clean slate boi!
    #astrometry(hdul)
    for i in files:
        progress(i)
                #life is a progress bar, from dust to dust....
        
        source_detection(i)
        #zzz.sleep(2)
        #detected_sources_plot(i)
        psf_photometry(i)
        Curve_of_growth(i)
        psf_fitting(i)
        psf_catalogue(i)
        psf_cross_matching(i)
        get_info(i)
        
#%%
def imageshow(self): #display the image
    plt.figure(figsize=(8,8)) #size
    plt.style.use(astropy_mpl_style)
    plt.imshow(self,cmap='gray')
    plt.colorbar()
    #plt.show()
    plt.close()

#%%
def LC_plot():
    plt.figure(figsize=(10,6))
    plt.scatter(bjd,magV,label='cr_mat Mag')
    plt.scatter(bjd,[zero_psfmed+i+.5 for i in data['EC instmagV']],
                label='inst mag')
    plt.scatter(bjd,[zero_psfmed+i+.5 for i in data['EC instmagC1']],
                label='inst mag C1')
    plt.legend()
    #plt.ylim(11,max(magV)+0.01)
    plt.gca().invert_yaxis()
    plt.savefig('LC.png',dpi=150)
    plt.show()
    plt.close()
    
    plt.figure(figsize=(10,6))
    plt.scatter(bjd,data['Diff mag V'],label='Differential Mag')
    #plt.ylim(-0.5,0.5)
    plt.legend()
    plt.savefig('LC_Diff_mag.png',dpi=150)
    plt.show()
    plt.close()

#%%
main()
record(Source_name)
LC_plot()
