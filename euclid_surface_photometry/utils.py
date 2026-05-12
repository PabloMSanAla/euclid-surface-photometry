# utils.py: Utility functions for the Euclid surface photometry pipeline

# (Add shared utility functions here, e.g., file I/O, config loading, logging helpers)

from .morphology import *

import os
import yaml
import json
import pickle
import copy
import collections.abc

from os.path import join, basename, abspath, dirname
import numpy as np
import pandas as pd
import cv2
from scipy.optimize import curve_fit, minimize
from scipy.special import gammaincinv
from scipy.signal import savgol_filter, medfilt
from scipy.signal import find_peaks, peak_prominences, peak_widths

from skimage.segmentation import  watershed

import astropy.units as u
from astropy.io import fits 
from astropy.table import Table, vstack
from astropy.stats import sigma_clipped_stats
from astropy.coordinates import SkyCoord
from astropy.visualization import PercentileInterval, LogStretch, make_lupton_rgb

from photutils.aperture import EllipticalAperture

from collections.abc import Iterable

from matplotlib import pyplot as plt
from matplotlib.colors import colorConverter, ListedColormap
from matplotlib.patches import Rectangle
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
from matplotlib.gridspec import GridSpec

from fabada import fabada

from astropipe.classes import Image, AstroGNU, MTObjects
from astropipe.profile import Profile, background_estimation_euclid, elliptical_radial_profile, autoprof_isophote_photometry, find_mode, measureImageNoise
from astropipe.utils import find, cutout, get_pixel_scale, derivative, mag_limit, morphology, center_hillclimb, redshift_to_kpc, arcsec_to_kpc
from astropipe.utils import closest
from astropipe.plotting import make_cmap, show

import warnings 
warnings.filterwarnings("ignore")


plt.rcParams['text.usetex']=True

## Global parameters of the pipeline
__SHAPE__ = (19200,19200)


def savefits(image, header, filename, keywords=None):
    ''' Save fits file adding new keywords in header'''

    header['COMMENT'] = 'astropipe Keywords'
    if not keywords is None:
        for key in keywords:
            header[key[0]]=(key[1],key[2])
    
    fits.PrimaryHDU(image, header).writeto(filename, overwrite=True)
    return os.path.isfile(filename)
    

def euclid_mask(outPath, center, out, fwhm=5):

    filename = find(outPath, f'*VIS*cutout.fits')[0]
    __keywords__ = [('INPUT',filename, 'Input filename'),
                    ('PIPE', 2, 'Euclid astropipe pipeline number'),
                    ('AP_TYPE', 'MASK', 'Type of product from astropipe pipeline'),
                    ('METHOD','MTO+NC', 'Methods to create product')]

    files = find(outPath, f'EUC_SEGMID*cutout.fits')
    segmFile = find(outPath, f'*segmap.fits')[0]
    __keywords__.append(('NFILES',len(files), 'Number of colours used to create the mask'))
    object_id = basename(filename).split('_')[2]

    # Initialize the image and segmentation map
    segmap = fits.getdata(segmFile)
    image = Image(filename)
    ra,dec = center
    image.obj(ra,dec)
    name = eucl_name_formater(ra,dec)
    
    # Reduce size of segment ID to improve the speed of the algorithm
    segmap = segmap - np.nanmin(segmap[segmap>0]) + 1
    segmap[segmap<0] = 0
    segmap = segmap.astype(np.uint16)
    maskObj = (segmap==segmap[image.y,image.x])

    # Measure the cosine similarity of the colours of the objects
    eucFile = [f for f in files if ('VIS' in f) or ('NIR')]
    cosine = cosine_similarity_image(eucFile, maskObj)

    ## Smoothing
    # _,_,std = sigma_clipped_stats(image.data)
    # smooth = fabada(image.data, std**2)
    
    # Run MTObjects for deblending
    mto = MTObjects()
    mto.move_factor = 0.3
    mto.run(image.data)
     
    # Run noisechisel for detection and dilate
    nc = AstroGNU(filename,loc='')
    nc.noisechisel(config='-N 1', keep=True)
    nc.segment()
    os.remove(nc.nc_file)
    nc.detections = cv2.dilate(nc.detections, np.ones((3,3),np.uint8))

    # Watershed segmentation map into noisechisel detection using colour image
    # markers = watershed(1/cosine, markers=segmap, mask=nc.detections)  # skimage watershed
    rgb = rgb_euclid(object_id, mask=nc.detections, abspath=os.path.dirname(filename))
    watermask = np.zeros_like(segmap) + segmap + 1
    watermask[(watermask==1)*(nc.detections==1)] = 0
    markers = cv2.watershed(rgb.astype(np.uint8), watermask.astype(np.int32))
    markers[markers==1] = 0           # Transform background value of cv2 and segment line to background 0
    markers[markers==-1] = 0          # Transform boundaries of watershed into background 0
    markers = markers.astype(np.uint16)

    # mask anything detected by noisechisel and not masked in watershed process
    nc_components = cv2.connectedComponentsWithAlgorithm(nc.detections.astype(np.uint8), connectivity=8, ltype=cv2.CV_32S, ccltype=cv2.CCL_WU)[1]
    nc_components[nc_components==nc_components[image.y,image.x]]=0
    markers[(nc_components>0)*(markers<1)] = np.nanmax(markers) + 1
    
    # Plotting
    plt.rcParams['axes.titlesize'] = 24
    qmin,qmax = np.nanpercentile(cosine[maskObj],[2,99])    
    fig, axes = plt.subplots(2,4,figsize=(20,10), sharex=True, sharey=True)
    ((ax1, ax2, ax3, ax7), (ax4, ax5, ax6,ax8)) = axes
    ax1.imshow(rgb, origin='lower')
    ax1.set_title(name)
    im = ax2.imshow(cosine, origin='lower',vmin=qmin,vmax=qmax)
    cax = inset_axes(ax2, width="100%", height="100%", loc='upper center',
                 bbox_to_anchor=(0.1, 0.91, 0.8, 0.05), bbox_transform=ax2.transAxes)
    cbar = fig.colorbar(im, cax=cax, orientation='horizontal')
    cbar.ax.tick_params(color='white', labelcolor='white',labelsize=14)
    ax2.set_title('Cosine Similarity')
    ax3.imshow(segmap, origin='lower', cmap=make_cmap(np.nanmax(markers)), interpolation='none')
    ax3.set_title(r'\textit{EUCLID} SEGMAP')
    ax4.imshow(nc.objects, origin='lower',cmap=make_cmap(np.nanmax(nc.objects)), interpolation='none')
    ax4.set_title(r'\texttt{NoiseChisel} + \texttt{Segment}')
    ax7.imshow(markers, origin='lower', cmap=make_cmap(np.nanmax(markers)), interpolation='none')
    ax7.set_title(r'SEGMAP + \texttt{NoiseChisel} + Watershed',fontsize=21)

    # Deblend sources of mtobjects according to their cosine similarity
    _,med_sim,std_sim = sigma_clipped_stats(cosine[(maskObj)*(mto.objects==mto.objects[image.y,image.x])])
    insiders = np.unique(mto.objects[(maskObj) * (mto.objects!=mto.objects[image.y,image.x]) * (mto.objects!=0)])
    galArea = np.sum(markers==markers[image.y,image.x])    
    for i in insiders:
        if i==0: continue
        area = np.sum(mto.objects==i)
        # Must be smaller than the main object [avoid larger objects that surrounds the small one]
        if area<galArea: 
            # If colour is compatible, then join to same ID as object, if not, add new ID
            if (np.nanmean(cosine[mto.objects==i]) < med_sim - 3*std_sim) or (area<fwhm**2):
                markers[mto.objects==i] = markers[image.y,image.x]
            else:
                markers[mto.objects==i] = np.nanmax(markers) + 1
    
    # Outliers:  Find the outliers in the galaxy region
    # Define region used to find outliers
    angle, sma, eps = morphology((segmap==segmap[image.y,image.x]).astype(np.uint8))
    aperture = EllipticalAperture((image.x,image.y), 1.2*sma, 1.2*sma*(1-eps), angle)
    mask = aperture.to_mask(method='center').to_image(segmap.shape)
    aperture = EllipticalAperture((image.x,image.y), 0.15*sma, 0.15*sma*(1-eps), angle) # Protect the center of the galaxy
    inner = aperture.to_mask(method='center').to_image(segmap.shape)
    mask = mask - inner

    # Avoid masking the object and leaving unmask the center
    inner_objects = np.unique(markers[(inner>0)*(maskObj)*(markers!=0)])
    if any(inner_objects): markers[np.isin(markers,inner_objects)] = markers[image.y,image.x]

    # Find outliers 
    outliers = (cosine < med_sim - 3*std_sim)*(cosine > med_sim - 15*std_sim)*mask
    kernel = np.ones((fwhm,fwhm))
    outliers = cv2.dilate(cv2.erode(outliers.astype(np.float32),kernel),kernel).astype(bool)
    outliers = cv2.connectedComponentsWithAlgorithm(outliers.astype(np.uint8), connectivity=8, ltype=cv2.CV_32S, ccltype=cv2.CCL_WU)[1]
    outliers = cv2.dilate(outliers.astype(np.float32), kernel).astype(np.int32)
    # markers[outliers.astype(bool)] = 0
    resetlab = np.nanmax(markers) + 1
    markers[outliers.astype(bool)] = (resetlab + outliers)[outliers.astype(bool)]
    # markers = markers.astype(np.int32) + outliers + resetlab
    # markers[markers==resetlab] = 0

    ax5.imshow(mto.objects, origin='lower', cmap=make_cmap(np.nanmax(mto.objects)), interpolation='none')
    ax5.set_title(r'\texttt{MTObjects}')
    ax6.imshow(outliers, origin='lower', cmap=make_cmap(np.nanmax(outliers)), interpolation='none')
    ax6.set_title('Outliers')
    ax8.imshow(markers, origin='lower', cmap=make_cmap(np.nanmax(markers)), interpolation='none')
    ax8.set_title('Our Method')
    # for ax in axes.flatten():
    #     ax.set_axis_off()
    fig.tight_layout()
    fig.savefig(out.replace('.fits','_segment.pdf'),format='pdf')
    fig.savefig(out.replace('.fits','_segment.jpg'),dpi=300)


    # Save result
    done = savefits(markers,image.header,out,keywords=__keywords__)
    plt.close('all')
    
    # Plotting of final mask with rgb image
    markers[markers==markers[image.y,image.x]]=0
    fig, ax = plt.subplots(figsize=(10,10))
    ax.imshow(rgb, origin='lower',interpolation='gaussian')
    ax.imshow(markers.astype(bool), origin='lower',cmap=mask_cmap(alpha=0.5))
    ax.plot(image.x, image.y, 'wx', ms=5)
    fig.tight_layout()
    fig.savefig(out.replace('.fits','.jpg'),dpi=300, bbox_inches='tight', pad_inches=0)
    fig.savefig(out.replace('.fits','.pdf'),format='pdf', bbox_inches='tight', pad_inches=0)

    plt.close('all')

    return done

def cosine_similarity_image(files, mask):
    
    shape = fits.getdata(files[0]).shape
    v_obj = np.zeros(len(files))
    colour = np.zeros((shape[0],shape[1],len(files)))

    for i in range(len(files)):
        colour[:,:,i] = fits.getdata(files[i])
        _,v_obj[i],_ = sigma_clipped_stats(colour[:,:,i][mask])

    v_obj /= np.linalg.norm(v_obj)
    norm = np.linalg.norm(colour, axis=2, keepdims=True)
    norm[norm == 0] = 1  # Avoid division by zero
    norm_colour = colour / norm
    dot_product = np.sum(norm_colour * v_obj, axis=2)

    return dot_product

def rgb_euclid(object_id, mask=None, verbose=False, abspath = ''):
    ''' cr
    Create an RGB image for an Euclid galaxy using their object ID

    PARAMETERS:
    ----------
        object_id: str
            ID of the object in the Euclid database
        verbose: bool
            Print information about the image
    
    RETURNS:
    -------
        enhanced_img: np.array
            RGB image of the galaxy
    '''

    y = join(abspath,f'EUC_SEGMID_{object_id}_NIR-Y_cutout.fits')
    j = join(abspath,f'EUC_SEGMID_{object_id}_NIR-J_cutout.fits')
    h = join(abspath,f'EUC_SEGMID_{object_id}_NIR-H_cutout.fits')
    vis = join(abspath,f'EUC_SEGMID_{object_id}_VIS_cutout.fits')

    vis_data = fits.getdata(vis)
    _,bkg,_ = sigma_clipped_stats(vis_data, mask=mask)
    vis_data = vis_data - bkg
    header_vis = fits.getheader(vis)
    zp_vis = header_vis['MAGZERO']
    pixscale_vis = get_pixel_scale(header_vis)

    zp_L = zp_vis + 5*np.log10(pixscale_vis)
    mag_L = zp_L - 2.5*np.log10(vis_data)

    rgb = np.zeros((vis_data.shape[0], vis_data.shape[1],3))
    zp = 30.0 
    for i,file in enumerate([h, j, y]):
        head = fits.getheader(file)
        data = fits.getdata(file) 
        zp_im = head['MAGZERO']
        _,med,std = sigma_clipped_stats(data, mask=mask)
        smooth = fabada(data-med,(1.2*std)**2)
        rgb[:,:,i] = smooth*10**((zp-zp_im)*0.4)
        if verbose: print(f'{med} - {std} - {np.nanmin(rgb[:,:,i])} - {np.nanmax(rgb[:,:,i])} - {np.nanmean(rgb[:,:,i])} - {np.nanstd(rgb[:,:,i])}')

    transform = LogStretch() + PercentileInterval(99.99)
    rgb = transform(rgb)

    rgb_lup = make_lupton_rgb(rgb[:,:,0], rgb[:,:,1], rgb[:,:,2], Q=1, stretch=0.9)
      

    ## SATURATION ADJUSTMENT
    satadj = 3.5
    # convert the RGB image to HSV format
    imghsv = cv2.cvtColor(rgb_lup, cv2.COLOR_BGR2HSV).astype("float32")
    # Apply saturation adjustment on the s channel
    (h, s, v) = cv2.split(imghsv)
    s = s*satadj
    s = np.clip(s,0,255)
    imghsv = cv2.merge([h,s,v])
    enhanced_img = cv2.cvtColor(imghsv.astype("uint8"), cv2.COLOR_HSV2BGR)

    ## LAB COLOR SPACE blue colour adjustment
    lab = cv2.cvtColor(enhanced_img, cv2.COLOR_BGR2LAB)
    l_channel, a, b = cv2.split(lab)
    # Applying CLAHE to L-channel
    clahe = cv2.createCLAHE(clipLimit=1, tileGridSize=(4,4))
    cl = clahe.apply(l_channel)
    # merge the CLAHE enhanced L-channel with the a and b channel
    limg = cv2.merge((l_channel,a,b))
    # Converting image from LAB Color model to BGR color spcae
    enhanced_img = cv2.cvtColor(limg, cv2.COLOR_LAB2BGR)

    _,med_vis,std_vis = sigma_clipped_stats(vis_data)
    smooth_vis = fabada(vis_data,(1.4*std_vis)**2,verbose=False)
    mag_L = zp_L - 2.5*np.log10(smooth_vis)
    vis_limit = mag_limit(std_vis, Zp=zp_vis, omega=8*pixscale_vis, scale=pixscale_vis, n=1)
    Lmask = mag_L[:,:,np.newaxis]*np.ones_like(rgb_lup) > 27.5
    enhanced_img[Lmask] = 0*np.ones_like(enhanced_img)[Lmask]

    return enhanced_img

def LSBrgb_euclid(object_id, mask=None, verbose=False, abspath = ''):
    ''' 
    Create an RGB image LSB style for an Euclid galaxy using their object ID

    PARAMETERS:
    ----------
        object_id: str
            ID of the object in the Euclid database
        mask: np.array
            Mask of the image to use to characterize the noise in the image
        abspath: str
            Absolute path to the directory where the FITS files are located
        verbose: bool
            Print information about the image
    
    RETURNS:
    -------
        enhanced_img: np.array
            RGB image of the galaxy
    '''
    y = join(abspath,f'EUC_SEGMID_{object_id}_NIR-Y_cutout.fits')
    j = join(abspath,f'EUC_SEGMID_{object_id}_NIR-J_cutout.fits')
    h = join(abspath,f'EUC_SEGMID_{object_id}_NIR-H_cutout.fits')
    vis = join(abspath,f'EUC_SEGMID_{object_id}_VIS_cutout.fits')

    vis_data = fits.getdata(vis)
    xc,yc,rad = vis_data.shape[1]//2, vis_data.shape[1]//2,vis_data.shape[0]/5
    rgb = np.zeros((vis_data.shape[0], vis_data.shape[1],3))

    zp = 30.0 
    luminance = 0
    for i,file in enumerate([h, j, y, vis]):
        head = fits.getheader(file)
        data = fits.getdata(file) 
        zp_im = head['MAGZERO']
        _,bkg,std = sigma_clipped_stats(data, mask=mask)
        smooth = fabada(data-bkg,std**2,verbose=True)
        if i<3: rgb[:,:,i] = smooth*10**((zp-zp_im)*0.4)
        luminance += smooth*10**((zp-zp_im)*0.4) 

    # Luminance
    luminance /= 4
    _,lumlim,_,_ = measureImageNoise(np.ma.masked_array(luminance,mask=mask),xc,yc,rad,0,0,20,30)
    _,_,stdlum = sigma_clipped_stats(luminance, mask=mask)
    smooth_lum = fabada(luminance, 5*stdlum**2,verbose=True)
    mag_L = 25 - 2.5*np.log10(smooth_lum)

    # Get limits for colour and grey images
    colormaglim = mag_limit(lumlim,Zp=30,scale=0.1,omega=0.1,n=4)
    greymaglim = mag_limit(lumlim,Zp=30,scale=0.1,omega=0.1,n=0.5)

    # Create the white bacground to increate contrast
    whites = mag_L > greymaglim+0.02
    whites = cv2.dilate(whites.astype(np.uint8), np.ones((3, 3), np.uint8), iterations=1)
    whites = cv2.GaussianBlur(whites, (7, 7), sigmaX=0).astype(bool)

    # Create colour image using magnitude scaling
    mag = 25 - 2.5*np.log10(rgb[:,:,:3])
    mag_max = np.nanpercentile(mag,0.05)
    rgb_mag = np.clip(mag,mag_max,colormaglim)
    rgb_mag = np.abs(rgb_mag - colormaglim)/mag_max
    rgb_mag[np.isnan(rgb_mag)]=0

    # Streach using asinh trasnformation
    rgb_lup = make_lupton_rgb(rgb_mag[:,:,0], rgb_mag[:,:,1], rgb_mag[:,:,2],Q=0.1, stretch=0.5)
        
    # Apply the white background to the RGB image
    for i in range(3):
        rgb_lup[:,:,i][whites] = 255

    ## SATURATION ADJUSTMENT
    satadj = 5
    # convert the RGB image to HSV format
    imghsv = cv2.cvtColor(rgb_lup, cv2.COLOR_BGR2HSV).astype("float32")
    # Apply saturation adjustment on the s channel
    (h, s, v) = cv2.split(imghsv)
    s = s*satadj
    s = np.clip(s,0,255)
    imghsv = cv2.merge([h,s,v])
    enhanced_img = cv2.cvtColor(imghsv.astype("uint8"), cv2.COLOR_HSV2BGR)

    ## LAB COLOR SPACE blue colour adjustment
    lab = cv2.cvtColor(enhanced_img, cv2.COLOR_BGR2LAB)
    l_channel, a, b = cv2.split(lab)
    # Applying CLAHE to L-channel
    clahe = cv2.createCLAHE(clipLimit=3, tileGridSize=(8,8))
    cl = clahe.apply(l_channel)
    # merge the CLAHE enhanced L-channel with the a and b channel
    limg = cv2.merge((l_channel,a,b))
    # Converting image from LAB Color model to BGR color spcae
    enhanced_img = cv2.cvtColor(limg, cv2.COLOR_LAB2BGR)

    return enhanced_img


def mask_cmap(alpha=0.5, color='red'):
    transparent = colorConverter.to_rgba('white',alpha = 0)
    gray = colorConverter.to_rgba(color,alpha = alpha)
    cmap = ListedColormap([transparent, gray])
    return cmap

def ellipse_points(center, a, b, angle, num_points=300):
    """
    Generate x and y coordinates for an ellipse.

    Parameters:
    center (tuple): (x, y) coordinates of the ellipse center.
    a (float): Semimajor axis.
    b (float): Semiminor axis.
    angle (float): Rotation angle in degrees.
    num_points (int): Number of points to generate.

    Returns:
    tuple: (x, y) coordinates of the ellipse.
    """
    theta = np.linspace(0, 2 * np.pi, num_points)
    x_circle = a * np.cos(theta)
    y_circle = b * np.sin(theta)

    angle_rad = np.deg2rad(angle)
    cos_angle = np.cos(angle_rad)
    sin_angle = np.sin(angle_rad)

    x_ellipse = cos_angle * x_circle - sin_angle * y_circle + center[0]
    y_ellipse = sin_angle * x_circle + cos_angle * y_circle + center[1]

    return x_ellipse, y_ellipse

def rectangle_add_patches(coords, width, height, ax, **kwargs):
    ''' Add rectangles patches to axis

    Parameters:
    ----------
        coord: array
            Array of coordinates (x,y) of the center of the rectangles
        width: float
            Width of the rectangles in pixels
        height: float
            Height of the rectangles in pixels
        ax: matplotlib axis
            Axis to add the patches
    
    Returns:
    --------
        True: bool
            True if the rectangles were added to the axis
    '''
    x = coords[0]
    y = coords[1]
    corner_x = x - width/2
    corner_y = y - height/2

    for cx,cy in zip(corner_x,corner_y):
        rect = Rectangle((cx,cy),width,height, **kwargs)
        ax.add_patch(rect)
    
    return True

def make_serializable(obj):
        """
        Recursively convert non-serializable objects to YAML-compatible types.
        """
        if isinstance(obj, np.ndarray):
            return obj.tolist()  # Convert numpy arrays to lists
        elif isinstance(obj, dict):
            return {key: make_serializable(value) for key, value in obj.items()}
        elif isinstance(obj, list):
            return [make_serializable(item) for item in obj]
        elif isinstance(obj, float):
            return float(obj)  # Ensure floats are properly handled
        else:
            return obj  # Return the object as is if it's already serializable

def write_pickle(data, file_path):
    with open(file_path, 'wb') as f:
        pickle.dump(data, f)
    return True

def write_yaml(data, file_path):
    """
    Save a dictionary to a YAML file, ensuring all data is YAML-serializable.

    Args:
        data (dict): The dictionary to save.
        file_path (str): The path to the YAML file.
    """
    # Convert the dictionary to a serializable format
    serializable_data = make_serializable(data)

    # Save the serializable dictionary to a YAML file
    with open(file_path, 'w') as yaml_file:
        yaml.dump(serializable_data, yaml_file, default_flow_style=False)

    return True

def eucl_name_formater(ra, dec):
    '''Given the RA and DEC of a source returns 
    the name formatted to Euclid standards, following IAU convention. 
    
    Parameters:
    -----------
        ra (float):
            Right Ascension in degrees
        dec (float): 
            Declination in degrees
    
    Returns:
    --------
        text (str): 
            Name formatted to Euclid standards in LaTeX
        '''
    while isinstance(ra, Iterable):
        ra = ra[0]
    while isinstance(dec, Iterable):
        dec = dec[0]
    sky_coord = SkyCoord(ra=ra*u.degree, dec=dec*u.degree)
    ra_hours = sky_coord.ra.hms
    dec_hours = sky_coord.dec.dms
    text = 'EUCL\,'
    text+=f'J{np.int8(ra_hours[0]):02d}{np.int8(ra_hours[1]):02d}'
    split = str(ra_hours[2]).split('.')
    text+=f'{int(split[0]):02d}.'+split[1][:2]
    text+='$-$' if dec_hours[0]<0 else '$+$'
    text+=f'{np.abs(dec_hours[0]).astype(np.int8):02d}{np.abs(dec_hours[1]).astype(np.int8):02d}'
    split = str(np.abs(dec_hours[2])).split('.')
    text+=f'{int(split[0]):02d}.'+split[1][:1]

    return text

def savitzky_golay_with_padding(array, window_length, polyorder):
    """
    Applies the Savitzky-Golay filter to a 1D array with padding to avoid edge effects.

    Parameters:
    - array: np.ndarray, the input 1D signal.
    - window_length: int, the length of the filter window (must be odd and >= polyorder + 2).
    - polyorder: int, the order of the polynomial to fit.

    Returns:
    - smoothed_array: np.ndarray, the smoothed signal with the same length as the input.
    """
    if window_length % 2 == 0 or window_length <= polyorder:
        raise ValueError("Window length must be odd and greater than polyorder.")
    
    # Reflect padding on both sides
    pad_size = window_length // 2
    padded_array = np.pad(array, pad_width=pad_size, mode='reflect')
    
    # Apply Savitzky-Golay filter
    smoothed_padded = savgol_filter(padded_array, window_length=window_length, polyorder=polyorder)
    
    # Remove padding to get the original signal length
    smoothed_array = smoothed_padded[pad_size:-pad_size]
    
    return smoothed_array


def redshift_surface_brightness_dimming(z):
    '''Calculate the surface brightness dimming to current redshift.
    Surface density scales as (1+z)^-3 so the surface brightness scales as
    7.5xlog10(1+z). 

    The rest frame mu_rest = mu_observed + dimming
    
    Parameters
    -----------
        z: float
            Redshift of the source
    
    Returns
    --------
        dimming: float
            Surface brightness dimming in magnitudes
    '''
    return -3*2.5*np.log10(1+z)

def break_type_string(number):
    '''Given the number of the classification
    convert to string format. 
    '''
    classify = {0: 'Type 0', 1: 'Type I',
        2: 'Type II', 3: 'Type III', 4: 'Type II+II',
        5: 'Type II+III',6: 'Type III+II',
        7: 'Type III+III',-1: 'Null'}
    out = classify[number] if number in classify.keys() else 'Null'
    return out


def euclid_to_besta(band):
    besta = {'CFIS-R': 'CFHT_MegaCam.r',
        'CFIS-U': 'CFHT_MegaCam.u',
        'DES-G': 'CTIO_DECam.g',
        'DES-I': 'CTIO_DECam.i',
        'DES-R': 'CTIO_DECam.r',
        'DES-Z': 'CTIO_DECam.z',
        'NIR-H': 'Euclid_NISP.H',
        'NIR-J': 'Euclid_NISP.J',
        'NIR-Y': 'Euclid_NISP.Y',
        'PANSTARRS-I': 'PANSTARRS_PS1.i',
        'VIS': 'Euclid_VIS.vis',
        'WISHES-G': 'Subaru_HSC.g',
        'WISHES-Z': 'Subaru_HSC.z'}
    return besta[band]

def deep_update(d, u):
    for k, v in u.items():
        if isinstance(v, collections.abc.Mapping):
            d[k] = deep_update(d.get(k, {}), v)
        else:
            d[k] = v
    return d

def create_yaml_config(filename='euclid_config.yaml', config_dict=None):
    """
    Writes the euclid_config.yaml file with comments and blank lines preserved.
    If config_dict is provided, its values will override the defaults.
    """
    # Default values
    defaults = {
        "pipeline_steps": {
            "M1": True,
            "M2": True,
            "M3": True,
            "M4": True,
            "M5": True,
        },
        "catalogs": {
            "objects": "/media/team_workspaces/euclid-surface-photometry-pipeline/code/test_euclid_surf/selected_test.csv",
            "ids": "/media/team_workspaces/euclid-surface-photometry-pipeline/code/test_euclid_surf/selected_test.csv",
        },
        "absolute_paths": {
            "data_path": "/data/user/euc_repository_idr_iqr1/Q1_R1/MER",
            "output_path": "/media/team_workspaces/euclid-surface-photometry-pipeline/storage/test_results",
        },
        "parameters": {
            "n_kron": 3,
            "max_size": 1000,
            "sb_threshold": 27.5,
        },
        "logging": {
            "level": "INFO",
            "file": "pipeline.log",
            "format": "%(asctime)s [%(levelname)-8s] %(processName)s - %(message)s",
            "datefmt": "%Y-%m-%d %H:%M:%S",
        }
    }

    # If config_dict is provided, update defaults
    if config_dict:
        config = deep_update(copy.deepcopy(defaults), config_dict)
    else:
        config = defaults

    yaml_content = f"""# Yaml file to run the Euclid surface photometry pipeline described 
# in Euclid Consortium: Sanchez-Alarcon et al. (2026)

pipeline_steps:
  M1: {config['pipeline_steps']['M1']}            # Module 1: Creates cutout of galaxy and mask
  M2: {config['pipeline_steps']['M2']}            # Module 2: Creates master profile from VIS cutout
  M3: {config['pipeline_steps']['M3']}            # Module 3: Measures the profile on NISP and External cutouts
  M4: {config['pipeline_steps']['M4']}            # Module 4: Measures parameters on the profiles
  M5: {config['pipeline_steps']['M5']}            # Module 5: Performs photometry and Disc modeling

catalogs:
  objects: "{config['catalogs']['objects']}"   # This is the catalog with the parameters of the objects to be processed. It should contain at least the columns: 'id', 'ra', 'dec', 'z'
  ids: "{config['catalogs']['ids']}"       # This is the catalog with the ids of the objects to be processed. It should contain at least the column 'id' with the same ids as in the 'objects' catalog. If not provided, all objects in the 'objects' catalog will be processed.

absolute_paths:
  data_path: "{config['absolute_paths']['data_path']}"                     # Path to the directory containing the data. 
  output_path: "{config['absolute_paths']['output_path']}"     # Path to the directory where the results will be stored.

parameters:
  n_kron: {config['parameters']['n_kron']}             ## Number of Kron radii to use for the cutout
  max_size: {config['parameters']['max_size']}        ## Maximum size of the cutout in pixels
  sb_threshold: {config['parameters']['sb_threshold']}    ## (29.75 at z=1) this is defined at z=0.

logging:
  # Logging level: DEBUG | INFO | WARNING | ERROR | CRITICAL
  level: {config['logging']['level']}
  # Log file path; set to null to disable file logging
  file: {config['logging']['file']}
  # Format string for log records
  format: "{config['logging']['format']}"
  datefmt: "{config['logging']['datefmt']}"
"""
    with open(filename, "w") as f:
        f.write(yaml_content)

    return os.path.isfile(filename), config