# pipelines.py: Main pipeline modules for Euclid surface photometry

from .utils import *
from .morphology import *

import os
import yaml
import json
import pickle

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


def euclid_module1(object_id, ra, dec, size, inPath, outPath, addMERback=True):
    '''Euclid Pipeline 1.
    Given an object ID of Euclid, create a cutout of 5xkron_radius, and 
    gerenates a mask for that source using Euclid SEGMENT, NoiseChisel and
    MTObjecst for deblending sources inside the galaxy.
    
    Fortmat of cutout: EUC_SEGMID-{object_id}_cutout.fits
    Fortmat of mask:   EUC_SEGMID-{object_id}_mask.fits
    
    PARAMETERS
    ----------
        object_id : str or int
            Euclid's object ID of the source
        ra : float
            Right Ascension of the source
        dec : float
            Declination of the source
        size : float
            Size of the cutout in degrees
        inPath : str
            Path to the input data (MER mosaics)
        outPath : str
            Path to save the cutout and mask
        addMerBack : bool
            If True, add the MER background to the cutout.
    
    RETURNS 
    -------
        cutoutFile : str
            Path to the cutout file
        maskFile : str
            Path to the mask file
    
    '''

    files = find(inPath, f'EUC_MER_BGSUB-MOSAIC*.fits')
    _cutout_ = join(outPath,f'EUC_SEGMID_{object_id}_XXXX_cutout.fits')

    # Create cutouts for each filter
    for f in files:
        inst = basename(f).split('_')[2].split('IC-')[-1]
        cutoutFile = _cutout_.replace('XXXX', inst)
        _ = cutout(f, (ra,dec), (size, size), hdu=0, mode='wcs', out=cutoutFile)
        if addMERback:
            f = find(inPath, f'EUC_MER_BGMOD-{inst}*fits')[0]
            _ = cutout(f, (ra,dec), (size, size), hdu=0, mode='wcs', out=cutoutFile.replace('.fits','_bg.fits'))
            cutoutData = fits.getdata(cutoutFile) + fits.getdata(cutoutFile.replace('.fits','_bg.fits'))
            cutoutHeader = fits.getheader(cutoutFile)
            savefits(cutoutData, cutoutHeader, cutoutFile, keywords=[('BGMOD',1,'MER Background added if 1')])
            del cutoutData, cutoutHeader

    # Create cutout of the segmentation map
    segFile = find(inPath, f'*FINAL-SEGMAP*.fits')[0]
    segCut = join(outPath,f'EUC_SEGMAP_{object_id}_segmap.fits')
    _ = cutout(segFile, (ra,dec), (size,size), hdu=0, mode='wcs', out=segCut)
    
    # Define cutout to output, i.e., VIS.
    cutoutFile = _cutout_.replace('XXXX','VIS')

    # Redefine center of the image
    image = Image(cutoutFile, zp=24.6)
    image.obj(ra,dec)
    new_center = center_hillclimb(image.data, image.pix, centeringring=np.int64(3600*size/2))
    sky = image.pixel_to_sky(*new_center)
    ra,dec  = sky.ra.value,sky.dec.value

    # Create the mask 
    maskFile = join(outPath,f'EUC_SEGMID_{object_id}_mask.fits') 
    _ = euclid_mask(outPath, (ra,dec), out=maskFile)

    return cutoutFile,maskFile,(ra,dec)

def euclid_module2(cutoutFile, maskFile, ra, dec, object_id):
    ''' Euclid Pipeline 2:
    Measure the profile of the galaxy in the VIS band to 
    have the master profile to measure in all images. 

    Result file: EUC_SEGMID-{object_id}_master-profile.fits
    
    PARAMETERS
    ----------
        cutoutFile : str
            Path to the VIS cutout file
        maskFile : str
            Path to the mask file
        ra : float
            Right Ascension of the source
        dec : float
            Declination of the source
    
    RETURNS
    -------
        result : Table
            Table with the master profile of the source
    '''
    __file__ = f"EUC_SEGMID_{object_id}_master-profile.csv"
    image = Image(cutoutFile, zp=24.6)
    image.obj(ra,dec)
    mask = fits.getdata(maskFile)
    _,bkg, bkgstd = sigma_clipped_stats(image.data, mask=mask) 
    mask[mask==mask[image.y,image.x]] = 0
    image.set_mask(mask)
    image.get_morphology(nsigma=3)
    result = autoprof_isophote_photometry(image.data - bkg, (image.x, image.y), image.pa, image.eps, growth=0.03,
                        fit_limit=0.3, smooth=2, bkgstd = bkgstd, psf=3)
    meta = {'bkg': bkg, 'bkgstd': bkgstd, 'file': cutoutFile, 'ra': ra, 'dec': dec, 
                    'pa': image.pa, 'eps': image.eps, 'x': image.x, 'y': image.y}
    for key,value in meta.items():
        result.meta[key] = value
    result.write(join(dirname(cutoutFile),__file__), format='ascii.ecsv', overwrite=True)
    return result

def euclid_module3(object_id, master, maskFile):
    '''Euclid Pipeline 3.
    Read master profile and measure it for the different bands.

    Result file: EUC_SEGMID_{object_id}_profiles.csv

    @TODO:
        Once the maximum radius if measured, crop the array
        allocated in memory prior to give it to the function
        to measure the profile. 
    
    PARAMETERS
    ----------
        object_id : int
            Euclid's object ID of the source
        master : Table
            Master profile of the source [from autoprof]
        maskFile : str
            Path to the mask file
        pixelscale : float
            Pixel scale of the image in arcsec/pixel
    
    RETURNS
    -------
        table : Table
            Table with the profiles of the source in different bands
    '''

    __file__ = f"EUC_SEGMID_{object_id}_profiles.csv"

    files = find(dirname(maskFile),f'*{object_id}*cutout.fits')
    
    
    mask = fits.getdata(maskFile)
    mask[mask==mask[np.int16(master['y'][0]), np.int16(master['x'][0])]] = 0
    
    tab,meta = {},{}

    size = len(master)
    pa =  np.nanmedian(master['fit pa'][-size//3:]) * 180/np.pi
    eps = np.nanmedian(master['fit ellip'][-size//3:])
    center = (master['x'][0], master['y'][0])
    
    # Measure backgrounds
    for i,file in enumerate(files):
        inst = basename(file).split('_')[3]
        data = np.ma.masked_array(fits.getdata(file), mask=mask)
        header = fits.getheader(file)  

        try: # New background estimation method
            bkgresults = background_estimation_euclid(data, center, pa, eps, init=0.8*master['fit R'][-1], 
                                     plot=file.replace('.fits','_bkg.jpg'))
            # write_yaml(bkgresults,file.replace('.fits','_bkg.yaml'))
            write_pickle(bkgresults,file.replace('.fits','_bkg.pickle'))
            bkg = bkgresults['ellip_bkg']
            bkgstd = bkgresults['rect_bkgstd']
            bkg_rad = bkgresults['bkgrad']

        except:
            _, med, std = sigma_clipped_stats(data)
            bkg,bkgstd,bkg_rad = med,std/np.sqrt(np.sum(mask)/10), np.nan
        
        meta[f'bkg_{inst}'] = bkg
        meta[f'bkgstd_{inst}'] = bkgstd
        meta[f'bkgrad_{inst}'] = bkg_rad
        meta[f'file_{inst}'] = file
        meta[f'zp_{inst}'] = header['MAGZERO']
        meta[f'pixscale_{inst}'] = get_pixel_scale(header)   
                        

    radmax = np.nanmax([meta[k] for k in meta.keys() if 'bkgrad_' in k])

    # Extend the master profile to the radius to the background measurement
    growth = np.nanmean(master['fit R'][1:]/master['fit R'][:-1])
    if radmax > growth*master['fit R'][-1]:
        alpha = np.log10(growth)
        n = np.log10(radmax/master['fit R'][-1])//alpha + 1
        new_rad = master['fit R'][-1]*10**(alpha*np.arange(1,n+1))
        new_rows = {}
        new_rows['fit R'] = new_rad
        for key in master.colnames:
            if 'fit R'!=key:
                new_rows[key] = np.interp(new_rad, master['fit R'], master[key])
        master = vstack([master, Table(new_rows)])


    # Measure the profiles in the different bands 
    for i,file in enumerate(files):
        inst = basename(file).split('_')[3]
        data = np.ma.masked_array(fits.getdata(file), mask=mask)
        header = fits.getheader(file)

        profile = elliptical_radial_profile(data, master['fit R'].value, (master['x'].value, master['y'].value), 
                                master['fit pa'].value*180/np.pi,master['fit ellip'].value)
        
        if i==0:
            tab['sma'] = profile.rad
            tab['pa'] = profile.pa
            tab['pa_err'] = np.append(master['fit pa_err'].value*180/np.pi,np.zeros(len(profile.rad)-len(master['fit pa_err'])))
            tab['eps'] = profile.eps
            tab['eps_err'] = np.append(master['fit ellip_err'].value,np.zeros(len(profile.rad)-len(master['fit ellip_err'])))
            tab['fmode_A2'] = np.append(master['fit Fmode A2'].value,np.zeros(len(profile.rad)-len(master['fit Fmode A2'])))
            tab['fmode_Phi2'] = np.append(master['fit Fmode Phi2'].value,np.zeros(len(profile.rad)-len(master['fit Fmode Phi2'])))
            tab['fmode_A4'] = np.append(master['fit Fmode A4'].value,np.zeros(len(profile.rad)-len(master['fit Fmode A4'])))
            tab['fmode_Phi4'] = np.append(master['fit Fmode Phi4'].value,np.zeros(len(profile.rad)-len(master['fit Fmode Phi4'])))
            tab['x'] = profile.x
            tab['y'] = profile.y
            tab['npixels'] = profile.npixels
            meta[f'mask'] = maskFile

        tab[f'int_{inst}'] = profile.int
        tab[f'intstd_{inst}'] = profile.intstd
        tab[f'flux_{inst}'] = profile.flux
        tab[f'fluxstd_{inst}'] = profile.fluxstd
        
    table = Table(tab)
    for m in meta.keys():
        table.meta[m] = meta[m]
    
    table.write(join(dirname(maskFile),__file__),format='ascii.ecsv',delimiter=',', overwrite=True)
    
    return table # Return success status and result dictionary


def euclid_module4(profiles, object_id, outPath, z=None):
    '''Euclid Pipeline 4.
    Measure the resolved photometry on the surface brightness profiles of different bands.

    Result file: EUC_SEGMID_{object_id}_photometry.csv

    TODO:
        - Add SB limit like when reaching 1 sigma error.
    
    PARAMETERS
    ----------
        profiles : Table
            Table with the profiles of the source in different bands
        object_id : int
            Euclid's object ID of the source
        outPath : str
            Path to save the photometry
    
    RETURNS
    -------
        done : bool
            Success status of the pipeline'''

    __file__ = f"EUC_SEGMID_{object_id}_photometry.csv"
    bands = np.unique([k.split('_')[-1] for k in profiles.colnames if 'int' in k])

    dimming = redshift_surface_brightness_dimming(z) if z != None else 0
    
    table = {}
    table = {'object_id':object_id}
    for band in bands:
        profile = Profile()
        profile.set_params(radii=profiles['sma'], intensity=profiles[f'int_{band}'], 
        instensity_err=profiles[f'intstd_{band}'], flux=profiles[f'flux_{band}'],
        fluxstd=profiles[f'fluxstd_{band}'], center=(profiles['x'], profiles['y']), npixels=profiles['npixels'],
        pa=profiles['pa'], pastd=profiles['pa_err'], eps=profiles['eps'], epsstd=profiles['eps_err'],
        bkg=profiles.meta[f'bkg_{band}'], bkgstd=profiles.meta[f'bkgstd_{band}'],
        zp=profiles.meta[f'zp_{band}'], pixscale=profiles.meta[f'pixscale_{band}'])

        sma, cog = profile.curveOfGrowth()
        try:
            table[f'mag_{band}'] = profile.totalMagnitude(sma, cog)
        except:
            table[f'mag_{band}'] = np.nanmin(cog)
        table[f'reff_{band}'] = profile.fractionalRadius(table[f'mag_{band}'], fluxFrac=0.5)  # Half-light radius
        table[f'sb_eff_{band}'] = profile.surfaceBrightness(table[f'reff_{band}'])
        table[f'c82_{band}'] = profile.concentration(table[f'mag_{band}'])                    # C82 = 5log(R80/R20)
        table[f'c31_{band}'] = profile.fractionalRadius(table[f'mag_{band}'], fluxFrac=0.75) / profile.fractionalRadius(table[f'mag_{band}'], fluxFrac=0.25)
        table[f'rpetro_{band}'] = profile.petrosianRadius()
        table[f'r25p5_{band}'],table[f'pa25p5_{band}'],table[f'eps25p5_{band}'] = profile.isophotalRadius(25.5, returnMorph=True)    # Radius at mu=25.5 mag/arcsec^2
        table[f'r26p5_{band}'],table[f'pa26p5_{band}'],table[f'eps26p5_{band}'] = profile.isophotalRadius(26.5, returnMorph=True)     # Radius at mu=26.5 mag/arcsec^2
        table[f'r28_{band}'],table[f'pa28_{band}'],table[f'eps28_{band}'] = profile.isophotalRadius(28, returnMorph=True)           # Radius at mu=28 mag/arcsec^2

        # Inclination corrected values
        mean_eps = np.nanmedian(profile.eps[profile.rad*profile.pixscale > table[f'r25p5_{band}']])
        axisRatio = 1 - mean_eps
        incl_corr = -2.5*np.log10(axisRatio)
        for corr,suffix in zip([incl_corr, dimming, incl_corr+dimming],['corr','dim','corr_dim']):
            profile.mu += corr
            table[f'r25p5_{band}_{suffix}'] = profile.isophotalRadius(25.5, returnMorph=False)    
            table[f'r26p5_{band}_{suffix}'] = profile.isophotalRadius(26.5, returnMorph=False)    
            table[f'r28_{band}_{suffix}'] = profile.isophotalRadius(28, returnMorph=False)
            profile.mu -= corr
        table[f'axisRatio_{band}'] = axisRatio  

    format_table = {}
    for key in table:
        format_table[key] = [table[key]]
    Table(format_table).write(join(outPath,__file__),format='csv', overwrite=True)

    return table  # Return success status and result dictionary


def euclid_module5(profiles, photometry, object_id, outPath, z=0, sb_threshold=None):
    '''Euclid Pipeline 5.
    Given the surface brightness profiles and the photometry of the source,
    classify the breaks in the profile and measure the parameters of the galaxy models.
    
    Result file: EUC_SEGMID_{object_id}_breaks.csv

    PARAMETERS
    ----------
        profiles : Table
            Table with the profiles of the source in different bands
        photometry : Table
            Table with the photometry of the source in different bands
        object_id : int
            Euclid's object ID of the source
        outPath : str
            Path to save the breaks
        z : float
            Redshift of the source to apply threshold
        sb_threshold : float
            Surface brightness threshold at z=0 to apply to the profile
    RETURNS
    -------
        breaks : Table
            Table with the parameters of the models fitted. 
    '''

    __file__ = join(outPath,f"EUC_SEGMID_{object_id}_breaks.csv")

    row = {'object_id':object_id,
        'sersic_mue':0.0,'sersic_re':0.0, 'sersic_n':0.0, 'sersic_chisq': 0.0,'sersic_chisq_min': 0.0, 'sersic_chisq_max':0.0, 'sersic_chisq_mean':0.0, 'sersic_BIC':0.0,
       'sersic_inner_mue':0.0,'sersic_inner_re':0.0, 'sersic_inner_n':0.0, 'sersic_inner_chisq': 0.0,'sersic_inner_chisq_min': 0.0, 'sersic_inner_chisq_max':0.0, 'sersic_inner_chisq_mean':0.0, 'sersic_inner_r_in':0.0,'sersic_inner_r_out':0.0,
       'disc1_broken_flag':0, 'disc1_mu0':0.0, 'disc1_h0':0.0, 'disc1_r_break':0.0, 'disc1_h1':0.0, 'disc1_chisq':0.0, 'disc1_chisq_min':0.0, 'disc1_chisq_max':0.0, 'disc1_chisq_mean':0.0, 'disc1_r_in':0.0, 'disc1_r_out':0.0, 'disc1_BIC':0.0, 'disc1_BIC_ratio':0.0,
       'disc2_broken_flag':0, 'disc2_mu0':0.0, 'disc2_h0':0.0, 'disc2_r_break':0.0, 'disc2_h1':0.0, 'disc2_chisq':0.0, 'disc2_chisq_min':0.0, 'disc2_chisq_max':0.0, 'disc2_chisq_mean':0.0, 'disc2_r_in':0.0, 'disc2_r_out':0.0, 'disc2_BIC':0.0, 'disc2_BIC_ratio':0.0,
       'final_mue':0, 'final_re':0, 'final_n':0,'final_mu0':0, 'final_h1':0, 'final_rbreak1':0, 'final_h2':0, 'final_rbreak2':0, 'final_h3':0, 'final_rbreak3':0, 'final_h4':0, 
       'final_ndiscs':0, 'final_r_sersic':0.0,  
       'final_chisq':0.0, 'final_chisq_min':0.0, 'final_chisq_max':0.0, 
       'final_chisq_mean':0.0, 'final_chisq_std':0.0, 'final_BIC':0.0, 'final_BIC_ratio':0.0,
       'final_classification':-1, 
       'sb_threshold':0, 'rad_cut':0 , 'sb_cut':0.0}
    
    verbose = False

    dimming = redshift_surface_brightness_dimming(z) 
    row['sb_threshold'] = sb_threshold if sb_threshold != None else -1
    
    band = 'VIS'
    profile = Profile()
    profile.set_params(radii=profiles['sma'], intensity=profiles[f'int_{band}'], 
        instensity_err=profiles[f'intstd_{band}'], flux=profiles[f'flux_{band}'],
        fluxstd=profiles[f'fluxstd_{band}'], center=(profiles['x'], profiles['y']), npixels=profiles['npixels'],
        pa=profiles['pa'], pastd=profiles['pa_err'], eps=profiles['eps'], epsstd=profiles['eps_err'],
        bkg=profiles.meta[f'bkg_{band}'], bkgstd=profiles.meta[f'bkgstd_{band}'],
        zp=profiles.meta[f'zp_{band}'], pixscale=profiles.meta[f'pixscale_{band}'])

    profile.brightness()
    profile.mu = profile.mu + dimming
    totmag = photometry[f'mag_{band}']
    reff = profile.fractionalRadius(totmag, fluxFrac=0.5)  # Half-light radius
    r25p5 = photometry[f'r25p5_{band}']

    # Set limit of the profile and define width of the breaks
    if sb_threshold != None:
        rad_cut_args = (profile.mu > sb_threshold)*np.isfinite(profile.mu)
        if any(rad_cut_args):
            rad_cut = profile.rad[rad_cut_args][0] 
        else: 
            rad_cut_args = np.isnan(profile.mu) * (profile.rad*profile.pixscale > r25p5)
            rad_cut = 0.9*profile.rad[rad_cut_args][0] if any(rad_cut_args) else profile.rad[-1]
    else:
        rad_cut_args = np.isnan(profile.mu) * (profile.rad*profile.pixscale > r25p5)
        rad_cut = 0.9*profile.rad[rad_cut_args][0] if any(rad_cut_args) else profile.rad[-1]

    rad_lim = 1.1*rad_cut
    width = np.int16(0.04*rad_lim)
    width = width+1 if width%2 == 0 else width
    width = 5 if width < 5 else width
    rad_lim = profile.rad[-1] if rad_lim > profile.rad[-1] else rad_lim

    # Define the region of the profile to analyze (cut nans and values outside the limit)
    truncInd = (profile.rad <= rad_cut) * np.isfinite(profile.mu)
    radius = profile.rad[truncInd].value
    mu = profile.mu[truncInd].value
    mu_uerr = profile.upperr[truncInd]
    mu_lerr = profile.lowerr[truncInd]

    # Smooth the profile and the errors
    arcsec = radius*profile.pixscale
    rad_lim *= profile.pixscale
    rad_cut *= profile.pixscale
    smooth = savitzky_golay_with_padding(mu, width, 1)
    smooth[arcsec<1.5] = mu[arcsec<1.5]                 # Recover the inner region
    mu_uerr = (mu_uerr - mu) + smooth
    mu_lerr = (mu_lerr - mu) + smooth
    error =  medfilt(mu_uerr - mu_lerr,  kernel_size=width) + 0.01*smooth  
    row['rad_cut'] = rad_cut
    row['sb_cut'] = smooth[closest(arcsec,rad_cut)]
    

    # Calculate the derivative, the mean derivative and the cumulative difference
    weights = np.append(arcsec[1]-arcsec[0], np.diff(arcsec))
    dmudrad = medfilt(derivative(radius,smooth), kernel_size=width)
    mean_der = np.nansum(dmudrad*weights)/np.nansum(weights)
    mean_der = np.nanmedian(dmudrad)
    cum_sum = np.nancumsum(dmudrad-mean_der)

    # Detect peaks
    peaks, _ = find_peaks(cum_sum,distance=width, prominence=0.15)
    # Detect minima (by analyzing the negative of y)
    minima, _ = find_peaks(-cum_sum, distance=width, prominence=0.15)
    # Convert to arcsecs
    peaks_arcsec = arcsec[peaks]
    minima_arcsec = arcsec[minima]
    breaks = np.sort(np.append(peaks_arcsec, minima_arcsec))

    # Change point analysis removing inner region [testing]
    mean_der_out = np.nansum(dmudrad[arcsec>peaks_arcsec[0]]*weights[arcsec>peaks_arcsec[0]])/np.nansum(weights[arcsec>peaks_arcsec[0]])
    dmudrad_cut = dmudrad + 0
    dmudrad_cut[arcsec<peaks_arcsec[0]] = mean_der_out
    cum_sum_outer = np.nancumsum(dmudrad_cut-mean_der_out)
    
    if len(breaks) == 0:
        breaks = [rad_lim]

    if verbose:
        print('Rad_lim:   ',rad_lim,' Width: ',width)
        print('Peaks at:  ',peaks_arcsec)
        print('Minima at: ',minima_arcsec)
        print('Breaks at: ',breaks)

    models = {}
    singles = {}

    ## Define boundaries for the models
    sersicBounds = np.array([(1,50),(1e-2,rad_cut/2),(1e-2,8)]).T    # mue, r_e, n
    discBounds = np.array([(1,50),(1e-1,50)]).T                # mu_0, h

    ## Sersic fit
    # Check if its a sersic model
    initial_guess = [mu[0], reff, 1]
    sersic_popt, _ = curve_fit(sersic_sb, arcsec, smooth, p0=initial_guess, sigma=error, bounds=sersicBounds)
    
    ## Compute model and chi-square
    sersic_model_all = sersic_sb(arcsec, *sersic_popt)
    chi2_sersic_all = (smooth - sersic_model_all)**2 / error**2

    singles[0] = {'model': 'Sersic',
            'params': sersic_popt,
            'chi2': chi2_sersic_all,
            'value': sersic_model_all,
            'arg': np.ones_like(arcsec).astype(bool)}

    row['sersic_mue'] = sersic_popt[0]
    row['sersic_re'] = sersic_popt[1]
    row['sersic_n'] = sersic_popt[2]
    row['sersic_chisq'] = np.nansum(chi2_sersic_all)
    row['sersic_chisq_min'] = np.nanmin(chi2_sersic_all)
    row['sersic_chisq_max'] = np.nanmax(chi2_sersic_all)
    row['sersic_chisq_mean'] = np.nanmean(chi2_sersic_all)

    ## Complex Modeling
    # Only inner part
    innerRad = 1.1*peaks_arcsec[0] if 1.1*peaks_arcsec[0] > 5*profile.pixscale else 5*profile.pixscale
    inner = arcsec < innerRad
    initial_guess = [mu[0], reff, 1]
    initial_guess = sersic_popt
    sersic_popt, _ = curve_fit(sersic_sb, arcsec[inner], smooth[inner], p0=initial_guess, bounds=sersicBounds)

    ## Compute model and chi-square
    sersic_model = sersic_sb(arcsec, *sersic_popt)
    chi2_sersic = (smooth - sersic_model)**2 / error**2

    # Redefine boundary 
    _,medChi,stdChi = sigma_clipped_stats(chi2_sersic[inner])
    threshold = medChi + 2*stdChi
    evaluate = (chi2_sersic > threshold) * (arcsec>innerRad)
    if np.sum(evaluate) > 2: 
        transitions = [arcsec[evaluate][2]]
    else:
        transitions = [innerRad]

    inner = arcsec < transitions[-1]

    # Define models:
    modelNum = 0
    models[0] = {'model': 'Sersic', 
                'params': sersic_popt,
                'chi2': chi2_sersic,
                'value': sersic_model,
                'arg': inner}

    row['sersic_inner_mue'] = sersic_popt[0]
    row['sersic_inner_re'] = sersic_popt[1]
    row['sersic_inner_n'] = sersic_popt[2]
    row['sersic_inner_chisq'] = np.nansum(chi2_sersic[inner])
    row['sersic_inner_chisq_min'] = np.nanmin(chi2_sersic[inner])
    row['sersic_inner_chisq_max'] = np.nanmax(chi2_sersic[inner])
    row['sersic_inner_chisq_mean'] = np.nanmean(chi2_sersic[inner])
    row['sersic_inner_r_in'] = 0.0
    row['sersic_inner_r_out'] = transitions[-1]

    if verbose:
        print(f'Inner: 0.0 - {transitions[-1]:1.2f}, InnerRad: {innerRad:1.2f}')
        print(f'Guess: {initial_guess}')
        print(f'Sersic: mu0: {sersic_popt[0]:1.2f}, r_eff: {sersic_popt[1]:1.2f}, n: {sersic_popt[2]:1.2f}')
        print(f'Chi2: {np.nanmean(chi2_sersic[inner]):1.2e}')
    
    # Fit exponentials
    ndisc = 0
    for i in range(len(peaks_arcsec)+1):
        # Define limits
        if i < len(peaks_arcsec):
            if peaks_arcsec[i] < transitions[-1]:
                continue
            else:
                upperlim = peaks_arcsec[i]
                lowerlim = transitions[-1] 
        else:
            upperlim = rad_lim
            lowerlim = transitions[-1]

        # Check if the last limits are too close        
        if (upperlim - lowerlim < 0.1*rad_cut) or np.sum((arcsec > lowerlim)*(arcsec < upperlim)) < 10:
            continue

        ndisc += 1
        if ndisc > 2: 
            continue
        
        ## Fit broken exponentials
        # Compute the guess parameters
        rhalf = 0.5*(lowerlim+upperlim)
        outer = (arcsec >= lowerlim) * (arcsec < upperlim)
        h1_guess = 2.5/(np.log(10)*np.nanmedian((dmudrad)[outer * (arcsec<rhalf)]))*np.nanmedian(np.diff(arcsec[outer * (arcsec<rhalf)]))
        h2_guess = 2.5/(np.log(10)*np.nanmedian(dmudrad[outer * (arcsec>rhalf) * (arcsec<rad_lim)]))*np.nanmedian(np.diff(arcsec[outer * (arcsec>rhalf)* (arcsec<rad_lim)]))
        initial_guess = [mu[outer][0],h1_guess,rhalf,h2_guess] #mu_1, h_1, r_break, h_2
        
        # Fit the broken exponential using loss function
        result = minimize(loss, initial_guess, args=(arcsec[outer], smooth[outer], error[outer]),
                            bounds=([1,50],[0.1,25],[lowerlim,upperlim],[0.1,25]))
        mu_1,h_1,r_break,h_2 = result.x
        broken_model = broken_exponential_sb(arcsec, mu_1, h_1, r_break, h_2)
        chi2_broken = (smooth - broken_model)**2 / error**2
        
        if verbose: 
            print(f'Outer: {lowerlim:1.2f} - {upperlim:1.2f} Rhalf: {rhalf:1.2f}')
            print(f'Guess: {initial_guess}')
            print(f'Broken: mu_1: {mu_1:1.2f}, h_1: {h_1:1.2f}, r_break: {r_break:1.2f}, h_2: {h_2:1.2f}')
            print(f'Chi2: {np.nanmean(chi2_broken[outer]):1.2e} +- {np.nanstd(chi2_broken[outer]):1.2e}')

        # Fit one component disc to test against the broken exponential
        h1_guess = np.nanmin([h_1, h_2])
        initial_guess = [mu_1,h1_guess]         #mu, h
        initial_guess = [19,2]
        disc_popt, _ = curve_fit(exponential_disc_sb, arcsec[outer], smooth[outer], p0=initial_guess, bounds=discBounds, sigma=error[outer])
        disc_model = exponential_disc_sb(arcsec, *disc_popt)
        chi2_disc = (smooth - disc_model)**2 / error**2
        
        # Compute BIC
        Nouter = np.sum(outer)
        bic_broken = Nouter*np.log(np.nansum(chi2_broken[outer])/Nouter) + 4*np.log(Nouter)
        bic_disc = Nouter*np.log(np.nansum(chi2_disc[outer])/Nouter) + 2*np.log(Nouter)
        
        # Define the flag for the broken exponential
        brokenflag = (r_break > lowerlim) * (r_break < upperlim) * ((h_1/h_2 > 1.1) + (h_2/h_1 > 1.1))
        brokenflag *= (r_break-lowerlim > 0.1*rad_cut) * (upperlim-r_break > 0.1*(upperlim-lowerlim))
        brokenflag *= (np.sum(arcsec[outer]>r_break) > 10)*(np.sum(arcsec[outer]<r_break) > 10)
        brokenflag *= (bic_broken < bic_disc)
        row[f'disc{ndisc}_broken_flag'] = brokenflag + 0

        if verbose:
            print(f'BIC: Broken: {bic_broken:1.2e} Disc: {bic_disc:1.2e}')

        if brokenflag:
            modelNum += 1
            models[modelNum] = {'model': 'BrokenExp', 
                'params': np.array([mu_1, h_1, r_break, h_2]),
                'chi2': chi2_broken,
                'value': broken_model,
                'arg': outer,
                'BIC': bic_broken,
                'BIC_ratio': bic_broken/bic_disc}
            
            _,medChi,stdChi = sigma_clipped_stats(chi2_broken[outer])
            threshold = medChi + 2*stdChi
            evaluate = (chi2_broken > threshold) * (arcsec > upperlim)
            if np.sum(evaluate) > 2:
                transitions += [arcsec[evaluate][2]]
            else:
                transitions += [upperlim]
            
            # Save 
            row[f'disc{ndisc}_mu0'] = mu_1
            row[f'disc{ndisc}_h0'] = h_1
            row[f'disc{ndisc}_r_break'] = r_break
            row[f'disc{ndisc}_h1'] = h_2
            row[f'disc{ndisc}_chisq'] = np.nansum(chi2_broken[outer])
            row[f'disc{ndisc}_chisq_min'] = np.nanmin(chi2_broken[outer])
            row[f'disc{ndisc}_chisq_max'] = np.nanmax(chi2_broken[outer])
            row[f'disc{ndisc}_chisq_mean'] = np.nanmean(chi2_broken[outer])
            row[f'disc{ndisc}_r_in'] = lowerlim
            row[f'disc{ndisc}_r_out'] = transitions[-1]
            row[f'disc{ndisc}_BIC'] = bic_broken
            row[f'disc{ndisc}_BIC_ratio'] = bic_broken/bic_disc
            

        else:
            modelNum += 1
            models[modelNum] = {'model': 'Exp', 
                'params':disc_popt,
                'chi2': chi2_disc,
                'value': disc_model,
                'arg': outer,
                'BIC': bic_disc,
                'BIC_ratio': bic_disc/bic_broken}
            
            _,medChi,stdChi = sigma_clipped_stats(chi2_disc[outer])
            threshold = medChi + 2*stdChi
            evaluate = (chi2_disc > threshold) * (arcsec > upperlim)
            if np.sum(evaluate) > 2:
                transitions += [arcsec[evaluate][2]]
            else:
                transitions += [upperlim]

            # Save 
            row[f'disc{ndisc}_mu0'] = disc_popt[0]
            row[f'disc{ndisc}_h0'] = disc_popt[1]
            row[f'disc{ndisc}_r_break'] = 0.0
            row[f'disc{ndisc}_h1'] = 0.0
            row[f'disc{ndisc}_chisq'] = np.nansum(chi2_disc[outer])
            row[f'disc{ndisc}_chisq_min'] = np.nanmin(chi2_disc[outer])
            row[f'disc{ndisc}_chisq_max'] = np.nanmax(chi2_disc[outer])
            row[f'disc{ndisc}_chisq_mean'] = np.nanmean(chi2_disc[outer])
            row[f'disc{ndisc}_r_in'] = lowerlim
            row[f'disc{ndisc}_r_out'] = transitions[-1]
            row[f'disc{ndisc}_BIC'] = bic_disc
            row[f'disc{ndisc}_BIC_ratio'] = bic_disc/bic_broken

                
    # Define final model and fit all components together again
    mu0 = row['disc1_mu0']
    h = [row[f'disc1_h0'], row['disc1_h1'], row[f'disc2_h0'], row['disc2_h1']]
    h = [x for x in h if x != 0]
    multiple_discs_flag = len(h) > 1
    rbreaks = [row['disc1_r_break'], row['disc2_r_in'], row['disc2_r_break']]
    rbreaks = [x for x in rbreaks if x != 0] if multiple_discs_flag else []
    final_bounds = np.append(sersicBounds, discBounds[:,0][:,np.newaxis],axis=1)
    if multiple_discs_flag:
        compBounds = np.append(discBounds[:,1][:,np.newaxis]*np.ones(len(h)), np.array([(0,rad_cut)]*len(rbreaks)).T, axis=1)
    else:
        compBounds = discBounds[:,1][:,np.newaxis]*np.ones(len(h))
    final_bounds = np.append(final_bounds, compBounds, axis=1)
    
    finalModel = MultiComponentModel(row['sersic_inner_mue'],
                                        row['sersic_inner_re'],
                                        row['sersic_inner_n'], mu0, h, rbreaks)
    
    if len(h) > 0:
        disc_popt, _ = curve_fit(finalModel.function_to_fitter, arcsec, smooth, 
                                p0=finalModel.param_list, bounds=final_bounds)#, sigma=error)
        finalModel.update_params(disc_popt)
        finaldict = finalModel.get_dictionary(prefix='final_')
        final_mu = finalModel(arcsec)
    else:
        finaldict = {}
        final_mu = sersic_model_all

    # Compare with sersic
    # Compute chi-square of the combined model and BICs
    final_chisq_array = ((smooth - final_mu )**2 / error**2)[arcsec<rad_cut]
    final_chisq = np.nansum(final_chisq_array)
    Nvarys = len(finalModel.param_list)
    Nelements = np.sum(arcsec<rad_cut)
    bic_sersic = Nelements*np.log(np.nansum(singles[0]["chi2"][arcsec<rad_cut])/Nelements) + 3*np.log(Nelements)
    bic_final = Nelements*np.log(final_chisq/Nelements) + Nvarys*np.log(Nelements)
    
    # Insert parameters into the row
    for key in finaldict.keys():
        row[key] = finaldict[key]
    row['sersic_BIC'] = bic_sersic
    row['final_chisq'] = final_chisq
    row['final_chisq_min'] = np.nanmin(final_chisq_array)
    row['final_chisq_max'] = np.nanmax(final_chisq_array)
    row['final_chisq_mean'] = np.nanmean(final_chisq_array)
    row['final_chisq_std'] = np.nanstd(final_chisq_array)
    row['final_BIC'] = bic_final
    row['final_BIC_ratio'] = bic_final/bic_sersic
    row['final_flag'] = bic_final<bic_sersic

    #######################
    #   Classify breaks   #
    #######################
    ndisc = row['final_ndiscs']
    h1_h2 = np.array(row['final_h1'])/row['final_h2']
    h2_h3 = np.array(row['final_h2'])/row['final_h3']
    
    classification = 0
    if ndisc==1: classification = 1
    
    # 2 Discs
    if ndisc==2:
       classification = 2
       if h1_h2<1:
           classification = 3

    # 3 Discs
    if ndisc==3:
        if (h1_h2>1)*(h2_h3>1):
            classification = 4  # II+II
        elif (h1_h2>1)*(h2_h3<1):
            classification = 5  # II+III
        elif (h1_h2<1)*(h2_h3>1):
            classification = 6  # III+II
        elif (h1_h2<1)*(h2_h3<1):
            classification = 7  # III+III

    # 4 Disc (unphysical?)
    if ndisc==4: classification = -1

    row['final_classification'] = classification
    

    # Some printing if verbose
    pr ='Single models: '
    for key in singles.keys():
        pr +=f'{singles[key]["model"]}: {np.nansum(singles[key]["chi2"][arcsec<rad_lim]):1.2e}, '
    pr +='\n'
    pr += 'Models: \n'
    for key in models.keys():
        # pr += ' '*8+f'{models[key]["model"]}: {np.nanmedian(models[key]["chi2"][models[key]["arg"]]):1.2e},  '
        pr += ' '*8+f'{models[key]["model"]}: ' +', '.join([f'{x:.2f}' for x in models[key]["params"]])
        reg = arcsec[models[key]["arg"]]
        pr += f' [{np.min(reg):3.2f},{np.max(reg):3.2f}] arcsec\n'
    pr +='\n'   
    pr +=f'Final model: Bulge+{finalModel.ndiscs} Discs '+ '\n' + ', '.join([f'{k}:{v:.2f}' for k,v in finalModel.get_dictionary().items()])
    pr += '\n'+ 'BIC: Sersic: ' + f'{bic_sersic:1.2e} Final: {bic_final:1.2e} %={bic_final/bic_sersic:1.2f}'
    pr +='\n'+'-'*50 + '\n'
    text = f' Models: Bulge+{finalModel.ndiscs} Discs Chi2: {final_chisq:1.2e}\n' + pr
    if verbose:print(pr)

    ####################################
    # Plotting
    ####################################
    
    fig, (ax1,axbar, ax2,ax3) = plt.subplots(4,1,figsize=(13,15),sharex=True)
    ax1.plot(profile.rad*profile.pixscale, profile.mu,'bo',lw=2.3,ms=5,alpha=0.4,zorder=1)
    ax1.plot(arcsec, smooth,'r-',lw=3,zorder=2)
    ax1.fill_between(arcsec, mu_lerr, mu_uerr, color='r', alpha=0.3)
    if sb_threshold != None:
        ax1.axhline(sb_threshold, c='k', ls='--')
        ax1.text(0, sb_threshold-0.1, '$\mu_{\mathrm{thres}}='+f'{sb_threshold:.2f}$', va='bottom', ha='left',fontsize=12)
    
    ax1.set_xlabel('Radius [arcsec]')
    ax1.set_ylabel('$\mu$ [mag$\,$arcsec$^{-2}$]')
    ax1.invert_yaxis()
    ax1.set_title(f'Morphological Decomposition (Bulge + {finalModel.ndiscs} Discs)')

    axbar.plot(arcsec, (smooth-final_mu)**2/error**2,'b-',lw=2.3, label='Multi-add')
    axbar.set_ylabel('$\chi^2$')
    axbar.axhline(0, c='k', ls='--')
    axerr = axbar.twinx()
    axerr.plot(arcsec, error/smooth, 'r-', lw=2.3, label='Error')
    axerr.set_ylabel('Relative error')

    # Models
    ax1.plot(arcsec, singles[0]['value'],'g-.',lw=2.3,zorder=5)
    ax1.text(0.98,0.96, text, transform=ax1.transAxes, va='top', ha='right')
    ax1.plot(arcsec, final_mu, 'k--', lw=2, zorder=10)
    discs = finalModel.get_components(arcsec)
    for disc in discs:
        ax1.plot(arcsec, disc, '*', color='darkred', markersize=2, zorder=10)
    
    ax2.plot(arcsec, dmudrad,'b-',lw=2)
    ax2.axhline(mean_der,c='k',ls='--')
    ax2.set_xlabel('Radius [arcsec]')
    ax2.set_ylabel('d$\mu$/dR')
    
    ax3.plot(arcsec, cum_sum,'b-',lw=1.5)
    ax3.plot(arcsec, cum_sum_outer,'b-.',lw=1.5, alpha=0.8)
    ax3.plot(arcsec[peaks], cum_sum[peaks], "ro", label="Peaks")
    ax3.plot(arcsec[minima], cum_sum[minima], "bo", label="Minima")

    for ax in [ax2,ax3]:
        for p in peaks:
            ax.axvline(arcsec[p],  c="r",ls='--')
        for m in minima:
            ax.axvline(arcsec[m], c="b",ls='--')

    ax3.set_xlabel('Radius [arcsec]')
    ax3.set_ylabel('Cumulative sum')
    ax1.set_xlim(-0.03*rad_lim,1.1*rad_lim)
    ax1.set_ylim(1.1*np.nanmax(mu),0.9*np.nanmin(mu))
    fig.tight_layout()
    fig.savefig(__file__.replace('.csv','.png'))
    plt.close('all')
        
    # Create a figure
    fig = plt.figure(figsize=(10, 5))  # Adjust the figure size as needed
    gs = GridSpec(2, 1, height_ratios=[2, 1])  # 2:1 height ratio
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex=ax1)
    profile_kpc = arcsec_to_kpc(profile.rad*profile.pixscale,z)
    arcsec_kpc = arcsec_to_kpc(arcsec, z)
    ax1.plot(profile_kpc, profile.mu,'r*',lw=2.3,ms=5,alpha=0.4,zorder=1)
    ax1.plot(arcsec_kpc, smooth,'r-',lw=3,zorder=2)
    ax1.plot(arcsec_kpc, final_mu, 'k--', lw=2, zorder=10)
    ax1.set_xlim([-0.5, 100])
    ax1.set_ylim([15,31])
    ax1.set_ylabel('$\mu$ [mag$\,$arcsec$^{-2}$]')
    ax1.invert_yaxis()

    ax2.plot(arcsec_kpc, smooth-final_mu,'b-',lw=2.3, label='Residual')
    ax2.set_ylabel('$\mu$-model')
    ax2.axhline(0, c='k', ls='--')
    ax2.set_xlabel('Radius [kpc]')
    ax3 = ax2.twinx()
    ax3.plot(arcsec_kpc, (smooth-final_mu)**2/error**2, 'r-.', lw=2.3, label='$\chi^2$')
    ax3.set_ylabel('$\chi^2$')
    fig.tight_layout()
    fig.savefig(__file__.replace('.csv','_kpc.png'), dpi=200, bbox_inches='tight')
    plt.close('all')

    ## Saving results
    new_row = {key:[row[key]] for key in row.keys()}
    table = Table(new_row)
    table.write(__file__,format='ascii.ecsv',delimiter=',',overwrite=True)

    ## Save Pickle
    savePickle = {'arcsec':arcsec, 'mu':mu, 'smooth':smooth, 
            'rad_lim':rad_lim, 'rad_cut':rad_cut, 'width':width,
            'singles':singles, 'models':models, 'finalModel': finalModel.get_dictionary() }

    with open(__file__.replace('.csv','.pkl'), 'wb') as file:
        pickle.dump(savePickle, file)

    return True  # Return success status and result dictionary



