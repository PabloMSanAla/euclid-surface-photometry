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



class MultiComponentModel:
    ''' Class to generate a multi-component model
    with a bulge modeled by a Sersic and multiple exponential discs
    
    Attributes:
    -----------
        mue : float
            Surface brightness at the effective radius
        re : float
            Effective radius in units of r
        n : float
            Sersic index
        mu0 : float
            Central surface brightness of the first disc component
        h : float or array-like
            The scale lengths of the disc components
        rbreaks : float or array-like
            The radii at which the scale length changes
        
    Methods:
    --------
        __call__(r) : Generates the multi-component model
           
        function_to_fitter(r, *params) : Funtction to provide to curve_fit
           
    '''

    def __init__(self, mue, re, n, mu0, h, rbreaks):
        ''' Initialize the model with the parameters 
        
        Parameters:
        ----------
            mue : float
                Surface brightness at the effective radius
            re : float
                Effective radius in units of r
            n : float
                Sersic index
            mu0 : float
                Central surface brightness of the first disc component
            h : array-like
                The scale lengths of the disc components
            rbreaks : array-like
                The radii at which the scale length changes
        '''
        self.mue = mue
        self.re = re
        self.n = n
        self.mu0 = mu0
        self.h = h
        self.rbreaks = rbreaks
        self.ndiscs = len(h)
        self.nbreaks = len(rbreaks)
        self.param_list = [mue, re, n, mu0] + [x for x in h] + [x for x in rbreaks]


    def __call__(self, r):
        ''' Call the model with the parameters '''
        return multicomponent_model_piecewise(r, self.mue, self.re, self.n, self.mu0, self.h, self.rbreaks)
    
    def function_to_fitter(self, r, *params):
        ''' Function to provide to curve fit to fit the 
        multi-component model to the data. The parameters
        are passed as a list, and the function returns the
        model surface brightness profile.
        
        The parameters are:
        [sersic_mu, sersic_r_e, sersic_n, mu0, 
                h1, h2, ..., rbreaks1, rbreaks2, ...]
        
        Parameters:
        -----------
            r : array-like
                The radii at which the surface brightness will be computed
            params : list
                The parameters of the model. The first four are the Sersic parameters,
                and the rest are the disc parameters.
        Returns:
        --------
            mu : array-like
                The surface brightness profile
        '''
        size = len(params)
        sersic_mu, sersic_re, sersic_n, mu0 = params[:4]
        d = size - 4
        ndiscs = 1 + d//2
        h = params[4:4+ndiscs]
        rbreaks = params[4+ndiscs:]
        mu = multicomponent_model_piecewise(r, sersic_mu, sersic_re, sersic_n, mu0, h, rbreaks)
        return mu
    
    def update_params(self, params):
        ''' Update the parameters of the model '''
        self.mue = params[0]
        self.re = params[1]
        self.n = params[2]
        self.mu0 = params[3]
        self.h = params[4:4+self.ndiscs]
        self.rbreaks = params[4+self.ndiscs:]
        self.param_list = [self.mue, self.re, self.n, self.mu0] + [x for x in self.h] + [x for x in self.rbreaks]
    
    def get_components(self, r):
        ''' Get the disc components of the model '''
        disc_comp = []
        mu0 = self.mu0
        for i in range(self.ndiscs):
            disc_comp += [exponential_disc_sb(r, mu0, self.h[i])]
            if i < len(self.rbreaks):
                mu0 = mu0 + (2.5/np.log(10))*self.rbreaks[i]*(self.h[i+1]-self.h[i])/(self.h[i]*self.h[i+1])
        return disc_comp

    def get_dictionary(self,prefix=''):
        ''' Get the parameters of the model as a dictionary '''
        params = {
            f"{prefix}mue": self.mue,
            f"{prefix}re": self.re,
            f"{prefix}n": self.n,
            f"{prefix}mu0": self.mu0
        }
        for i in range(self.ndiscs):
            params[f"{prefix}h{i+1}"] = self.h[i]
            if i < self.nbreaks:
                params[f"{prefix}rbreak{i+1}"] = self.rbreaks[i]
        params[f"{prefix}ndiscs"] = self.ndiscs
        return params
    
    def __str__(self):
        ''' String representation of the model '''
        params = self.get_dictionary()
        return f"Multi-component model with {self.ndiscs} discs:\n" + str(params)


def exponential_disc_sb(r, mu_0, h):
    ''' Exponential disc in a surface brightness profile
    
    Parameters:
    -----------
        r : float
            Radius in arcsec
        mu_0 : float
            Central surface brightness in mag/arcsec^2
        h : float
            Scale length in arcsec
    
    Returns:
    --------
        mu : float
            Surface brightness at radius r'''
    a = 2.5*(1/np.log(10))
    return mu_0 + a*r*(1/h)

def bar_profile_sb(r, mu_1, h_1, r_bar, h_2, beta):
    ''' Broken exponential in a  surface brightness profile 
    following  Sánchez-Menguiano, L + 2017 
    (DOI: 10.1051/0004-6361/201731486)
    
    Parameters:
    -----------
        r : float
            Radius in arcsec
        mu_1 : float
            Surface brightness at r=0 of first profile
        h_1 : float
            Scale length at r=0 (slope of first profile)
        r_break : float
            Position of the break in arcsec
        h_2 : float
            Scale length at r=0 (slope of second profile)
        beta : float
            Transition parameter (default 1e-8)
    
    Returns:
    --------
        mu : float
            Surface brightness at radius r
    '''
    a = 2.5*(1/np.log(10))
    mu_2 = mu_1 + a*r_bar*(h_2-h_1)/(h_1*h_2)
    W = (1/np.pi)*(np.pi/2 + np.arctan2(r-r_bar, beta))
    return (mu_1 + a*r*(1/h_1))*(1-W) + (mu_2 + a*r*(1/h_2))*W
 

def broken_exponential_sb(r, mu_1, h_1, r_break, h_2):
    ''' Broken exponential in a  surface brightness profile 
    following  Sánchez-Menguiano, L + 2017 
    (DOI: 10.1051/0004-6361/201731486)
    
    Parameters:
    -----------
        r : float
            Radius in arcsec
        mu_1 : float
            Surface brightness at r=0 of first profile
        h_1 : float
            Scale length at r=0 (slope of first profile)
        r_break : float
            Position of the break in arcsec
        h_2 : float
            Scale length at r=0 (slope of second profile)
        beta : float
            Transition parameter (default 1e-8)
    
    Returns:
    --------
        mu : float
            Surface brightness at radius r
    '''
    beta=1e-8
    a = 2.5*(1/np.log(10))
    mu_2 = mu_1 + a*r_break*(h_2-h_1)/(h_1*h_2)
    W = (1/np.pi)*(np.pi/2 + np.arctan2(r-r_break, beta))
    return (mu_1 + a*r*(1/h_1))*(1-W) + (mu_2 + a*r*(1/h_2))*W
    

def sersic_sb(r, mu_e, r_e, n):
    ''' Sersic profile in surface brightness units
    
    Parameters:
    -----------
        r : float
            Radius in arcsec
        mu_e : float
            Surface brightness at r_e
        r_e : float
            Effective radius in arcsec
        n : float
            Sersic index
    
    Returns:
    --------
        mu : float
            Surface brightness at radius r'''
    bn = gammaincinv(2 * n, 0.5)
    return mu_e + 2.5*(1/np.log(10))*bn*((r/r_e)**(1/n) - 1)

def loss(p, r, y, error):
    ''' Loss function for fitting a broken exponential profile
    to a surface brightness profile
    
    Parameters:
    -----------
        p : list
            Parameters of the model [mu_1, h_1, r_break, h_2]
        r : array
            Radii in arcsec
        y : array
            Surface brightness in mag/arcsec^2
    
    Returns:
    --------
        loss : float
            Loss value'''
    rmin,rmax = np.nanmin(r), np.nanmax(r)
    alpha, delta = 1e-2,1e-4
    gamma = len(r)
    a = 2.5*(1/np.log(10))
    mu_1, h_1, r_break, h_2 = p
    npoints_end = np.sum(r>r_break)
    npoints_start = np.sum(r<r_break)
    # Get the model and measure the chi2
    y_pred = broken_exponential_sb(r, mu_1, h_1, r_break, h_2)
    chisq = np.sum((y-y_pred)**2/error**2)
    # Regularization term:
    regularization = alpha*(1/((h_1-h_2)**2 + delta))   
    regularization +=  gamma*((1/(1+np.exp(3*(r_break-rmin))))+(1/(1+np.exp(2*(rmax-r_break))))) # r_break inside
    regularization += 5*np.exp(-(npoints_end-7)/3)
    regularization += 5*np.exp(-(npoints_start-7)/3)
    return chisq + regularization


def disc_weights(r, rbreaks, beta):
    """
    Function to generate the softening function for the breaks in the disc

    Parameters
    ----------
    r : array-like
        The radii at which the softening function will be computed
    rbreaks : array-like
        The radii at which the scale length changes
    beta : float
        The softening parameter
    
    Returns
    -------
        Wb : array-like
            The weights for each disc components
    """
    rbreaks = np.array(rbreaks)
    r = np.array(r)
    Wb = (1/np.pi)*(np.pi/2 + np.arctan2(r[:,np.newaxis]-rbreaks[np.newaxis,:], beta))
    return Wb

def multiple_exponential_discs(r, mu0, h, rbreaks, modelsOut=False):
    """
    Function to generate multiple exponentials discs
    in surface brightness units (mag arcsec^-2)

    Parameters
    ----------
    r : array-like
        The radii at which the surface brightness will be computed
    mu0 : float
        The central surface brightness of the first component disc
    h : array-like
        The scale lengths of the disc components
    rbreaks : array-like
        The radii at which the scale length changes
    modelsOut : bool, optional
        If True, the function will return the surface brightness of each disc 
        component as well as the total surface brightness. Default is False.
    
    Returns
    -------
        mu : array-like
            The surface brightness profile

    """

    # Make sure that the inputs are numpy arrays
    r = np.array(r)
    h = np.array(h)
    rbreaks = np.array(rbreaks)

    # Check how many discs, and raise error if incompatible
    components = len(h)
    number_breaks = len(rbreaks)
    if components - number_breaks != 1:
        raise ValueError(f"The number of scale lengths ({components}) must be one more than the number of breaks ({number_breaks}) ")
    
    # Check that the number of components is at least 2
    if components <= 1:
        raise ValueError("At least two disc component must be provided")
    
    # Generate weights maps for each disc component
    Wb = disc_weights(r, rbreaks, 1e-6)  
    Wb = np.append(np.ones_like(r)[:,np.newaxis],np.append(Wb, np.zeros_like(r)[:,np.newaxis], axis=1),axis=1)

    # Compute the mu0 for each component of the disc
    a = 2.5*(1/np.log(10))
    factor = a*rbreaks*(h[1:]-h[:-1])/(h[:-1]*h[1:])
    mu0s = np.append(mu0, mu0 * np.ones_like(factor) + np.cumsum(factor))


    # Compute the surface brightness profile
    mu_pre = mu0s[np.newaxis,:] + a*r[:,np.newaxis]/h[np.newaxis,:]
    mu = np.sum(mu_pre*Wb[:,:-1]*(1-Wb[:,1:]), axis=1)
    output = (mu, mu_pre) if modelsOut else mu
    return output

def multicomponent_model_piecewise(r, mue, re, n, mu0, h, rbreaks):
    ''' Function to generate a multi-component model 
    with Sersic bulge and multiple exponential discs

    Parameters:
    -----------
    
        r : array-like
            The radii at which the surface brightness will be computed
        mue : float
            Surface brightness at the effective radius
        re : float
            Effective radius in units of r
        n : float
            Sersic index
        mu0 : float
            Central surface brightness of the first disc component
        h : float or array-like
            The scale lengths of the disc components
        rbreaks : float or array-like
            The radii at which the scale length changes

    Returns:
    --------
        mu : array-like
            The surface brightness
    '''
    # Define the different components of the model
    if isinstance(h, Iterable) and len(h) > 1:
        disc_comp = multiple_exponential_discs(r, mu0, h, rbreaks)
    else:
        disc_comp = exponential_disc_sb(r, mu0, h[0])
    bulge = sersic_sb(r, mue, re, n)
    mu_pre = np.append(bulge[:,np.newaxis], disc_comp[:,np.newaxis], axis=1)
    
    # Find the indices where the difference changes sign
    diff =  disc_comp - bulge
    cross_indices = np.where(np.diff(np.sign(diff)))[0]

    # Calculate the crossing points using linear interpolation
    x1, x2 = r[cross_indices], r[cross_indices + 1]
    y1, y2 = diff[cross_indices], diff[cross_indices + 1]
    m = (y2 - y1) / (x2 - x1)
    b = y1 - m * x1
    cross_points = -b / m
    
    # Define the radius of the bulge
    rbulge = np.nanmin(cross_points)


    # Generate weights maps for each disc component
    Wb = disc_weights(r, [rbulge], 1e-3)  
    Wb = np.append(np.ones_like(r)[:,np.newaxis],np.append(Wb, np.zeros_like(r)[:,np.newaxis], axis=1),axis=1)
    
    # Combine the components
    mu = np.sum(mu_pre*Wb[:,:-1]*(1-Wb[:,1:]), axis=1)

    return mu

def multicomponent_model(r, mue, re, n, mu0, h, rbreaks):
    ''' Function to generate a multi-component model 
    with a bulge modeled by a Sersic and multiple exponential discs

    Parameters:
    -----------
    
        r : array-like
            The radii at which the surface brightness will be computed
        mue : float
            Surface brightness at the effective radius
        re : float
            Effective radius in units of r
        n : float
            Sersic index
        mu0 : float
            Central surface brightness of the first disc component
        h : float or array-like
            The scale lengths of the disc components
        rbreaks : float or array-like
            The radii at which the scale length changes

    Returns:
    --------
        mu : array-like
            The surface brightness
    '''
    # Define the different components of the model
    if isinstance(h, Iterable) and len(h) > 1:
        disc_comp = multiple_exponential_discs(r, mu0, h, rbreaks)
    else:
        disc_comp = exponential_disc_sb(r, mu0, h[0])
    
    # Compute the bulge surface brightness
    bulge = sersic_sb(r, mue, re, n)

    # Find the indices where the difference changes sign
    diff =  disc_comp - bulge
    cross_indices = np.where(np.diff(np.sign(diff)))[0]

    # Calculate the crossing points using linear interpolation
    x1, x2 = r[cross_indices], r[cross_indices + 1]
    y1, y2 = diff[cross_indices], diff[cross_indices + 1]
    m = (y2 - y1) / (x2 - x1)
    b = y1 - m * x1
    cross_points = -b / m
    
    # Define the radius of the bulge
    rbulge = np.nanmin(cross_points) if len(cross_points)>0 else np.nanmax(r)

    # Combine the components by adding them in linear normalized_sigmoid(r,x0=2*rbulge)
    total_linear = 10**(-0.4*bulge)+ 10**(-0.4*disc_comp)
    mu = -2.5*np.log10(total_linear) 
    
    return mu

def normalized_sigmoid(x, x0=0):
    """
    Normalized sigmoid function.
    
    Parameters:
    -----------
    x : float or array-like
        Input value(s).
    
    Returns:
    --------
    float or array-like
        Normalized sigmoid output between 0 and 1.
    """
    return 1 - 1 / (1 + np.exp(-(x-x0)))