#plotter code
import numpy as np
from spectral_cube import SpectralCube
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider
from astropy.wcs import WCS
from mpl_toolkits.axes_grid1 import make_axes_locatable, axes_size
from astropy.io import fits 
import os

__all__ = ['Plotter', 'Load_fits', 'Interactive']

class Plotter():
    def __init__(self, cube_path):
        '''
        Initializes the class and inputs the 3D cube header that you wish to plot
        
        Parameters
        ----------
        cube_path: string
            the path to your 3D cube of data
            
        Returns
        -------
        print(self.cube): spectral cube
            returns the 
        
        '''
        cube = fits.open(cube_path)
        cube = SpectralCube.read(cube)
        self.cube = cube.with_spectral_unit(u.km / u.s)
        header = self.cube.header
        header1 = self.cube[0].header
        self.wcs_cube = WCS(header1)
        return print(self.cube)
        
    def spectrum(self):
        '''
        Plots the spectrum of your cube
        
        Returns
        -------
        fig, ax: matplotlib subplot
            the fig, ax for your spectrum plot
        
        '''
        Vsh, Ysh, Xsh = self.cube.shape
        chans = np.linspace(0,Vsh,num=Vsh)
        velo = np.zeros(Vsh)
        for i in range(0,Vsh):
            velo[i]= (self.cube.header['CRVAL3']) + ((i-self.cube.header['CRPIX3'])*self.cube.header['CDELT3'])
        spec=0.
        pix =0
        for xx in range(0,Xsh):
            for yy in range(0,Ysh):
                if np.all(np.isnan(self.cube[:,yy,xx]))==False:
                    spec = spec + self.cube[:,yy,xx]
                    pix = pix+1
                    
        fig, ax = plt.subplots(figsize=(8,6))
        ax.plot(velo, spec/pix, color='black')
        ax.grid()
        ax.set_xlabel('Velocity [km/s]')
        ax.set_ylabel('Flux')
        
        #return chans, velo, spec/pix
        return fig, ax
    
    def map_plot(self, data, contour_path=None, cont_levels=None, figsize=(10,10), cmap='jet', scale=0.8, colorbar=True, 
           grid=True, clabel=None, fraction=0.05, location='top', **kwargs):
        '''
        This can plot any of the 2D maps produced in the Mapper class
        
        Parameters
        ----------
        data: array-like
            the data you wish to plot
        contour_path: string (optional)
            the path to the contours you wish to overlay on your plot
        cont_levels: tuple (optional)
            the list of levels on your contour plot
        figsize: tuple (optional)
            the size of the figure
        cmap: string (optional)
            the colormap for the figure
        scale: float (optional)
            the scale of the figure, between 0 and 1
        colorbar: Boolean (optional)
            if True, a colorbar will appear with the figure
        grid: Boolean (optional)
            if True, a grid will appear on the figure
        clabel: string (optional)
            the label for the colorbar
        fraction: float (optional)
            the fraction of the figure the colorbar takes up
        location: string (optional)
            the location of the colorbar
        
        Retruns
        -------
        fig, ax: matplotlib attributes
            the fig, ax for the plot you created
        
        '''
        fig, ax = plt.subplots(figsize=figsize, subplot_kw={'projection':self.wcs_cube})
        ax.set_xlabel('Right Ascension (J2000) [hms]', fontsize=13)
        ax.set_ylabel('Declination (J2000) [degrees]', fontsize=13)

        med = np.nanmedian(data)
        sigma = np.nanstd(data)
        low = med - scale*sigma
        high = med + scale*sigma
    
        for key, value in kwargs.items():
            print(key, value)
        if not kwargs:
            im = ax.imshow(data, origin='lower', cmap=cmap, vmin=low, vmax=high)
        else:
            im = ax.imshow(data, origin='lower', cmap=cmap, **kwargs)
        if colorbar == True:
            plt.colorbar(im, pad=0.06, fraction=fraction, location=location,label=clabel)
            
        ax.tick_params(direction='in', length=7, labelsize=12, color='white')
    
        if grid == True:
            ax.coords.grid(color='white', ls='dotted', linewidth=1.3)
        
        if contour_path != None:
            contours = fits.open(contour_path)
            contour = contours[0].data
            wcs_contour = WCS(contours[0].header)
            plt.contour(contour, transform=ax.get_transform(wcs_contour), levels=cont_levels, colors='black', alpha=1.0)

        return fig, ax


def Load_fits(loc, ext=0):
    '''
    A function that opens and reads a fits file. It requires one string input, 
    which is the path to the fits file, and an optional integer input, which can be used 
    to change the extension of the fits file. 
    
    Parameters
    ----------
    loc: string
        the location of your fits file that you want to open and read, which has
        been identified as a string.
    ext: int, optional
        the extension, which has been identified as an integer, the user wishes to set 
        when they are opening the fits file.
    
    Returns
    -------
    header, data: fits header, numpy array
        This function returns a tuple one as your fits header and the other as a numpy array of your data. 
    '''
    with fits.open(loc) as hdu:
        header = hdu[ext].header
        data = hdu[ext].data
    return header, data 



def Interactive(base_path, cont_path):
    '''
    generates an interactive plot of the contours and vmin & vmax

    Parameters
    ----------
    base_path: string
        path to the 2D fits image, on which you wish to overlay the contours
    cont_path: string
        path to the 2D fits image that will be used to produce contours 

    '''
    header, data = load_fits(base_path) ##user input

    header_cont, data_cont = load_fits(cont_path)
    wcs_contour = WCS(header_cont)


    fig, ax = plt.subplots(figsize = (6,6))
    ax.imshow(data, origin = 'lower', cmap = 'gray')
    ax.set_title('Object')
    ax.set_xlabel('Pixels X')
    ax.set_ylabel('Pixels Y')

    contour_axis = plt.gca()

    divider = make_axes_locatable(ax)
    aspect = 25 #aspect ratio of the slider (tall and skinny)
    pad_fraction = 0.5 # padding
    width = axes_size.AxesY(ax, aspect = 1./aspect) # this will be the width of each slider
    pad = axes_size.Fraction(pad_fraction, width) # this is the true padding in fig coords

    axmax = divider.append_axes("bottom", size = width, pad = pad + 4*width) # add axes to the bottom
    axmin= divider.append_axes("bottom", size = width, pad = pad + 2.5*width) #add axes to the bottom
    smax = Slider(axmax, 'Cont Max', 0, 3*np.nanmean(data_cont), 50, color = 'blue') 
    smin = Slider(axmin, 'Cont Min', 0, 0.33*np.nanmean(data_cont), 0, color = 'blue')


    mean_ax = divider.append_axes("right", size = width, pad = pad) # add axes to the right
    scale_ax = divider.append_axes("right", size = width, pad = pad + 2*width) #add axes to the right
    mean = Slider(ax = mean_ax, 
              label = 'Mean', 
              valmin = 0.4*np.mean(data), #im
              valmax = 1.5*np.mean(data), #im
              valinit = np.mean(data), #im
              valstep = 0.5,
              color = 'red',
              orientation = 'vertical')

    scale = 0.5
    scale = Slider(ax = scale_ax, 
               label ='Scale',
               valmin = .33*(np.std(data)), #im
               valmax = 3*(np.std(data)), #im
               valinit = np.std(data), 
               valstep = 0.5,
               color = 'red',
               orientation = 'vertical')

    ct = contour_axis.contour(data_cont, levels = np.linspace(0.33*np.nanmean(data_cont),3*np.nanmean(data_cont), 10), colors='blue')
    contour_axis.clabel(ct, fontsize = 10, inline = 1, fmt = '%0.2f')

    def update(val):
        contour_axis.clear() 
        contour_axis.imshow(data, origin='lower', cmap='gray')
        contour_axis.set_title('Object')
        contour_axis.set_xlabel('Pixels X')
        contour_axis.set_ylabel('Pixels Y')
        image = contour_axis.get_images()[0]
        CS = contour_axis.contour(data_cont, levels=np.linspace(smin.val, smax.val, 10), colors='blue')
        contour_axis.clabel(CS, fontsize=10, inline=1, fmt = '%0.2f')
        vmin = mean.val - scale.val
        vmax = mean.val + scale.val
        image.set_clim(vmin, vmax)
        fig.canvas.draw()

    smax.on_changed(update)
    smin.on_changed(update)
    mean.on_changed(update)
    scale.on_changed(update)

    plt.show()






