#mapper code
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.wcs import WCS
from spectral_cube import SpectralCube
import astropy.units as u
from astropy.utils import data
import os
from reproject import reproject_interp

def Qrot(T):
    '''
    The rotational partition function, which will be called in the column density function
    
    Parameters
    ----------
    T: float
        the excitation temperatures
    
    Returns
    -------
    zeta: float
        the rotational partition function evaluated at the input excitation temperature
    '''
    zeta=0.0
    leng = 100000
    ar = np.arange(leng)
    
    for j in range (0,leng):
        f1 = 2*j +1
        Ej = j*(j+1)*h*B
        f2 = np.exp(-Ej/(k*T))
        zeta = zeta + (f1*f2)
        
    return zeta

def Ntot(freq, Tex, Aul, gup, Elow, area): 
    '''
    Column Density in Optically Thin Regime
    
    Parameters
    ----------
    freq: float
        the frequency of your transition
    Tex: float
        the excitation temperature
    Aul: float
        Einstein coefficient
    gup: float
        the statistical weight (or degeneracy) of the upper state of the given transition
    Elow: float
        the energy of the lower state of the given transition
    area: float
        the channel width * intensity at the given pixel
    
    Retruns
    -------
    Ntot*1e5: float
        Column Density in cm^2
    '''
    #s= time.time()
    #print(f'Starting F1 at {s}')
    
    #f1 = 8.*np.pi*freq**3/(Aul*gup*c**3)
    #f2 = (h*freq)/k
    #f3 = 1./((1./(np.exp(f2/Tex)-1))-(1./(np.exp(f2/Tbkg))))
    #f4 = 1.-np.exp(-f2/Tex)
    #f5 = np.exp(Elow/Tex)
    
    #s5 = time.time() - s 
    #const = (f1/f2)*f3*f5*area/f4 #Qrot(Tex)
    #print(f'all steps before Qrot {s5}')

    #Ntot = (f1/f2)*f3*f5*Qrot(Tex)*area/f4
    
    #s6 = time.time() - s5
    #print(f'Qrot took {s6}')

    f1 = 8.*np.pi*freq**3/(Aul*gup*c**3)
    f2 = (h*freq)/k
    f3 = 1./((1./(np.exp(f2/Tex)-1))-(1./(np.exp(f2/Tbkg))))
    f4 = 1.-np.exp(-f2/Tex)
    f5 = np.exp(Elow/Tex)
    Ntot = (f1/f2)*f3*f5*Qrot(Tex)*area/f4
        
    return Ntot*1e5

def Ntot_t(freq, Tex, Aul, gup, Elow,tau):
    '''
    Column Density in non-Optically Thin Regime
    
    Parameters
    ----------
    freq: float
        the frequency of your transition
    Tex: float
        the excitation temperature
    Aul: float
        the Einstein coefficient
    gup: float
        the statistical weight (or degeneracy) of the upper state of the given transition
    Elow: float
        the energy of the lower state of the given transition
    tau: float
        the optical depth
    
    Retruns
    -------
    Ntot*1e5: float
        Column Density in cm^2
    '''
    f1 = 8.*np.pi*freq**3/(Aul*gup*c**3)
    f2 = (h*freq)/k
    f4 = 1.-np.exp(-f2/Tex)
    f5 = np.exp(Elow/Tex)
    Ntot = f1*f5*Qrot(Tex)*tau/f4

    return Ntot*1e5
    
def tau(freq,Tex,Tmb):
    '''
    The optical depth
    
    Parameters
    ----------
    freq: float
        the frequency of the given transition
    Tex: float
        the excitation temperature of a given element in the dataset
    Tmb: float
        the main beam brightness temperature
    
    Returns
    -------
    -np.log(f4): float
        the optical line depth for an element of the dataset, given the input parameters 
    
    '''
    f1 = Tmb*1.38e-23/(freq*6.62e-34)
    f2 = (6.62e-34*freq)/1.38e-23
    f3 = 1./((1./(np.exp(f2/Tex)-1))-(1./(np.exp(f2/2.73))))
    f4 = 1.+(f1*f3)
    return -np.log(f4) 

def RMS_Map(cube, chan_low, chan_high):
    '''
    Generating the RMS Map from a subcube of channels that have no signal 
    
    Parameters
    ----------
    cube: 3D array-like 
        fits data 
    chan_low: int
        the lowest channel without signal you want for your RMS calculation
    chan_high: int
        the highest channel without signal you want for your RMS calculation
    
    Returns
    -------
    RMS: array-like
        the RMS data for each pixel 
    
    '''
    Dsquare = np.square(cube[0].data)
    P1 = np.sum(Dsquare[chan_low:chan_high,:,:],axis=0)

    rms = np.sqrt(P1/(chan_high - chan_low))
    RMS = rms[:,:]
    return RMS
 
def Pixel_Regrid(path_template, path_toRegrid, One_Dim3=False):
    '''
    Regridding two data sets into the same shape and RA/DEC
    
    Parameters
    ----------
    path_template: string
        the path to the data set you would like to make your template
        
    path_toRegrid: string
        the path to the data set you would like to regrid to the template data
        
    One_Dim3: Boolean (optional)
        Sets whether or not your template data is 3D, while your regridding data is 2D
        If these are both 3D, you can run One_Dim3 = False
        
    Returns
    -------
    regridded: array-like
        The regridded data set
    '''
    template = fits.open(path_template)
    toRegrid = fits.open(path_toRegrid)
    if One_Dim3 == True:
        cube1_data = fits.open(path_template) #MUST MAKE 2D, so do moment map
        cubesigma = SpectralCube.read(cube1_data)  # Initiate a SpectralCube
        cubesigma = cubesigma.with_spectral_unit(u.km / u.s)
        moment_0 = cubesigma.with_spectral_unit(u.km/u.s).moment(order=0)

        template_new = moment_0
        regridded, _ = reproject_interp(toRegrid, template_new.header)
        return regridded
    else:
        regridded, _ = reproject_interp(toRegrid, template[0].header)
        return regridded
    
class Mappers():
    def __init__(self, path_trans1, path_trans2):
        '''
        Inputs the two transitions of the given CO isotopologue (i.e. 13CO(1-0) and 13CO(2-1)) as well as their
        channel width and shape and saves them as class attributes
        
        Parameters
        ----------
        path_trans1: string
            the path to your first transition
        path_trans2: string
            the path to your second transition
        
        '''
        self.CO1 = fits.open(path_trans1)
        self.CO2 = fits.open(path_trans2)
        self.dv = np.abs(self.CO2[0].header['CDELT3'])/1e3
        
        self.Vsh, self.Ysh, self.Xsh = self.CO2[0].shape

    def run_RMS(self, chan_low, chan_high):
        '''
        Runs the RMS maps of the two transitions. At this point, the data is assumed to be 
        gridded to the same shape, so the chan_low and chan_high will match across data sets
        
        Parameters
        ----------
        chanlow_l: int
            the lowest channel below your signal from which you would like to take the rms
        chanlow_h: int
            the highest channel below your signal from which you would like to take the rms
        chanhigh_l: int
            the lowest channel above your signal from which you would like to take the rms
        chanhigh_h: int
            the highest channel above your signal from which you would like to take the rms
            
        Returns
        -------
        self.RMS2, self.RMS1: array-like
            this returns the RMS data of the two input transitions
        '''
        self.RMS2 = RMS_Map(self.CO2, chan_low, chan_high)
        self.RMS1 = RMS_Map(self.CO1, chan_low, chan_high)
        return self.RMS1, self.RMS2

    def Tex_thin(self, RMS1, RMS2, Tex_low = 3., Tex_high = 63., steps = 150):
        '''
        Runs the excitation temperature in the optically thin limit, given your two transitions
        
        Parameters:
        -----------
        Tex_low: float (optional)
            The lowest possible excitation temperature to try
        Tex_high: float (optional)
            The highest possible excitation temperature to try
        steps: int (optional)
            The number of steps you want between your Tex_low and Tex_high
            
        Returns
        -------
        self.tempMap_thin: array-like
            The excitation temperature data in the optically thin limit
        '''
        Tex = np.linspace(Tex_low, Tex_high, num = steps ,endpoint=True)
        
        results = np.zeros((len(Tex),2))
        tempMap = np.zeros((self.Vsh, self.Ysh, self.Xsh))
        NCO1 = Ntot(freq10, Tex, Aul10, gup10, Elow10,1.)
        NCO2 = Ntot(freq21, Tex, Aul21, gup21, Elow21,1.)

        tempMapMean = 0.0
        weight = 0.0
        for xx in range(0, self.Xsh):
            for yy in range(0, self.Ysh):
                for v in range(0, self.Vsh):
                    if self.CO1[0].data[v,yy,xx]>=3*RMS1[yy,xx] and self.CO2[0].data[v,yy,xx]>=3*RMS2[yy,xx]:
                        for i in range(0,len(Tex)):
                            N1 = NCO1[i]*self.CO1[0].data[v,yy,xx]*np.abs(self.dv)
                            N2 = NCO2[i]*self.CO2[0].data[v,yy,xx]*np.abs(self.dv)
                            results[i]=[np.abs((N2/N1)-1.), Tex[i]]
                        tempMap[v,yy,xx]= Tex[np.argmin(results[:,0])]
                    else:
                        tempMap[v,yy,xx]='nan'

        for xx in range(0, self.Xsh):
            for yy in range(0, self.Ysh):
                if np.any(np.isnan(tempMap[:,yy,xx])==False):
                    Tmean = np.nanmean(tempMap[:,yy,xx])
                    for v in range(0, self.Vsh):
                        if np.isnan(tempMap[v,yy,xx])==True:
                            tempMap[v,yy,xx]=Tmean

        meanTex = np.nanmean(tempMap)
        for xx in range(0, self.Xsh):
            for yy in range(0, self.Ysh):
                if np.all(np.isnan(tempMap[:,yy,xx])==True):
                    for v in range(0, self.Vsh):
                        tempMap[v,yy,xx]=meanTex
        self.tempMap_thin = tempMap
        return self.tempMap_thin
    
    def Tex_thick(self, RMS1, RMS2, Tex_low = 3., Tex_high = 63., steps = 150):
        '''        
        Runs the excitation temperature in the non-optically thin limit, given your two transitions
        
        Parameters:
        -----------
        Tex_low: float (optional)
            The lowest possible excitation temperature to try
        Tex_high: float (optional)
            The highest possible excitation temperature to try
        steps: int (optional)
            The number of steps you want between your Tex_low and Tex_high
            
        Returns
        -------
        self.tempMap_thick: array-like
            The excitation temperature data in the non-optically thin limit
        '''
        tempMap = np.zeros((self.Vsh, self.Ysh, self.Xsh))
        
        Tex = np.linspace(Tex_low, Tex_high, num=steps, endpoint=True)
        results = np.zeros((len(Tex), 4))
        NCO1 = Ntot_t(freq10, Tex, Aul10, gup10, Elow10,1.)
        NCO2 = Ntot_t(freq21, Tex, Aul21, gup21, Elow21,1.)

        tempMapMean = 0.0
        weight = 0.0
        for xx in range(0, self.Xsh):
            for yy in range(0, self.Ysh):
                for v in range(0, self.Vsh):
                    if self.CO1[0].data[v,yy,xx]>=3*RMS1[yy,xx] and self.CO2[0].data[v,yy,xx]>=3*RMS2[yy,xx]:
                        for i in range(0,len(Tex)):
                            tau1 = tau(freq10,Tex[i],self.CO1[0].data[v,yy,xx])
                            tau2 = tau(freq21,Tex[i],self.CO2[0].data[v,yy,xx])
                            N1 = NCO1[i]*tau1
                            N2 = NCO2[i]*tau2
                            results[i]=[np.abs((N2/N1)-1.), Tex[i], tau1,tau2]
                        tempMap[v,yy,xx]= Tex[np.argmin(results[:,0])]
                    else:
                        tempMap[v,yy,xx]='nan'

        for xx in range(0, self.Xsh):
            for yy in range(0, self.Ysh):
                if np.any(np.isnan(tempMap[:,yy,xx])==False):
                    Tmean = np.nanmean(tempMap[:,yy,xx])
                    for v in range(0, self.Vsh):
                        if np.isnan(tempMap[v,yy,xx])==True:
                            tempMap[v,yy,xx]=Tmean
    
        meanTex = np.nanmean(tempMap)

        for xx in range(0, self.Xsh):
            for yy in range(0, self.Ysh):
                if np.all(np.isnan(tempMap[:,yy,xx])==True):
                    for v in range(0, self.Vsh):
                        tempMap[v,yy,xx]=meanTex
                        
        self.tempMap_thick = tempMap
        return self.tempMap_thick

    def Col_Den_thin(self, TexMap_thin):
        '''
        Runs the column density in the optically thin limit, given the optically thin 
        excitation temperature map that should've been run previously.
        
        Parameters
        ----------
        TexMap_thin: array-like
            This is the map of excitation temperature in the optically thin limit. 
            This map can be produced in the previous steps.
            
        Returns
        -------
        self.Ncol_thin: array-like
            The column density data in the optically thin limit
        
        '''
        
        self.Ncol_thin = Ntot(freq21, TexMap_thin, Aul21, gup21, Elow21,1.)*self.CO2[0].data*self.dv
        return self.Ncol_thin
    
    def Col_Den_thick(self, TexMap_thick):
        '''
        Runs the column density in the non-optically thin limit, given the non-optically thin 
        excitation temperature map that should've been run previously.
        
        Parameters
        ----------
        TexMap_thick: array-like
            This is the map of excitation temperature in the non-optically thin limit. 
            This map can be produced in the previous steps.
        
        Returns
        -------
        self.Ncol_thin: array-like
            The column density data in the non-optically thin limit
        '''
        self.Ncol_thick = Ntot_t(freq21, TexMap_thick, Aul21, gup21, Elow21,1.)*self.CO2[0].data*self.dv
        return self.Ncol_thick
    
    def Dim2_Map(self, tempMap, Ncol): 
        '''
        Generates 2D maps from cubes. For the tempMap, we produced a column density weighted map, 
        and for the Ncol, we produced a summed map for every position-position-velocity element.
        
        Parameters
        ----------
        TempMap: array-like
            This is the 3D cube of excitation temperature to be weighted and collapsed using the input Ncol.
        
        Ncol: array-like
            This is the 3D cube of column density, but the user can also put in any other 3D map they wish
            to sum into a 2D map.
        
        Returns
        -------
        self.TexMap2d, self.NcoMap2d: array-like, array-like
            The 2D maps of excitation temperature and column density (or any other characteristic 
            the user decided to make 2D)
        '''
        TexMap2d = np.zeros((self.Ysh, self.Xsh))
        NcoMap2d = np.zeros((self.Ysh, self.Xsh))
        for xx in range(0, self.Xsh):
            for yy in range(0, self.Ysh):
                TexMap2d[yy,xx] = np.nansum(tempMap[:,yy,xx]*Ncol[:,yy,xx])/np.nansum(Ncol[:,yy,xx])
                NcoMap2d[yy,xx] = np.nansum(Ncol[:,yy,xx])
        self.TexMap2d = TexMap2d
        self.NcoMap2d = NcoMap2d
        return self.TexMap2d, self.NcoMap2d
    
    def sigma_map(self, NCOMap2d, C18O=False):
        '''
        Makes the mass surface density map 
        
        Parameters
        ----------
        NCOMap2d: array-like
            the 2D column density map, generated in previous steps
        C18O: Boolean
            if False, the code will run for 13CO mass surface density
            if True, the code will run for C18O mass surface density
            
        Returns
        -------
        self.sigma = array-like
            the mass surface density data 
        '''
        for xx in range(0,self.Xsh):
            for yy in range(0,self.Ysh):
                if NCOMap2d[yy,xx] <= 0:
                    NCOMap2d[yy,xx] = 'nan'
        if C18O == False:
            self.sigma = (1.24*10**-2) * NCOMap2d / 10**16 ##13CO 
            return self.sigma
        else: 
            self.sigma = (7.652*10**-2) * NCOMap2d / (10**16) ##C18O
            return self.sigma
    
    def CODep_Map(self, sigma, known_sigma_path):
        '''
        Makes the CO Depletion Map
        
        Parameters
        ----------
        sigma: array-like
            the mass surface density map that the user created above 
        known_sigma_path: string
            the path to the mass surface density that you would like to compare our calculated 
            mass surface density to
        
        Returns
        -------
        self.CO_Dep: array-like
            the CO Depletion Map 
        
        '''
        known = fits.open(known_sigma_path)
        known_sigma = known[0].data

        for xx in range(0,self.Xsh):
            for yy in range(0,self.Ysh):
                if known_sigma[yy,xx] <= 0:
                    known_sigma[yy,xx] = 'nan'
    
        self.CO_Dep = known_sigma / sigma
        return self.CO_Dep
    

    
