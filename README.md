# AstroCart
Astrophysical Cartography (AstroCart) is a package that can produce RMS, excitation temperature, column density, 2D, mass surface density, and CO Depletion maps. It also has a spectrum and map plotting tool, which can be put into interactive mode to easily change the vmin and vmax of the base image and overlay varying levels of contours from another image. All images should be in fits format. 

## Installation
- Git Clone https://github.com/cng24/AstroCart.git
- Navigate to the AstroCart directory
- pip install -e . 

The package should now be importable as 
```ruby
import astrocart
```

## Mapping
The following will introduce what each main feature of the mapper code does 

### Pixel Regridding
Regrids a fits image to a template fits image—both of which are input by the user. In addition to regridding to the same RA/DEC, this function will also regrid dimensions of the input data to the template data.
- N.B. the template image should be the one with the worst resolution of the two data sets

run with...
```ruby
temp = 'path_to_template_image.fits'
reg = 'path_to_image_tobe_regridded.fits'

#One_Dim3=True if template data is 3D, while your regridding data is 2D
regridded = astrocart.Pixel_Regrid(path_template=temp, path_toRegrid=reg, One_Dim3=True) 
```

### Mappers
The class that runs RMS, excitation temperature, column density, mass surface density, and CO depletion maps. 
initialize with...
```'ruby
path1 = 'path_to_first_data.fits'
path2 = 'path_to_second_data.fits '
pipe = astrocart.Mappers(path1, path2)
```

#### RMS Map 
Produces a 2D RMS map, given the 3D cube of data
run with...
```ruby
RMS1, RMS2 = pipe.run_RMS(10, 50) #between channels 10 and 50
```

#### Excitation Temperature Map 
Produces a 3D cube of excitation temperature, taking values only over 3rms at each position-position-velocity (or frequency) element. This excitation temperature mapper requires the user to have two different transitions of the same CO isotopologue, e.g. 13CO(J=1-0) and 13CO(J=2-1). 
- N.B. the user is required to have the following variables defined before running excitation temperature and column density map. In the file constants.py, the user can input these values. The default values are for 13CO(J=1-0) and 13CO(J=2-1). These constants can be found at [Spatalogue](https://splatalogue.online//) 
   - B = rotational constant for the CO isotopologue
   - h = planck's constant
   - k = boltzmann constant
   - c = speed of light
   - Tbkg = cosmic microwave background temperature
   - freq10 = frequency of first transition
   - Aul10  = einstein coefficient of first transition
   - gup10  = the statistical weight (or degeneracy) of the upper transition of the second transition
   - Elow10 = the energy of the lower state of the first transition
   - freq21 = frequency of second transition
   - Aul21  = einstein coefficient of second transition
   - gup21  = the statistical weight (or degeneracy) of the upper transition of the second transition
   - Elow21 = the energy of the lower state of the second transition

run with...
```ruby
tempMap_thin = pipe.Tex_thin(RMS1, RMS2)
```

#### Column Density Map
Produces a 3D cube of column density, which the user can choose to produce in the optically thin or non-optically thin regime.
- N.B. the column density calculations in the non-optically thin regime corrects for optical depth but does NOT correct for opacity.  

run with...
```ruby
ColDen_thin = pipe.Col_Den_thin(tempMap_thin)
```

#### 2D Maps
Produces a 2D map from a 3D cube. The 2D collapse of the 3D excitation temperature map is weighted with the column density map.
- N.B. to sum over the 3D data cube and produce a 2D map, put the data cube in for the Ncol argument
- N.B. to produce a weighted 2D map from two 3D cubes of data, put the data you want to collapse in the tempMap argument with the weight in the Ncol argument. 

run with...
```ruby
Tex2, Col2 = pipe.Dim2_Map(tempMap_thin, ColDen_thin)
```

#### Mass Surface Density Map
Produces a mass surface density map for the input 2D column density map. 
- N.B. this currently works for only 13CO and C18O

run with...
```ruby
sig13CO = pipe.sigma_map(Col2, C18O=False) #if running 13CO
```

#### CO Depletion Map
Produces a CO Depletion map from a 2D mass surface density map (also referred to as a sigma map) produced with CO and a known mass surface density map produced through some other method e.g. a Herschel mass surface density map.

run with...
```ruby
ksig_path = 'path_to_knownsigma_map.fits'
CODep = pipe.CODep_Map(sig13CO, ksig_path)
```


## Plotting
The following will introduce what each main feature of the plotter code does 

### Plotter
The class that with all non-interactive plotting code. 
initialize with...
```ruby
path1 = 'path_to_first_data.fits' #either data one or two if they are regridded to the same RA/DEC and shape
plots = astrocart.Plotter(path1)
```

#### Spectrum
Produces a spectrum from a 3D cube of data with km/s on the x-axis
```ruby
fig, ax = plots.spectrum()
```

#### Map_Plot
Produces a plot of any 2D map with the option to add a colorbar, contours, and more. 
```ruby
fig, ax = plots.map_plot(Tex2, vmin=5, vmax=10) #no contours
```

if the user would like to add contours, then the add...
```ruby
cont = 'pathto_image_thatwillbe_contours.fits'

levels=np.linspace(60, 120, 5)
fig, ax = plots.map_plot(Tex2, contour_path=cont, cont_levels=levels, vmin=5, vmax=10) #no contours
```

### Interactive
Produces an interactive map of 2D data with a contour plot overlay. The user is able to vary the minimum and maximum value on their contour levels using 'Cont Max' and 'Cont Min' as well as the vmin and vmax of their base image using 'Mean' and 'Scale'.
- N.B. The interactive component of this map will not run in Jupyter Lab

```ruby
base = 'path_to_base_image.fits'
cont = 'pathto_image_thatwillbe_contours.fits'

astrocart.Interactive(base, cont)
```





