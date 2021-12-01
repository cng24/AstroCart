#setup.py
import setuptools

setuptools.setup(
    name="AstroCart",
    version="0.1",
    author="Caitlin Gainey",
    author_email="caitlin.gainey@yale.edu",
    description="A map generating and plotting application for CO isotopologue data, specifically 3D cubes",
    #packages=["exampy","exampy/plotting"]
    packages=setuptools.find_packages(include=['AstroCart','AstroCart.*']),
    python_requires='>=3',
    install_requires=["numpy", "matplotlib.pyplot", "astropy.io", "astropy.wcs", 
    "spectral_cube", "astropy.units", "reproject.reproject_interp", "matplotlib.widgets.Slider", "mpl_toolkits.axes_grid1"]
)
