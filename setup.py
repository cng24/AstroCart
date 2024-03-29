#setup.py
import setuptools

setuptools.setup(
    name="AstroCart",
    version="0.1",
    author="Caitlin Gainey",
    author_email="caitlin.gainey@yale.edu",
    description="A map generating and plotting application for CO isotopologue data, specifically 3D cubes",
    packages=setuptools.find_packages(include=['astrocart','astrocart/*']),
    python_requires='>=3',
    install_requires=["numpy", "matplotlib", "astropy", "spectral_cube", "reproject"]
)
