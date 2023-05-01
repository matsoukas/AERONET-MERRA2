# AERONET_SSA
Comparison of Single Scattering Albedo (SSA) and Aerosol Absorption Optical Depth (AAOD) between all AERONET stations and the MERRA-2 reanalysis.

AERONET product: "All points” Almucantar V3 Inversion Data, Level 2.0
MERRA-2 product: tavg1_2d_aer_Nx
The comparision resulted in the manuscript with DOI: 10.2139/ssrn.4385670

Files

LICENCE: A GNU GPL v.3 licence

README.md: The current file

mylib.py: A collection of various functions used in the other python programs

read_AERONET.py: Code reading relevant data from all the files downloaded from https://aeronet.gsfc.nasa.gov/new_web/download_all_v3_inversions.html,     Almucantar Products (Except Phase Functions) 2.0, All points

MERRA550.py: Code reading MERRA2 550 nm data downloaded from https://disc.gsfc.nasa.gov/datasets/M2T1NXAER_5.12.4/summary, Subset/Get Data

graphs.py: Code producing the figures in DOI: 10.2139/ssrn.4385670

environment.yml: Environment file for the python libaries. On linux, one can run "> conda env create -f environment.yml" which will create the same environment as the one on whcich we perfomed the analysis.

sites_and_continents.csv: A list of AERONET stations used, their coordinates and their continents

We also used Peter A. Rochford (2016) SkillMetrics: A Python package for calculating the skill of model predictions against observations, http://github.com/PeterRochford/SkillMetrics. This skill_metrics directory was copied in the local directory with all the above files.
