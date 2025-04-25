# Description
This project simulates semi-LASER basis sets using shaped refocusing pulses, producing .basis files compatible with conventional MRS fitting software. Simulated spectra are also generated and exported in PDF format. Core simulation routines leverage functions from FID-A and Osprey for accurate spectroscopic analysis.

It utilizes functions from [FID-A](https://www.opensourceimaging.org/project/fid-a-advanced-processing-and-simulation-of-mr-spectroscopy/) and [Osprey](https://schorschinho.github.io/osprey/).

# Modifications
The sLASER_makebaseset.m script, this code has been extensively modified by Niklaus Zoelch and Jessica Archibald (2020–2025) to add:

  -A 0-ppm reference peak

  -An automated loop over all specified metabolites

  -Saving of raw data (.raw), figure snapshots (.png), and MATLAB workspace files (.mat)

  -Output of both a compiled .pdf report and a .basis file, using modified Osprey routines (Georg Oeltzschner, Johns Hopkins University, 2019)

# Installation
To clone the repository, run the following command:

`git clone https://github.com/arcj-hub/BasisSetSimulation.git`

`cd BasisSetSimulation`

Required Dependencies: FID-A

# Directory Structure

BasisSetSimulation/
├── dependencies/   # Needed functions
├── my_basis/       # Output folder: simulated .basis files
├── my_mets/        # Input folder: custom metabolite spin-systems
├── my_pulse/       # Input folder: vendor-specific RF pulse files
├── run_generate_my_mets.m   # Script to add new metabolite spin systems
├── sLASER_makebaseset.m   # Main simulation script
├── README.md       # This documentation
└── LICENSE         # Project license

my_basis/
Stores all generated .basis files

my_mets/
Place any new metabolite spin-system definitions here (automatically from run_generate_my_mets.m) 

my_pulse/
Add the shaped refocusing pulse files provided by your scanner vendor

# Usage
To simulate a basis set, run the following script after setting all the required input parameters.

`sLASER_makebasiset.m`

Note: You will need the refocusing pulse file specific to your vendor added to the path.

# Output

Review outputs in my_basis/:

.basis files for each metabolite

.pdf report visualizing simulated spectra

.raw, .png, and .mat files for further inspection

# Reference

This project utilizes functions from [FID-A](https://www.opensourceimaging.org/project/fid-a-advanced-processing-and-simulation-of-mr-spectroscopy/) and [Osprey](https://schorschinho.github.io/osprey/). If you find this work helpful or use any part of it in your research or publications, please consider citing it { xx NeuroLibre link } 


