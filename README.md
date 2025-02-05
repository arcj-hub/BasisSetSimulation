# Description
This project simulates semi-LASER basis sets using shaped refocusing pulses, creating the .basis files necessary for conventional fitting software. This simulation helps in the accurate analysis of spectroscopic data.

It utilizes functions from [FID-A](https://www.opensourceimaging.org/project/fid-a-advanced-processing-and-simulation-of-mr-spectroscopy/) and [Osprey](https://schorschinho.github.io/osprey/).

# Installation
To clone the repository, run the following command:

`git clone https://github.com/arcj-hub/BasisSetSimulation.git`

`cd BasisSetSimulation`

Required Dependencies: FID-A

# Usage
To simulate a basis set, run the following script after setting all the required input parameters.

`sLASER_makebasiset.m`

Note: You will need the refocusing pulse file specific to your vendor added to the path.

# Output

The output is a .basis file compatible with conventional fitting software. The pdf file includes simulated spectra based on the specified parameters and pulse sequences.


# Reference

This project utilizes functions from [FID-A](https://www.opensourceimaging.org/project/fid-a-advanced-processing-and-simulation-of-mr-spectroscopy/) and [Osprey](https://schorschinho.github.io/osprey/). If you find this work helpful or use any part of it in your research or publications, please consider citing it { xx NeuroLibre link } 


