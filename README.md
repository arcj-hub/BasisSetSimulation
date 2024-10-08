# Basis set simulation

Contact: Dr. Jessica Archibald, Weill Cornell Medicine ([jea4025@med.cornell.edu](mailto:jea4025@med.cornell.edu))

## Description

This code creates MRS metabolite basis sets based on the semi-LASER sequence and produces the .basis files necessary for conventional MRS spectral fitting software. The simulated basis sets are integral to the accurate analysis of MR spectroscopic data.

The code utilizes functions from [FID-A](https://github.com/CIC-methods/FID-A) and [Osprey](https://schorschinho.github.io/osprey/).

## Installation

To clone the repository, run the following command:

`git clone https://github.com/arcj-hub/BasisSetSimulation.git`

Required dependencies:

- [FID-A](https://github.com/CIC-methods/FID-A)

## Usage

To simulate a basis set, run the following script after setting all the required input parameters.

`run_sLASER_makebasiset.m`

Note: You will need the refocusing pulse file specific to your vendor added to the path.

## Output

The output is a .basis file compatible with conventional fitting software. The PDFs include simulated spectra based on the specified parameters and pulse sequences.

## Acknowledgments

This project utilizes functions from [FID-A](https://github.com/CIC-methods/FID-A) and [Osprey](https://schorschinho.github.io/osprey/). If you find this work helpful or use any part of it in your research or publications, please consider citing it { xx NeuroLibre link }.
