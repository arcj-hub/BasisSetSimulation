# Basis Set Simulation

This project simulates semi-LASER basis sets using shaped refocusing pulses, producing .basis files compatible with conventional MRS fitting software. Simulated spectra are also generated and exported in PDF format. Core simulation routines leverage functions from FID-A and Osprey for accurate spectroscopic analysis.

Contact: Dr. Jessica Archibald ([jess.archibald9@proton.me](mailto:jess.archibald9@proton.me))

## Description

This code creates MRS metabolite basis sets based on the semi-LASER sequence and produces the .basis files necessary for conventional MRS spectral fitting software. The simulated basis sets are integral to the accurate analysis of MR spectroscopic data.

The code utilizes functions from [FID-A](https://github.com/CIC-methods/FID-A) and [Osprey](https://schorschinho.github.io/osprey/).

## Installation

To clone the repository, run the following command:

`git clone https://github.com/arcj-hub/BasisSetSimulation.git`

Required dependencies:

- [FID-A](https://github.com/CIC-methods/FID-A) (this is included in the repository)

## Directory Structure

### `dependencies/`

Stores required functions

### `my_basis/`

Stores all generated .basis files

### `my_mets/`

Place any new metabolite spin-system definitions here (automatically from run_generate_my_mets.m)

### `my_pulse/`

Add the shaped refocusing pulse files provided by your scanner vendor

### `sLASER_makebaseset.m`

Main simulation script

### `run_generate_my_mets.m`

Script to add new metabolite spin systems

## Usage

To simulate a basis set, run the following script after setting all the required input parameters.

`run_sLASER_makebasisset.m`

Note: You will need the refocusing pulse file specific to your vendor added to the path.

## Output

Review outputs in `my_basis/`:

.basis files for each metabolite

.pdf report visualizing simulated spectra

.raw, .png, and .mat files for further inspection

## Citing the Work

If you find this work helpful or use any part of it in your research or publications, please consider citing it { xx NeuroLibre link }.

## Acknowledgments

This project utilizes functions from [FID-A](https://github.com/CIC-methods/FID-A) and [Osprey](https://schorschinho.github.io/osprey/).
