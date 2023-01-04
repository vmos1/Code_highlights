# Fit the antipodal 4 pt function to the expected form of a tower of exponentials
The goal is to fit the antipodal four-point function to the expect form.
For the test data given here, which is for the free theory, the form is a series of integer exponents. 

These exponential fits are quite challenging and unstable.


This folder is self-contained with all the required data. All notebooks can be run directly after cloning the repository.

## Description

The Table below describes the location of codes and data.

| Name | Description |
| --- | ---|
| [cGAN_3d_pytorch/code/main.py](https://github.com/vmos1/Code_highlights/blob/main/3_cond_GANs_cosmology/cGAN_3d_pytorch/code/main.py) | Notebook to fit a single 4pt function |
|[2_combined_fit_4pt_fcn.ipynb]() | Notebook to perform a combined fit of 4pt functions with different 'l's |
| |Location of the free theory data |



## Conda environment 
This code requires the packages `gvar` and `lsqfit` that can be installed using pip.
