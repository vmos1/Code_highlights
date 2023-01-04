# Code_highlights
Highlighting codes from a few projects


## 1. Lattice Multigrid solver in 2D for Laplace and Wilson operators
Most lattice gauge theory efforts to study strongly coupled field theories such as Quantum Chromodynamics (QCD) suffer from the problem of `critical slowing down` i.e. conventional solvers slown dramatically as one explores the physics region of small lattice spacing. Lattice multigrid methods have shown remarkable success in ameliorating this problem for different formulations such as Wilson, Staggered and Domain-wall fermions. Here we explore the potential of improved Multigrid methods with the `non-telescoping` method, involving using multiple blocking schemes of the system at lower levels of the system. In this work, we discuss the performance for the Laplace and Wilson operators.

## 2. Multi-exponential fit of the 4-point function in R x S^2
We are studying conformal field theory.
(computed on  $ R \times S^2 $ manifold )
The goal is fit the anti-podal four-point function to a series of exponentials to extract the coeffecients and exponents.
Using these we can look at the behavior 
For simiplicity, we use the data from the free thoery, without any interactions.

## 3. Conditional GANs to generate mass distributions in the early universe
We develop conditional Generative Adversarial Neural networks (cGANs) to produce maps of the universe, conditioned on the cosmological paramters $\sigma$.
The models are trained with images obtained from N-body simulations. The goal is to use cGANs to obtain novel images in regions of parameter space not explored by N-body simulations.


## Setting up the environment
A number of files in this repository use interactive widgets in the jupyterlab environment. These codes also make use of custom packages such as [`lsqfit`](https://pypi.org/project/lsqfit/). The easiest way to run these codes is to setup a custom conda environment. The folder [conda_env](https://github.com/vmos1/Code_highlights/tree/main/conda_env) has files to set this up and [README](https://github.com/vmos1/Code_highlights/tree/main/conda_env/README.md) gives more details on setting this up.
