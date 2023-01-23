# Code_highlights
This is a repository highlighting codes from a few of my current projects. A brief description is provided below. Further details and illustrations are provided in README files in each respective folder.


## 1. Lattice Multigrid solver in 2D for Laplace and Wilson operators
Most lattice gauge theory efforts to study strongly coupled field theories such as Quantum Chromodynamics (QCD) suffer from the problem of `critical slowing down` i.e. conventional solvers slown dramatically as one explores the physics region of small lattice spacing. Lattice multigrid methods have shown remarkable success in ameliorating this problem for different formulations such as Wilson, Staggered and Domain-wall fermions. 

Here we explore the potential of improved Multigrid methods with the `non-telescoping` method, involving using multiple blocking schemes of the system at lower levels of the system. In this work, we discuss the performance for the Laplace and Wilson operators.

## 2. Multi-exponential fit of the 4-point function in $ R x S^2 $

The study of Conformal field theories (CFTs) is challenging on the lattice, due to their property of scale invariance. The quantum finite element method [arxiv:2006.15636](https://arxiv.org/abs/2006.15636) has shown great promise in studying CFTs on the lattice. 

Here, we are study the $\phi^4$ theory using the lattice field theory approach on curved manifolds (such as $R \times S^2$) with the intention of exploring the behavior near the conformal fixed point.

The goal is fit the anti-podal four-point function to a series of exponentials to extract the coeffecients and exponents.
Using these we can look at the behavior as we approach the continuum limit.
For simiplicity, we use the data from the free theory, without any interactions.

## 3. Conditional GANs to generate mass distributions in the early universe
We develop conditional Generative Adversarial Neural networks (cGANs) to produce maps of the universe, conditioned on the cosmological paramters $\sigma$.
The models are trained with images obtained from N-body simulations. Once trained, these ML models can very quickly generate novel images belonging to the distribution. Conditional GANs offer the opportunity to obtain novel images in regions of parameter space not explored by N-body simulations.


## Setting up the environment
A number of files in this repository use interactive widgets in the jupyterlab environment. These codes also make use of custom packages such as [`lsqfit`](https://pypi.org/project/lsqfit/). The easiest way to run these codes is to setup a custom conda environment. The folder [conda_env](https://github.com/vmos1/Code_highlights/tree/main/conda_env) has files to set this up and [README](https://github.com/vmos1/Code_highlights/tree/main/conda_env/README.md) gives more details on setting this up.
