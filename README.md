# stochastic-schroedinger-equation

## Description
This code numerically solves the linear stochastic Schroedinger equation with additive or multiplicative Gaussian noise. For the spatial discretisation, a spectral Galerkin approach is used, and for the temproal discretisation, the exponential Euler, implicit Euler, and Crank-Nicolson method are compared. The numerical convergence rates obtained correspond to the analytical convergence rates proven in 'Pathwise Uniform Convergence of Time Discretisation Schemes for SPDEs' by Katharina Klioba and Mark Veraar, see [arXiv 2303.00411](https://arxiv.org/abs/2303.00411).

## How to reproduce the data and plots
After cloning this repository, run the following programs in Matlab. Version R2022b or newer is advised.
- To reproduce the numerical convergence rates for the case of additive noise, run `schroedinger_spectral_additive.m`.
- To reproduce the numerical convergence rates for the case of multiplicative noise, run `schroedinger_spectral_multiplicative.m`.
- To recreate Figure 1 from the paper, run `plotpaperfigure.m`. This requires the files `ErrEXP_additive_samples100_M1024_dt25912.mat` and `ErrEXP_multiplicative_samples100_M1024_dt25912.mat` from this repository.

## Requirements
Installation of [Matlab](https://de.mathworks.com/products/matlab.html) version R2022b or newer is required.

## Authors and acknowledgment
Authors: Katharina Klioba (Hamburg University of Technology) and Mark Veraar (Delft University of Technology)

The second author is supported by the VICI subsidy VI.C.212.027 of the Netherlands Organisation for Scientific Research (NWO).

## License
GNU GENERAL PUBLIC LICENSE                       Version 3, 29 June 2007
