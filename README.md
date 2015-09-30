GeneralizedSampling
===================

[![Build Status](https://travis-ci.org/robertdj/GeneralizedSampling.jl.svg?branch=master)](https://travis-ci.org/robertdj/GeneralizedSampling.jl)

This package implements the generalized sampling framework by Anders Hansen, Ben Adcock et al. between Fourier measurements and wavelets.

The theory behind generalized sampling can be found in e.g.:

- Ben Adcock and Anders C. Hansen, "A Generalized Sampling Theorem for Stable Reconstructions in Arbitrary Bases", Journal of Fourier Analysis and Applications, 2012, vol. 18, p. 685--716.
[DOI](https://dx.doi.org/10.1007/s00041-012-9221-x) [arXiv](http://arxiv.org/abs/1007.1852)
- Ben Adcock, Anders C. Hansen and Milana Gataric, "On stable reconstructions from univariate nonuniform Fourier measurements", SIAM Journal on Imaging Sciences, 2014, vol. 7, pp. 690--1723.
[DOI](https://dx.doi.org/10.1137/130943431) [arXiv](http://arxiv.org/abs/1310.7820)
- Ben Adcock, Anders C. Hansen and Milana Gataric, "Weighted frames of exponentials and stable recovery of multidimensional functions from nonuniform Fourier samples", 2015.
[arXiv](http://arxiv.org/abs/1405.3111)
-  Milana Gataric, Clarice Poon, "A practical guide to the recovery of wavelet coefficients from Fourier measurements", 2015,
[arXiv](http://arxiv.org/abs/1505.05308)


## Dependencies

The following Julia packages are required for GeneralizedSampling:

- [NFFT](https://github.com/tknopp/NFFT.jl)
- [Wavelets](https://github.com/JuliaDSP/Wavelets.jl)

In order to use the Kaczmarz algorithm for solving least squares problems you also need

- [Distributions](https://github.com/JuliaStats/Distributions.jl)
- [ArrayViews](https://github.com/JuliaLang/ArrayViews.jl)

