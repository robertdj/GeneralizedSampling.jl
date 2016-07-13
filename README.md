GeneralizedSampling
===================

[![Build Status](https://travis-ci.org/robertdj/GeneralizedSampling.jl.svg?branch=master)](https://travis-ci.org/robertdj/GeneralizedSampling.jl)
[![codecov.io](https://codecov.io/github/robertdj/GeneralizedSampling.jl/coverage.svg?branch=master)](https://codecov.io/github/robertdj/GeneralizedSampling.jl?branch=master)
[![Documentation Status](https://readthedocs.org/projects/generalizedsamplingjl/badge/?version=latest)](http://generalizedsamplingjl.readthedocs.io/en/latest/?badge=latest)

A Julia package implementing the generalized sampling framework by Anders Hansen, Ben Adcock et al. from Fourier measurements to Daubechies wavelets.


## Installation

The package is not yet registered, so install using

```julia
Pkg.clone("https://github.com/robertdj/GeneralizedSampling.jl")
```

Note that one of the necessary packages [VoronoiCells](https://github.com/robertdj/VoronoiCells.jl) currently needs special attention;
follow the instructions on the link.


## Resources

- **Documentation**: <http://generalizedsamplingjl.readthedocs.org/en/latest/>


## Examples

The example file `examples/brain.jl` uses precomputed data from this [MRI phantom Matlab package](http://bigwww.epfl.ch/algorithms/mriphantom) and the data is *not* included in the package as it takes up about 21 MB.

The data is available for download in the native Julia format `jld` from <http://people.math.aau.dk/~robert/software/brain.zip>.
Unzip the three `jld` files in the `examples` folder.

