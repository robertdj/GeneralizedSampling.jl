GeneralizedSampling
===================

[![Build Status](https://travis-ci.org/robertdj/GeneralizedSampling.jl.svg?branch=master)](https://travis-ci.org/robertdj/GeneralizedSampling.jl)
[![AppVeyor](https://ci.appveyor.com/api/projects/status/github/robertdj/GeneralizedSampling.jl?branch=master&svg=true)](https://ci.appveyor.com/project/robertdj/generalizedsampling-jl)
[![codecov.io](https://codecov.io/github/robertdj/GeneralizedSampling.jl/coverage.svg?branch=master)](https://codecov.io/github/robertdj/GeneralizedSampling.jl?branch=master)
[![Documentation Status](https://readthedocs.org/projects/generalizedsamplingjl/badge/?version=latest)](http://generalizedsamplingjl.readthedocs.io/en/latest/?badge=latest)

A Julia package implementing the generalized sampling framework by Anders Hansen, Ben Adcock et al. from Fourier measurements to Daubechies wavelets.


## Installation

In Julia, run

```julia
Pkg.add("GeneralizedSampling")
```


## Resources

- **Documentation**: <http://generalizedsamplingjl.readthedocs.org/en/latest/>


## Examples

The example file `examples/brain.jl` uses precomputed data from this [MRI phantom Matlab package](http://bigwww.epfl.ch/algorithms/mriphantom) and the data is *not* included in the package as it takes up about 21 MB.

The data is available for download in the native Julia format `jld` from <http://people.math.aau.dk/~robert/software/brain.zip>.
Unzip the three `jld` files in the `examples` folder.

