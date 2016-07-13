GeneralizedSampling
===================

[![Build Status](https://travis-ci.org/robertdj/GeneralizedSampling.jl.svg?branch=master)](https://travis-ci.org/robertdj/GeneralizedSampling.jl)

A Julia package implementing the generalized sampling framework by Anders Hansen, Ben Adcock et al. between Fourier measurements and wavelets.


## Installation

The package is not yet registered, so install using

```julia
Pkg.clone("https://github.com/robertdj/GeneralizedSampling.jl")
```

Note that one of the necessary package [VoronoiCells](https://github.com/robertdj/VoronoiCells.jl) currently needs special attention.
Follow the instructions on *VoronoiCells'* the above link.


## Resources

- **Documentation**: <http://generalizedsamplingjl.readthedocs.org/en/latest/>


## Examples

The example file `examples/brain.jl` uses precomputed data that is *not* included in the package as it takes up about 21 MB.
These files can be downloaded from <http://people.math.aau.dk/~robert/software/brain.zip>.
Unzip the three `jld` files in the `examples` folder.

