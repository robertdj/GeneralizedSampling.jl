module GeneralizedSampling

using ArrayViews
using NFFT
using RCall # Used for Voronoi computations

import Wavelets: wavelet

export
	# Types
	Freq2wave1D,
	Freq2wave2D,
	freq2wave,

	# Fourier transforms
	FourHaarScaling,
	FourHaarWavelet,
	FourDaubScaling,

	# Linear equation solvers
	cgnr,

	# Special multiplication
	mul!,
	mulT!,
	collect,

	# misc
	had!,
	isuniform,
	weights,
	wscale,
	wsize,
	density,
	frac,
	frac!,
	grid


include("types.jl")
include("misc.jl")
include("FourierTransforms.jl")
include("CGNR.jl")

end # module
