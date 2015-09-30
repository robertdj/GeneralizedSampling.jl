module GeneralizedSampling

using NFFT
import Distributions: Categorical, sampler, rand
import Wavelets: wavelet
import ArrayViews: view

# package code goes here
include("types.jl")
include("misc.jl")
include("FourierTransforms.jl")
include("Kaczmarz.jl")
include("CGNR.jl")

export
	# Types
	Freq2wave1D,
	freq2wave,

	# Fourier transforms
	FourHaarScaling,
	FourHaarWavelet,
	FourDaubScaling,

	# Linear equation solvers
	REK,
	cgnr,

	# misc
	had!,
	isuniform,
	weights,
	wscale,
	density,
	frac,
	frac!,
	mul!,
	mulT!

end # module
