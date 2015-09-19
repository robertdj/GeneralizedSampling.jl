module GeneralizedSampling

# TODO: qmf no longer available from Wavelets
#import Wavelets: qmf, wavelet
import Distributions: Categorical, sampler, rand
using NFFT
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
	density,
	frac,
	frac!,
	mul!,
	mulT!

end # module
