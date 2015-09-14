module GeneralizedSampling

# TODO: qmf no longer available from Wavelets
#import Wavelets: qmf, wavelet
import Distributions: Categorical, sampler, rand
using NFFT

# package code goes here
include("types.jl")
include("misc.jl")
include("FourierTransforms.jl")
include("CoB.jl")
include("Kaczmarz.jl")

export
	# Types
	Freq2wave1D,
	freq2wave,
	H,

	# Fourier transforms
	freq2Haar,
	FourHaarScaling,
	FourHaarWavelet,
	FourDaubScaling,

	# Linear equation solvers
	REK,

	# misc
	had!,
	isuniform,
	weights,
	density,
	frac,
	frac!

end # module
