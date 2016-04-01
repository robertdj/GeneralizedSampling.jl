module GeneralizedSampling

using ArrayViews
using NFFT
using WaveletPlot

import Deldir: deldir, voronoiarea
import Wavelets: wavelet, WT
#= import ArrayViews: view, reshape_view =#

const sqrt2 = sqrt(2)

export
	# Types
	CoB,
	Freq2Wave,
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
	grid,
	wavefilter


include("types.jl")
include("misc.jl")
include("FourierTransforms.jl")
#= include("CGNR.jl") =#

end # module
