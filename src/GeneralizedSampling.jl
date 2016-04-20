module GeneralizedSampling

using NFFT
using WaveletPlot

import Deldir: deldir, voronoiarea
import Wavelets: wavelet, WT
import ArrayViews: reshape_view

const sqrt2 = sqrt(2)
const SMALL_PREC = eps()
const LARGE_PREC = sqrt(eps())

export
	# Types
	CoB,
	Freq2Wave,

	# Fourier transforms
	FourHaarScaling,
	FourHaarWavelet,
	FourDaubScaling,

	# Linear equation solvers
	cgnr,

	# Special CoB functions
	freq2wave,
	collect,

	# misc
	isuniform,
	weights,
	wscale,
	wsize,
	density,
	frac,
	frac!,
	grid,


include("Types.jl")
include("Misc.jl")
include("FourierTransforms.jl")
include("CGNR.jl")

end # module
