module GeneralizedSampling

using NFFT
using WaveletPlot
import WaveletPlot: isuniform, van_moment

import Deldir: deldir, voronoiarea
import Wavelets: wavelet, WT
import ArrayViews: reshape_view

export
	# Types
	CoB,
	Freq2Wave,
	Freq2BoundaryWave,
	Freq2NoBoundaryWave,

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
	had!,
	hadc!,
	yphad!,
	isuniform,
	grid,
	weights,
	density,
	van_moment,
	isdaubechies,
	ishaar,
	hasboundary,
	dim,
	split,
	wscale,
	wsize,
	frac,
	frac!

const sqrt2 = sqrt(2)
# TODO: Rename: SMALL_EPS, LARGE_EPS
const SMALL_PREC = eps()
const LARGE_PREC = sqrt(eps())
const ComplexOne = one(Complex{Float64})

include("Types.jl")
include("Misc.jl")
include("FourierTransforms.jl")
include("CGNR.jl")

end # module
