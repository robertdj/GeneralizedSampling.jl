module GeneralizedSampling

using Compat
import Compat: view

using Base.Cartesian

using NFFT

using IntervalWavelets
import IntervalWavelets: isuniform, van_moment

import VoronoiCells: density, voronoiarea
import Wavelets: wavelet, WT
import ArrayViews: flatten_view, reshape_view
import Base: A_mul_B!, Ac_mul_B!

export
	# Types
	CoB,
	Freq2Wave,
	Freq2BoundaryWave,
	Freq2NoBoundaryWave,

	# Fourier transforms
	FourScalingFunc,
	FourHaarScaling,
	FourDaubScaling,

	# Linear equation solvers
	cgnr,

	# Special CoB functions
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
const SMALL_EPS = eps()
const LARGE_EPS= sqrt(eps())
const ComplexOne = one(Complex{Float64})
const ComplexZero = zero(Complex{Float64})
const twoπ = 2.0*pi

include("Freq2Wave/Types.jl")
include("Misc.jl")
include("FourierTransforms.jl")
include("NFFT.jl")
include("CGNR.jl")

end # module
