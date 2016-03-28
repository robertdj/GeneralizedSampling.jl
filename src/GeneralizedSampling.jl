module GeneralizedSampling

using ArrayViews
using NFFT

import Deldir: deldir, voronoiarea
import Wavelets: wavelet, WT

export
	# Types
	CoB,
	Freq2wave,
	freq2wave,
	BoundaryFilter,
	ScalingFilters,

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
include("Filters.jl")
include("FourierTransforms.jl")
include("CGNR.jl")

end # module
