module GeneralizedSampling

import Wavelets: qmf, wavelet

# package code goes here
include("CoB.jl")
include("FourierTransforms.jl")

export
	freq2Haar,
	FourHaarScaling,
	FourHaarWavelet,
	FourDaubScaling

end # module
