@doc """
`CoB` is an abstract *c*hange *o*f *b*asis super-type 
"""->
abstract CoB


# ------------------------------------------------------------
# Frequency to wavelets

@doc """
`Freq2wave` is a change of basis type between frequency samples and wavelets. 
"""->
immutable Freq2wave{D} <: CoB
	# Sampling
	samples::Array{Float64,D}
	weights::Nullable{Vector{Float64}}

	# Reconstruction
	wave::AbstractString
	column1::Vector{Complex{Float64}}
	# TODO: J is redundant; remove?
	J::Int

	# Multiplication with FFT: T*x = diag*NFFT(x)
	diag::Vector{Complex{Float64}}
	FFT::NFFT.NFFTPlan{D,Float64}
end


# ------------------------------------------------------------
# Load functions for types

include("CoB/common.jl")
include("CoB/freq2wave.jl")

