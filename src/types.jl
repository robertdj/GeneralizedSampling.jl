@doc """
`CoB` is an abstract *c*hange *o*f *b*asis super-type 
"""->
abstract CoB


# ------------------------------------------------------------
# Frequency to wavelets

@doc """
`Freq2wave1D` is a change of basis type between 1D frequency samples and 1D wavelets. 
"""->
immutable Freq2wave1D <: CoB
	# Sampling
	samples::Vector{Float64}
	weights::Nullable{Vector{Float64}}

	# Reconstruction
	wave::String
	column1::Vector{Complex{Float64}}
	# TODO: J is redundant; remove?
	J::Int

	# Multiplication with FFT: T*x = diag*NFFT(x)
	diag::Vector{Complex{Float64}}
	FFT::NFFT.NFFTPlan{1,Float64}
end

@doc """
`Freq2wave2D` is a change of basis type between 2D frequency samples and 2D wavelets. 
"""->
immutable Freq2wave2D <: CoB
	# TODO: If samples *are* uniform, this should be *two* matrices
	# Sampling
	samples::Matrix{Float64}
	weights::Nullable{Vector{Float64}}

	# Reconstruction
	wave::String
	column1::Vector{Complex{Float64}}
	# TODO: J is redundant; remove?
	J::Int

	# Multiplication with FFT: T*x = diag*NFFT(x)
	diag::Vector{Complex{Float64}}
	FFT::NFFT.NFFTPlan{2,Float64}
end

@doc """
`Freq2wave` is the union of `Freq2wave1D` and `Freq2wave2D`.
"""->
Freq2wave = Union{Freq2wave1D, Freq2wave2D}


# ------------------------------------------------------------
# Load functions for types

include("CoB/common.jl")
include("CoB/freq2wave.jl")

