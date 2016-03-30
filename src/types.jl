@doc """
`CoB` is an abstract **c**hange **o**f **b**asis super-type 
"""->
abstract CoB


# ------------------------------------------------------------
# Frequency to wavelets

@doc """
`Freq2Wave` is a change of basis type between frequency samples and wavelets. 

There are sub types for wavelets with and without boundary correction.
"""->
abstract Freq2Wave <: CoB

macro common_freq2wave()
	esc(quote
		# Sampling
		# TODO: Are the parts of samples saved in NFFT sufficient for collect?
		samples::Array{Float64, D}
		weights::Nullable{Vector{Float64}}

		# Reconstruction
		J::Int
		wavename::AbstractString

		# In 1D: T = [left diag*NFFT right]
		diag::Vector{Complex{Float64}}
		NFFT::NFFT.NFFTPlan{D,Float64}
	end)
end

# Uniform samples, no boundary correction
immutable Freq2NoBoundaryWave{D} <: Freq2Wave
	@common_freq2wave()
end

# Uniform samples, boundary correction
immutable Freq2BoundaryWave{D} <: Freq2Wave
	@common_freq2wave()
	left::Matrix{Complex{Float64}}
	right::Matrix{Complex{Float64}}
end

function isuniform(T::Freq2Wave)
	!isnull(T.weights)
end

function withboundary(T::Freq2Wave)
	isdefined(T, :left)
end


# ------------------------------------------------------------
# Load functions for types

include("CoB/common.jl")
include("CoB/freq2wave.jl")

