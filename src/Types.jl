@doc """
`CoB` is an abstract **c**hange **o**f **b**asis super-type 
"""->
abstract CoB


# ------------------------------------------------------------
# Frequency to wavelets

@doc """
`Freq2Wave` is a change of basis type between frequency samples and wavelets. 

There are sub types for wavelets with and without boundary correction.

The `weights` entry is a `Nullable` type and

- If the samples *are* uniform, `weights` is `Null`.
- If the samples are *not* uniform, `weights` contains the weights as a `Nullable` vector and `diag` are scaled with `weights`.
"""->
abstract Freq2Wave{D} <: CoB

macro common_freq2wave()
	esc(quote
		# Sampling
		internal::Array{Complex{Float64}, D}
		weights::Nullable{Vector{Complex{Float64}}}

		# Reconstruction
		J::Int
		wavename::AbstractString

		# In 1D: T = [left diag*NFFT right]
		# In 2D, diagonals for each dimension is needed
		diag::Array{Complex{Float64}, D}
		NFFT::NFFT.NFFTPlan{D,Float64}
	end)
end

# Uniform samples, no boundary correction
immutable Freq2NoBoundaryWave{D} <: Freq2Wave{D}
	@common_freq2wave()
end

macro bound(D)
	if D == 1
		left::DenseMatrix{Complex{Float64}}
		right::DenseMatrix{Complex{Float64}}
	else
		left::DenseArray{Complex{Float64}, 3}
		right::DenseArray{Complex{Float64}, 3}
	end
end

# Uniform samples, boundary correction
immutable Freq2BoundaryWave{D} <: Freq2Wave{D}
	@common_freq2wave()
	left::DenseArray{Complex{Float64}}
	right::DenseArray{Complex{Float64}}
	#= if D == 1 =#
	#= 	left::DenseMatrix{Complex{Float64}} =#
	#= 	right::DenseMatrix{Complex{Float64}} =#
	#= else =#
	#= 	left::DenseArray{Complex{Float64}, 3} =#
	#= 	right::DenseArray{Complex{Float64}, 3} =#
	#= end =#

	# TODO: Include arrays for temporary results in multiplication like NFFT.tmpVec?
	innery::Vector{Complex{Float64}}
end

function Freq2BoundaryWave{D}(internal::Array{Complex{Float64},D}, weights, J, wavename, diag, NFFT, left, right)
	innery = Array{Complex{Float64}}( size(internal,1) )

	Freq2BoundaryWave( internal, weights, J, wavename, diag, NFFT, left, right, innery )
end


# ------------------------------------------------------------
# Load methods for types

include("CoB/common.jl")
include("CoB/freq2wave.jl")

