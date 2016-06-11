@doc """
`CoB` is an abstract **c**hange **o**f **b**asis super-type 
"""->
abstract CoB


# ------------------------------------------------------------
# Frequency to wavelets

@doc """
`Freq2Wave` is a change of basis type between frequency samples and wavelets. 

There are sub types for wavelets with and without boundary correction.
To initialize a `Freq2Wave` type, run

	Freq2Wave(samples, wavename::String, J::Int, B; ...)

- `samples` are the sampling locations as a vector for 1D and M-by-2 matrix for 2D.
- `wave` is the name of the wavelet; see documentation for possibilities.
- `J` is the scale of the wavelet transform to reconstruct.
- `B` is the bandwidth of the samples; only needed if `samples` are non-uniform.
- Optional arguments (if needed) are passed to the functions computing Fourier transforms.
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
		diag::Matrix{Complex{Float64}}
		NFFT::NFFT.NFFTPlan{D,Float64}
		NFFTx::NFFT.NFFTPlan{1,Float64}
		NFFTy::NFFT.NFFTPlan{1,Float64}
	end)
end

# Uniform samples, no boundary correction
immutable Freq2NoBoundaryWave{D} <: Freq2Wave{D}
	@common_freq2wave()
end

# Uniform samples, boundary correction
immutable Freq2BoundaryWave{D} <: Freq2Wave{D}
	@common_freq2wave()
	# TODO: Make type stable
	left::Matrix{Any}
	right::Matrix{Any}

	# TODO: Include arrays for temporary results in multiplication like NFFT.tmpVec?
	innery::Vector{Complex{Float64}}
end

function Freq2BoundaryWave{D}(internal::Array{Complex{Float64},D}, weights, J, wavename, diag, NFFT, left, right)
	innery = Array{Complex{Float64}}( size(internal,1) )

	if D == 2
		NFFTx = NFFTPlan( NFFT.x[1,:], (NFFT.N[1],) )
		NFFTy = NFFTPlan( NFFT.x[2,:], (NFFT.N[2],) )
	else
		NFFTx = NFFT
		NFFTy = NFFT
	end

	Freq2BoundaryWave( internal, weights, J, wavename, diag, NFFT, NFFTx, NFFTy, left, right, innery )
end


# ------------------------------------------------------------
# Load methods for types

#= include("CoB/common.jl") =#
include("CoB/freq2wave.jl")

