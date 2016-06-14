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
abstract Freq2Wave <: CoB

abstract Freq2Wave1D <: Freq2Wave
abstract Freq2Wave2D <: Freq2Wave
# The 2D Freq2BoundaryWave's need a lot that the 1D's do not.
# Therefore, I don't use the dimension as a TypePar.

macro common_freq2wave(D)
	esc(quote
		# Sampling
		internal::Array{Complex{Float64}, D}
		weights::Nullable{ Vector{Complex{Float64}} }

		# Reconstruction
		J::Int64
		wavename::AbstractString

		# Multiplication
		diag::Vector{Complex{Float64}}
		NFFT::NFFT.NFFTPlan{D,Float64}
	end)
end

# No boundary correction
immutable Freq2NoBoundaryWave1D <: Freq2Wave1D
	# Sampling
	internal::Vector{Complex{Float64}}
	weights::Nullable{ Vector{Complex{Float64}} }

	# Reconstruction
	J::Int64
	wavename::AbstractString

	# Multiplication
	diag::Vector{Complex{Float64}}
	NFFT::NFFT.NFFTPlan{1,Float64}
end

immutable Freq2NoBoundaryWave2D <: Freq2Wave2D
	internal::Matrix{Complex{Float64}}
	weights::Nullable{ Vector{Complex{Float64}} }

	J::Int64
	wavename::AbstractString

	diag::Matrix{Complex{Float64}}
	NFFT::NFFT.NFFTPlan{2,Float64}
end

# Boundary correction
immutable Freq2BoundaryWave1D <: Freq2Wave1D
	internal::Vector{Complex{Float64}}
	weights::Nullable{Vector{Complex{Float64}}}

	J::Int64
	wavename::AbstractString

	diag::Vector{Complex{Float64}}
	NFFT::NFFT.NFFTPlan{1,Float64}

	# TODO: Make type stable
	left::Matrix{Complex{Float64}}
	right::Matrix{Complex{Float64}}

	innery::Vector{Complex{Float64}}
end

immutable Freq2BoundaryWave2D <: Freq2Wave2D
	internal::Matrix{Complex{Float64}}
	weights::Nullable{Vector{Complex{Float64}}}

	J::Int64
	wavename::AbstractString

	diag::Matrix{Complex{Float64}}
	NFFT::NFFT.NFFTPlan{2,Float64}
	NFFTx::NFFT.NFFTPlan{1,Float64}
	NFFTy::NFFT.NFFTPlan{1,Float64}

	# TODO: Make type stable
	left::Matrix{Any}
	right::Matrix{Any}

	innery::Vector{Complex{Float64}}
end


# ------------------------------------------------------------------------
# Constructors

function Freq2BoundaryWave1D(internal::Vector{Complex{Float64}}, weights, J, wavename, diag, NFFT, left, right)
	innery = Array{Complex{Float64}}( length(internal) )
	Freq2BoundaryWave1D( internal, weights, J, wavename, diag, NFFT, left, right, innery )
end

function Freq2BoundaryWave2D(internal::Matrix{Complex{Float64}}, weights, J, wavename, diag, NFFT, left, right)
	#= size(internal,1) == 2 || throw(DimensionMismatch()) =#
	innery = Array{Complex{Float64}}( size(internal,1) )

	# TODO: Different for uniform samples
	NFFTx = NFFTPlan( NFFT.x[1,:], (NFFT.N[1],) )
	NFFTy = NFFTPlan( NFFT.x[2,:], (NFFT.N[2],) )

	Freq2BoundaryWave2D( internal, weights, J, wavename, diag, NFFT, NFFTx, NFFTy, left, right, innery )
end


# ------------------------------------------------------------------------
# Load methods

include("CoB/freq2wave.jl")
#= include("Freq2Wave/freq2wave.jl") =#

