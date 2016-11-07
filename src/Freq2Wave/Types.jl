@doc """
`CoB` is an abstract **c**hange **o**f **b**asis super-type 
"""->
abstract CoB


# ------------------------------------------------------------
# "N"FFT plan for equidistant points 

type FFTPlan{D}
    forwardFFT::Base.DFT.FFTW.cFFTWPlan{Complex{Float64},-1,true,D}
    backwardFFT::Base.DFT.FFTW.cFFTWPlan{Complex{Float64},1,true,D}
	pre_phaseshift::Array{Complex{Float64}, D}
	post_phaseshift::Vector{Complex{Float64}}
	x::Matrix
	M::NTuple{D, Int64}
    N::NTuple{D, Int64}
	q::NTuple{D, Int64}
	FFTvec::Array{Complex{Float64}, D}
end

function Base.show(io::IO, P::FFTPlan)
    print(io, "FFTPlan with ", P.M, " sampling points for ", P.N, " array")
end

FFTPlan(samples::Vector, J::Integer, N::Integer) = FFTPlan( samples', J, (N,) )

function FFTPlan{D}(samples::AbstractMatrix, J::Integer, N::NTuple{D, Int})
    size(samples, 1) == D || throw(DimensionMismatch())

    M = Array{Int64}(D)
    for d in 1:D
        M[d] = length( unique(samples[d, :]) )
    end
    prod(M) == size(samples, 2) || throw(AssertionError())

    myeps = samples[D,2] - samples[D,1]
    inv_eps = 1 / myeps
    if isapprox(inv_eps, round(inv_eps))
        inv_eps = round(Int, inv_eps)
    else
        throw(DomainError())
    end

    pre = Array{Complex{Float64}}(N)
    dilation = 2^J
    for idx in CartesianRange(size(pre))
        pre[idx] = cis( pi*myeps/dilation*dot([idx.I...]-1, M) )
    end

    xi = samples / dilation
    log_post = [N...]'* xi
    post = cis( pi*vec(log_post) )

    q = div(M, dilation*inv_eps) + 1
    FFTvec = Array{Complex{Float64}}( tuple(q*dilation*inv_eps...) )

    fP = plan_fft!(FFTvec)#; flags=FFTW.PATIENT)
    bP = plan_bfft!(FFTvec)#; flags=FFTW.PATIENT)

    FFTPlan{D}(fP, bP, pre, post, xi, (M...), N, (q...), FFTvec)
end


# ------------------------------------------------------------
# Frequency to wavelets

@doc """
`Freq2Wave` is a change of basis type between frequency samples and wavelets. 

There are sub types for wavelets with and without boundary correction.
To initialize a `Freq2Wave` type, run

	Freq2Wave(samples, wavename::String, J::Int, B; ...)

- `samples` are the sampling locations as a vector for 1D and a matrix with 2 columns for 2D.
- `wave` is the name of the wavelet; see documentation for possibilities.
- `J` is the scale of the wavelet transform to reconstruct.
- `B` is the bandwidth of the samples; only needed if `samples` are non-uniform.
- Optional arguments (if needed) are passed to the functions computing Fourier transforms.
"""->
abstract Freq2Wave <: CoB

abstract Freq2Wave1D <: Freq2Wave
abstract Freq2Wave2D <: Freq2Wave
# The 2D Freq2BoundaryWave's need entries the 1D's do not, so don't
# use the dimension as a TypePar.

# No boundary correction
immutable Freq2NoBoundaryWave1D{P} <: Freq2Wave1D
	# Sampling
	internal::Vector{Complex{Float64}}
	weights::Nullable{ Vector{Complex{Float64}} }

	# Reconstruction
	J::Int64
	wavename::AbstractString

	# Multiplication
	NFFT::P

	tmpMulVec::Vector{Complex{Float64}}
end

immutable Freq2NoBoundaryWave2D{P} <: Freq2Wave2D
	internal::Dict{ Symbol, Vector{Complex{Float64}} }
	weights::Nullable{ Vector{Complex{Float64}} }

	J::Int64
	wavename::AbstractString

	NFFT::P

	tmpMulVec::Vector{Complex{Float64}}
end

# Boundary correction
immutable Freq2BoundaryWave1D{P} <: Freq2Wave1D
	internal::Vector{Complex{Float64}}
	weights::Nullable{ Vector{Complex{Float64}} }

	J::Int64
	wavename::AbstractString

	NFFT::P

	left::Matrix{Complex{Float64}}
	right::Matrix{Complex{Float64}}

	tmpMulVec::Vector{Complex{Float64}}
end

immutable Freq2BoundaryWave2D{P} <: Freq2Wave2D
	internal::Dict{ Symbol, Vector{Complex{Float64}} }
	weights::Nullable{ Vector{Complex{Float64}} }

	J::Int64
	wavename::AbstractString

	NFFT::P
	NFFTx::NFFT.NFFTPlan{2,1,Float64}
	NFFTy::NFFT.NFFTPlan{2,2,Float64}

	left::Dict{ Symbol, Matrix{Complex{Float64}} }
	right::Dict{ Symbol, Matrix{Complex{Float64}} }

	tmpMulVec::Matrix{Complex{Float64}}
	tmpMulVecT::Matrix{Complex{Float64}} # Serves as the transpose of tmpMulVec
	tmpMulcVec::Vector{Complex{Float64}}
	weigthedVec::Vector{Complex{Float64}}
end


# ------------------------------------------------------------------------
# Constructors

function Freq2NoBoundaryWave1D(internal, weights, J, wavename, NFFT)
	tmpMulVec = Array{Complex{Float64}}( NFFT.M )
	Freq2NoBoundaryWave1D( internal, weights, J, wavename, NFFT, tmpMulVec )
end

function Freq2NoBoundaryWave2D(internal, weights, J, wavename, NFFT)
    tmpMulVec = Array{Complex{Float64}}( prod(NFFT.M) )
	Freq2NoBoundaryWave2D( internal, weights, J, wavename, NFFT, tmpMulVec )
end

function Freq2BoundaryWave1D(internal, weights, J, wavename, NFFT, left, right)
	tmpMulVec = Array{Complex{Float64}}( NFFT.M )
	Freq2BoundaryWave1D( internal, weights, J, wavename, NFFT, left, right, tmpMulVec )
end

function Freq2BoundaryWave2D(internal, weights, J, wavename, NFFT, left, right)
	tmpMulVec = similar(left[:x])
	tmpMulVecT = tmpMulVec.'

    tmpMulcVec = Array{Complex{Float64}}( prod(NFFT.M) )
	weigthedVec = similar(tmpMulcVec)

	vm = van_moment(wavename)
	NFFTx = NFFTPlan( frac(NFFT.x[1,:]), 1, (NFFT.N[1], vm) )
	NFFTy = NFFTPlan( frac(NFFT.x[2,:]), 2, (vm, NFFT.N[2]) )

	Freq2BoundaryWave2D( internal, weights, J, wavename, NFFT,
	NFFTx, NFFTy, left, right, tmpMulVec, tmpMulVecT, tmpMulcVec, weigthedVec )
end


# ------------------------------------------------------------------------
# Load methods

include("freq2wave.jl")

