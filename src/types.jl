type Freq2wave1D
	# Sampling
	samples::Vector{Float64}
	weights::Union{Bool, Vector{Float64}}

	# Reconstruction
	wave::String
	column::Vector{Complex{Float64}}
	# TODO: J is redundant; remove?
	J::Int

	# Multiplication with FFT: T*x = diag*NFFT(x)
	diag::Vector{Complex{Float64}}
	FFT::NFFT.NFFTPlan

	#Freq2wave1D() = new()
end


@doc """
	freq2wave(samples::Vector, wave::String, J::Int)

Make change of basis for switching between frequency responses and wavelets.

- `samples` are the sampling locations.
- `wave` is the name of the wavelet; see documentation for possibilities.
- `J` is the scale of the wavelet transform to reconstruct.
"""->
function freq2wave(samples::Vector, wave::String, J::Int)
	M = length(samples)
	# TODO: Warning if J is too big

	# NFFTPlans: Frequencies must be in the torus [-1/2, 1/2)
	N = 2^J
	xi = scale(samples, 1/N)
	xi_frac = frac(xi)
	p = NFFTPlan(xi_frac, N)

	# Evaluate the first column in change of basis matrix
	# TODO: Parse strings such as "db4"
	func = string("Four", wave, "Scaling")
	basis = eval(parse(func))( samples, J, 0 )
	diag = basis .* cis(-pi*N*xi_frac)

	if isuniform(samples)
		weights = false
	else
		# TODO: Bandwidth?
		weights = weights(samples, B)
	end

	return Freq2wave1D(samples, weights, wave, basis, J, diag, p)
end


function Base.show(io::IO, T::Freq2wave1D)
	println(io, "1D change of basis matrix")

	typeof(T.weights) == Bool ?  U = " " : U = " non-"

	M, N = size(T)
	println(io, "From: ", M, U, "uniform frequency samples")
	print(io, "To: ", N, " ", T.wave, " wavelets")
end


# ------------------------------------------------------------
# Basic operations for Freq2wave1D

@doc """
	size(Freq2wave) -> (M,N)
	size(Freq2wave, d) -> M or N

`M` is the number of samples and `N` is the number of reconstruction functions.
"""->
function Base.size(T::Freq2wave1D)
	(length(T.samples), 2^T.J)
end

function Base.size(T::Freq2wave1D, d::Int)
	if d == 1
		length(T.samples)
	elseif d == 2
		2^T.J
	else
		error("Dimension does not exist")
	end
end


@doc """
Compute `Freq2wave * vector`
"""->
function Base.(:(*))(T::Freq2wave1D, x::Vector{Complex{Float64}})
	@assert size(T,2) == length(x)

	y = nfft(T.FFT, x)
	had!(y, T.diag)

	return y
end

function Base.(:(*)){T<:Number}(C::Freq2wave1D, x::Vector{T})
	z = complex(float(x))
	return C*z
end


# TODO: Find a better name for multiplication with the adjoint
@doc """
	H(T::Freq2wave1D, x)

Compute `T'*x`.
"""->
function H(T::Freq2wave1D, x::Vector{Complex{Float64}})
	@assert size(T,1) == length(x)

	D = conj(T.diag)
	y = had!(D, x)
	z = nfft_adjoint(T.FFT, y)
end

