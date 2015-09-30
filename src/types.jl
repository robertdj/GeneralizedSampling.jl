type Freq2wave1D
	# Sampling
	samples::Vector{Float64}
	weights::Union{Bool, Vector{Float64}}

	# Reconstruction
	wave::String
	column1::Vector{Complex{Float64}}
	# TODO: J is redundant; remove?
	J::Int

	# Multiplication with FFT: T*x = diag*NFFT(x)
	diag::Vector{Complex{Float64}}
	FFT::NFFT.NFFTPlan

	#Freq2wave1D() = new()
end


@doc """
	freq2wave(samples::Vector, wave::String, J::Int; B=0)

Make change of basis for switching between frequency responses and wavelets.

- `samples` are the sampling locations.
- `wave` is the name of the wavelet; see documentation for possibilities.
- `J` is the scale of the wavelet transform to reconstruct.
- `B` is the bandwidth of the samples; only needed if `samples` are non-uniform.

If the samples *are* uniform, `weights` in the output is `false`.

If the samples are *not* uniform, `weights` contains the weights as a vector and `basis` and `diag` are scaled with `sqrt(weights)`.
"""->
function freq2wave(samples::Vector, wave::String, J::Int; B::Float64=0.0)
	M = length(samples)
	# TODO: Warning if J is too big

	# Evaluate the first column in change of basis matrix
	# TODO: Parse strings such as "db4"
	func = string("Four", wave, "Scaling")
	column1 = eval(parse(func))( samples, J )

	if isuniform(samples)
		W = false
	else
		if B <= 0.0
			error("Samples are not uniform; supply bandwidth")
		end

		W = weights(samples, B)
		mu = complex( sqrt(W) )
		had!(column1, mu)
	end

	# NFFTPlans: Frequencies must be in the torus [-1/2, 1/2)
	N = 2^J
	xi = scale(samples, 1/N)
	frac!(xi)
	p = NFFTPlan(xi, N)
	diag = column1 .* cis(-pi*samples)

	return Freq2wave1D(samples, W, wave, column1, J, diag, p)
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
	(size(T,1), size(T,2))
end

function Base.size(T::Freq2wave1D, d::Int)
	if d == 1
		length(T.samples)
	elseif d == 2
		2^wscale(T)
	else
		error("Dimension does not exist")
	end
end

@doc """
	wscale(T::Freq2wave1D)

Return the scale of the wavelet coefficients.
"""->
function wscale(T::Freq2wave1D)
	T.J
end


@doc """
	mul!(T::Freq2wave1D, x::Vector, y::Vector)
	
Replace `y` with `T*x`.
"""->
function mul!(T::Freq2wave1D, x::Vector{Complex{Float64}}, y::Vector{Complex{Float64}})
	@assert size(T,1) == length(y)
	@assert size(T,2) == length(x)

	# TODO: nfft! requires that y contains only zeros. Change in NFFT package?
	fill!(y, 0.0+0.0*im)

	NFFT.nfft!(T.FFT, x, y)
	had!(y, T.diag)
end

@doc """
	mulT!(T::Freq2wave1D, v::Vector, z::Vector)
	
Replace `z` with `T'*v`.
"""->
function mulT!(T::Freq2wave1D, v::Vector{Complex{Float64}}, z::Vector{Complex{Float64}})
	@assert size(T,1) == length(v)
	@assert size(T,2) == length(z)

	D = conj(T.diag)
	had!(D, v)
	# TODO: As in mul!
	fill!(z, 0.0+0.0*im)
	NFFT.nfft_adjoint!(T.FFT, D, z)
end


@doc """
	*(T::Freq2wave, x::vector)

Compute `T*x`.
"""->
function Base.(:(*))(T::Freq2wave1D, x::Vector)
	y = Array(Complex{Float64}, size(T,1))
	mul!(T, complex(x), y)
	return y
end

@doc """
	Ac_mul_B(T::Freq2wave1D, x::vector)
	'*(T::Freq2wave1D, x::vector)

Compute `T'*x`.
"""->
function Base.Ac_mul_B(T::Freq2wave1D, x::Vector{Complex{Float64}})
	y = Array(Complex{Float64}, size(T,2))
	mulT!(T, complex(x), y)
	return y
end


@doc """
	collect(Freq2wave1D)

Return the full change of basis matrix.
"""->
function Base.collect(T::Freq2wave1D)
	# TODO: Check if the matrix fits in memory

	# Fourier matrix
	J = wscale(T)
	xk = T.samples*[0:2^J-1;]'
	scale!(xk, -2*pi*2.0^(-J))
	F = cis(xk)

	broadcast!(*, F, T.column1, F)
end

