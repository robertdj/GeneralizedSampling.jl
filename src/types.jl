#=
@doc """
	CoB

Change of basis
"""->
abstract CoB
Freq2wave <: CoB
=#

# ------------------------------------------------------------
# 1D change of basis matrix

immutable Freq2wave1D
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
	FFT::NFFT.NFFTPlan{1,Float64}

	#Freq2wave1D() = new()
end

@doc """
	freq2wave(samples, wave::String, J::Int; B=0)

Make change of basis for switching between frequency responses and wavelets.

- `samples` are the sampling locations.
- `wave` is the name of the wavelet; see documentation for possibilities.
- `J` is the scale of the wavelet transform to reconstruct.
- `B` is the bandwidth of the samples; only needed if `samples` are non-uniform.

If the samples *are* uniform, `weights` in the output is `false`.

If the samples are *not* uniform, `weights` contains the weights as a vector and `basis` and `diag` are scaled with `sqrt(weights)`.
"""->
function freq2wave(samples::Vector, wave::String, J::Int; B::Float64=0.0)
	# TODO: Warning if J is too big

	# Evaluate the first column in change of basis matrix
	# TODO: Parse strings such as "db4"
	func = string("Four", wave, "Scaling")
	column1 = eval(parse(func))( samples, J )

	# NFFTPlans: Frequencies must be in the torus [-1/2, 1/2)
	# TODO: Should window width m and oversampling factor sigma be changed for higher precision?
	N = 2^J
	xi = scale(samples, 1/N)
	frac!(xi)
	p = NFFTPlan(xi, N)
	diag = column1 .* cis(-pi*samples)

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
	mul!(T::Freq2wave, x::Vector, y::Vector)
	
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
	mulT!(T::Freq2wave, v::Vector, z::Vector)
	
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
	Ac_mul_B(T::Freq2wave, x::vector)
	'*(T::Freq2wave1D, x::vector)

Compute `T'*x`.
"""->
function Base.Ac_mul_B(T::Freq2wave1D, x::Vector{Complex{Float64}})
	y = Array(Complex{Float64}, size(T,2))
	mulT!(T, complex(x), y)
	return y
end


@doc """
	collect(Freq2wave)

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


# ------------------------------------------------------------
# 2D change of basis matrix

immutable Freq2wave2D
	# Sampling
	samples::Matrix{Float64}
	weights::Union{Bool, Vector{Float64}}

	# Reconstruction
	wave::String
	column1::Vector{Complex{Float64}}
	# TODO: J is redundant; remove?
	J::Int

	# Multiplication with FFT: T*x = diag*NFFT(x)
	diag::Vector{Complex{Float64}}
	FFT::NFFT.NFFTPlan{2,Float64}
end

Freq2wave = Union(Freq2wave1D, Freq2wave2D)

@doc """
	size(Freq2wave) -> (M,N)
	size(Freq2wave, d) -> M or N

`M` is the number of samples and `N` is the number of reconstruction functions.
"""->
function Base.size(T::Freq2wave)
	(size(T,1), size(T,2))
end

@doc """
	wscale(T::Freq2wave)

Return the scale of the wavelet coefficients.
"""->
function wscale(T::Freq2wave)
	T.J
end

function freq2wave(samples::Matrix, wave::String, J::Int; B::Float64=0.0)
	M, D = size(samples)
	@assert D == 2
	# TODO: Warning if J is too big

	# Evaluate the first column in change of basis matrix
	# TODO: Parse strings such as "db4"
	func = string("Four", wave, "Scaling")
	samplesx = view(samples, :, 1)
	samplesy = view(samples, :, 2)
	column1 = eval(parse(func))( samplesx, J )
	column1y = eval(parse(func))( samplesy, J )
	had!(column1, column1y)

	diag = column1 .* cis( -pi*(samplesx + samplesy) )

	# NFFTPlans: Frequencies must be in the torus [-1/2, 1/2)
	# TODO: This is for functions on the unit square. Should rectangles be allowed?
	N = 2^J
	xi = scale(samples, 1/N)
	frac!(xi)
	p = NFFTPlan(xi', (N,N))

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

	return Freq2wave2D(samples, W, wave, column1, J, diag, p)
end


function Base.show(io::IO, T::Freq2wave2D)
	println(io, "2D change of basis matrix")

	typeof(T.weights) == Bool ?  U = " " : U = " non-"

	M, N = size(T)
	N = Int(sqrt(N))
	println(io, "From: ", M, U, "uniform frequency samples")
	print(io, "To: ", N, "-by-", N, " ", T.wave, " wavelets")
end

# TODO: Should size(T,2) return total number or for one dimension?
function Base.size(T::Freq2wave2D, d::Int)
	if d == 1
		size(T.samples, 1)
	elseif d == 2
		4^wscale(T)
	else
		error("Dimension does not exist")
	end
end

function mul!(T::Freq2wave2D, x::Vector{Complex{Float64}}, y::Vector{Complex{Float64}})
	@assert size(T,1) == length(y)
	@assert (N = size(T,2)) == length(x)

	# TODO: nfft! requires that y contains only zeros. Change in NFFT package?
	fill!(y, 0.0+0.0*im)
	N = Int(sqrt(N))
	X = reshape_view(x, (N,N))

	NFFT.nfft!(T.FFT, X, y)
	had!(y, T.diag)

	return y
end

function Base.(:(*))(T::Freq2wave2D, x::Vector)
	y = Array(Complex{Float64}, size(T,1))
	mul!(T, complex(x), y)
end

function mulT!(T::Freq2wave2D, v::Vector{Complex{Float64}}, z::Vector{Complex{Float64}})
	@assert size(T,1) == length(v)
	@assert (N = size(T,2)) == length(z)

	D = conj(T.diag)
	had!(D, v)
	# TODO: As in mul!
	fill!(z, 0.0+0.0*im)
	N = Int(sqrt(N))
	Z = reshape_view(z, (N,N))
	NFFT.nfft_adjoint!(T.FFT, D, Z)

	return vec(z)
end

function Base.Ac_mul_B(T::Freq2wave2D, x::Vector)
	N = size(T,2)
	y = Array(Complex{Float64}, N)
	mulT!(T, complex(x), y)
end


function Base.collect(T::Freq2wave2D)
	# TODO: Check if the matrix fits in memory

	# Fourier matrix
	M, NN = size(T)
	F = Array(Complex{Float64}, M, NN)
	N = Int(sqrt(NN))
	for m = 1:M
		samplex = T.samples[m,1]
		sampley = T.samples[m,2]
		x = cis(-2*pi*samplex*[0:N-1;]/N)
		y = cis(-2*pi*sampley*[0:N-1;]/N)

		F[m,:] = kron(x,y)
	end

	broadcast!(*, F, T.column1, F)
end

