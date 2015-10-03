# ------------------------------------------------------------
# 1D change of basis matrix

@doc """
	freq2wave(samples, wave::String, J::Int; B=0)

Make change of basis for switching between frequency responses and wavelets.

- `samples` are the sampling locations.
- `wave` is the name of the wavelet; see documentation for possibilities.
- `J` is the scale of the wavelet transform to reconstruct.
- `B` is the bandwidth of the samples; only needed if `samples` are non-uniform.

`weights` in the output is a `Nullable` type and

- If the samples *are* uniform, `weights` is `Null`.
- If the samples are *not* uniform, `weights` contains the weights as a `Nullable` vector and `basis` and `diag` are scaled with `sqrt(weights)`.
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

	if isuniform(samples)
		W = Nullable{Vector{Float64}}()
	else
		if B <= 0.0
			error("Samples are not uniform; supply bandwidth")
		end

		W = weights(samples, B)
		mu = complex( sqrt(W) )
		had!(column1, mu)
		W = Nullable(W)
	end

	diag = column1 .* cis(-pi*samples)

	return Freq2wave1D(samples, W, wave, column1, J, diag, p)
end

function isuniform(T::Freq2wave; prec::Float64=eps())
	isnull(T.weights)
end

function Base.show(io::IO, T::Freq2wave1D)
	println(io, "1D change of basis matrix")

	isuniform(T) ?  U = " " : U = " non-"

	M, N = size(T)
	println(io, "From: ", M, U, "uniform frequency samples")
	print(io, "To: ", N, " ", T.wave, " wavelets")
end


# ------------------------------------------------------------
# Basic operations for Freq2wave1D

function Base.size(T::Freq2wave, ::Type{Val{1}})
	length(T.samples)
end

function Base.size(T::Freq2wave1D, ::Type{Val{2}})
	T.FFT.N[1]
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
	*(T::Freq2wave, x::vector)

Compute `T*x`.
"""->
function Base.(:(*))(T::Freq2wave1D, x::Vector)
	y = Array(Complex{Float64}, size(T,1))
	mul!(T, complex(x), y)
	return y
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
	M, N = size(T)
	F = Array(Complex{Float64}, M, N)
	for n = 1:N
		for m = 1:M
			@inbounds F[m,n] = T.column1[m]*cis( -2*pi*T.samples[m]*(n-1)/N )
		end
	end

	return F
end


# ------------------------------------------------------------
# 2D change of basis matrix

@doc """
	wscale(Freq2wave)

Return the scale of the wavelet coefficients.
"""->
function wscale(T::Freq2wave)
	T.J
end

@doc """
	wside(T::Freq2wave2D)

The number of wavelet coefficients in each dimension.
"""->
function wside(T::Freq2wave2D)
	T.FFT.N[1]
end


function freq2wave(samples::AbstractMatrix, wave::String, J::Int; B::Float64=0.0)
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

	# NFFTPlans: Frequencies must be in the torus [-1/2, 1/2)
	# TODO: This is for functions on the unit square. Should rectangles be allowed?
	N = 2^J
	xi = scale(samples, 1/N)
	frac!(xi)
	p = NFFTPlan(xi', (N,N))

	if isuniform(samples)
		W = Nullable{Vector{Float64}}()
	else
		if B <= 0.0
			error("Samples are not uniform; supply bandwidth")
		end

		W = weights(samples, B)
		mu = complex( sqrt(W) )
		had!(column1, mu)
		W = Nullable(W)
	end

	diag = column1 .* cis( -pi*(samplesx + samplesy) )

	return Freq2wave2D(samples, W, wave, column1, J, diag, p)
end


function Base.show(io::IO, T::Freq2wave2D)
	println(io, "2D change of basis matrix")

	isuniform(T) ?  U = " " : U = " non-"

	M, N = size(T)
	N = wside(T)
	println(io, "From: ", M, U, "uniform frequency samples")
	print(io, "To: ", N, "-by-", N, " ", T.wave, " wavelets")
end

function Base.size(T::Freq2wave2D, ::Type{Val{2}})
	prod( T.FFT.N )
end

function mul!(T::Freq2wave2D, x::Vector{Complex{Float64}}, y::Vector{Complex{Float64}})
	@assert size(T,1) == length(y)
	@assert size(T,2) == length(x)

	# TODO: nfft! requires that y contains only zeros. Change in NFFT package?
	fill!(y, 0.0+0.0*im)
	N = wside(T)
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
	@assert size(T,2) == length(z)

	D = conj(T.diag)
	had!(D, v)
	# TODO: As in mul!
	fill!(z, 0.0+0.0*im)
	N = wside(T)
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
	# Fourier matrix
	M, N = size(T)
	F = Array(Complex{Float64}, M, N)
	K = wside(T)

	samplesx = view(T.samples, :, 1)
	samplesy = view(T.samples, :, 2)

	for n = 1:N
		# Manual kron for each row
		r = rem(n,K)
		if r == 0
			r = K
		end
		q = div(n-r,K)
		r -= 1

		for m = 1:M
			@inbounds F[m,n] = T.column1[m]*cis( -2*pi*(samplesx[m]*q + samplesy[m]*r)/K )
		end
	end

	return F
end

