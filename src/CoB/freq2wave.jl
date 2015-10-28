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
function freq2wave(samples::Vector, wave::AbstractString, J::Int; B::Float64=0.0)
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

@doc """
	isuniform(Freq2wave)

Are the samples used in change of basis matrix uniform.
"""->
function isuniform(T::Freq2wave)
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

function Base.size(T::Freq2wave1D, ::Type{Val{1}})
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


function freq2wave(samples::AbstractMatrix, wave::AbstractString, J::Int; B::Float64=0.0)
	M = size(samples, 1)
	@assert size(samples,2) == 2
	# TODO: Warning if J is too big

	# Evaluate the first column in change of basis matrix
	# TODO: Parse strings such as "db4"
	func = string("Four", wave, "Scaling")
	samplesx = view(samples, :, 1)
	samplesy = view(samples, :, 2)
	column1 = eval(parse(func))( samplesx, J )
	column1y = eval(parse(func))( samplesy, J )
	had!(column1, column1y)

	# NFFTPlans: Frequencies must be in the torus [-1/2, 1/2)^2
	N = 2^J
	xi = samples'
	scale!(xi, 1/N)
	frac!(xi)
	p = NFFTPlan(xi, (N,N))

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

	M = size(T,1)
	N = wside(T)
	println(io, "From: ", M, U, "uniform frequency samples")
	print(io, "To: ", N, "-by-", N, " ", T.wave, " wavelets")
end

function Base.size(T::Freq2wave2D, ::Type{Val{1}})
	size(T.samples, 1)
end

function Base.size(T::Freq2wave2D, ::Type{Val{2}})
	prod( T.FFT.N )
end

function mul!(T::Freq2wave2D, X::DenseArray{Complex{Float64},2}, y::Vector{Complex{Float64}})
	@assert size(T,1) == length(y)
	N = wside(T)
	@assert (N,N) == size(X)

	# TODO: nfft! requires that y contains only zeros. Change in NFFT package?
	fill!(y, 0.0+0.0*im)
	NFFT.nfft!(T.FFT, X, y)
	had!(y, T.diag)

	return y
end

function mul!(T::Freq2wave2D, x::Vector{Complex{Float64}}, y::Vector{Complex{Float64}})
	@assert size(T,1) == length(y)
	@assert size(T,2) == length(x)

	N = wside(T)
	X = reshape_view(x, (N,N))
	mul!(T, X, y)

	return y
end

@doc """
	*(T::Freq2wave, x::vector)

Compute `T*x`.
"""->
function Base.(:(*))(T::Freq2wave, x::Vector)
	y = Array(Complex{Float64}, size(T,1))
	mul!(T, complex(x), y)
	return y
end


function mulT!(T::Freq2wave2D, v::Vector{Complex{Float64}}, z::Vector{Complex{Float64}})
	@assert size(T,1) == length(v)
	@assert size(T,2) == length(z)

	N = wside(T)
	Z = reshape_view(z, (N,N))
	mulT!(T, v, Z)

	return z
end

function mulT!(T::Freq2wave2D, v::Vector{Complex{Float64}}, Z::DenseArray{Complex{Float64},2})
	@assert size(T,1) == length(v)
	N = wside(T)
	@assert (N,N) == size(Z)

	# TODO: Save conj(T.diag) ?
	D = conj(T.diag)
	had!(D, v)
	# TODO: As in mul!
	fill!(Z, 0.0+0.0*im)
	NFFT.nfft_adjoint!(T.FFT, D, Z)

	return Z
end

@doc """
	Ac_mul_B(T::Freq2wave, x::vector)
	'*(T::Freq2wave1D, x::vector)

Compute `T'*x`.
"""->
function Base.Ac_mul_B(T::Freq2wave, x::Vector)
	y = Array(Complex{Float64}, size(T,2))
	mulT!(T, complex(x), y)
	return y
end


@doc """
	NDFT(xsample::Vector{Float64}, ysample::Vector{Float64}, N::Int)
	
Compute the non-uniform Fourier matrix `F` from an `N-by-N` grid centered
around the origin (as in `grid`) to samples `(xsample,ysample)`.
The grid is assumed sorted by the `y` coordinate, i.e., the order is
(0,0),
(1,0),
(2,0),
(0,1),
(1,1),
(2,1)
etc
"""->
function NDFT(xsample::Vector{Float64}, ysample::Vector{Float64}, N::Int)
	# TODO: Check if the matrix fits in memory
	@assert (M = length(xsample)) == length(ysample)
	F = Array(Complex{Float64}, M, N^2)

	transl = [0:N-1;] - N/2
	idx = 0
	for ny = 1:N
		transly = transl[ny]
		for nx = 1:N
			translx = transl[nx]
			idx += 1
			for m = 1:M
				@inbounds F[m,idx] = cis( -2*pi*(translx*xsample[m] + transly*ysample[m]) )
			end
		end
	end

	return F
end

function Base.collect(T::Freq2wave2D)
	N = wside(T)
	xsample = T.samples[:,1]/N
	ysample = T.samples[:,2]/N

	F = NDFT(xsample, ysample, N)
	broadcast!(*, F, F, T.diag)
end

