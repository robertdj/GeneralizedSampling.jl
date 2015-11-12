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
function freq2wave(samples::Vector, wave::AbstractString, J::Int; B::Float64=NaN)
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
		if isnan(B)
			error("Samples are not uniform; supply bandwidth")
		end

		W = weights(samples, B)
		mu = complex( sqrt(W) )
		had!(column1, mu)
		W = Nullable(W)
	end

	diag = column1 .* cis(-pi*samples)

	return Freq2wave(samples, W, wave, column1, J, diag, p)
end

@doc """
	isuniform(Freq2wave)

Are the samples used in change of basis matrix uniform.
"""->
function isuniform(T::Freq2wave)
	isnull(T.weights)
end


# ------------------------------------------------------------
# Basic operations for 1D Freq2wave

function Base.size(T::Freq2wave{1}, ::Type{Val{1}})
	length(T.samples)
end


function Base.collect(T::Freq2wave{1})
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
	dim(Freq2wave{D})

Return the dimension `D`.
"""->
function dim{D}(::Freq2wave{D})
	D
end

@doc """
	wscale(Freq2wave)

Return the scale of the wavelet coefficients.
"""->
function wscale(T::Freq2wave)
	T.J
end

@doc """
	wsize(T::Freq2wave{D})

The size of the reconstructed wavelet coefficients.

- When `D` == 1, the output is (Int,)
- When `D` == 2, the output is (Int,Int)
"""->
function wsize(T::Freq2wave)
	T.FFT.N
end


function freq2wave(samples::AbstractMatrix, wave::AbstractString, J::Int; B::Float64=NaN)
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
		if isnan(B)
			error("Samples are not uniform; supply bandwidth")
		end

		W = weights(samples, B)
		mu = complex( sqrt(W) )
		had!(column1, mu)
		W = Nullable(W)
	end

	diag = column1 .* cis( -pi*(samplesx + samplesy) )

	return Freq2wave(samples, W, wave, column1, J, diag, p)
end


function Base.show{D}(io::IO, T::Freq2wave{D})
	println(io, D, "D change of basis matrix")

	isuniform(T) ?  U = " " : U = " non-"

	M = size(T,1)
	D == 1 ? N = size(T,2) : N = wsize(T)
	println(io, "From: ", M, U, "uniform frequency samples")
	print(io, "To: ", N, " ", T.wave, " wavelets")
end

function Base.size(T::Freq2wave{2}, ::Type{Val{1}})
	size(T.samples, 1)
end

function Base.size(T::Freq2wave, ::Type{Val{2}})
	prod( T.FFT.N )
end

@doc """
	mul!(T::Freq2wave, x::Vector, y::Vector)
	
Replace `y` with `T*x`.
"""->
function mul!{D}(T::Freq2wave{D}, X::DenseArray{Complex{Float64},D}, y::Vector{Complex{Float64}})
	@assert size(T,1) == length(y)
	@assert wsize(T) == size(X)

	# TODO: nfft! requires that y contains only zeros. Change in NFFT package?
	fill!(y, 0.0+0.0*im)
	NFFT.nfft!(T.FFT, X, y)
	had!(y, T.diag)

	return y
end

function mul!(T::Freq2wave{2}, x::Vector{Complex{Float64}}, y::Vector{Complex{Float64}})
	@assert size(T,1) == length(y)
	@assert size(T,2) == length(x)

	N = wsize(T)
	X = reshape_view(x, N)
	mul!(T, X, y)

	return y
end

@doc """
	*(T::Freq2wave, x::vector)

Compute `T*x`.
"""->
function Base.(:(*))(T::Freq2wave, x::Vector)
	@assert size(T,2) == length(x)
	y = Array(Complex{Float64}, size(T,1))
	mul!(T, complex(x), y)
	return y
end


@doc """
	mulT!(T::Freq2wave, v::Vector, z::Vector)
	
Replace `z` with `T'*v`.
"""->
function mulT!{D}(T::Freq2wave{D}, v::Vector{Complex{Float64}}, Z::DenseArray{Complex{Float64},D})
	@assert size(T,1) == length(v)
	@assert wsize(T) == size(Z)

	# TODO: Save conj(T.diag) ?
	Tdiag = conj(T.diag)
	had!(Tdiag, v)
	# TODO: As in mul!
	fill!(Z, 0.0+0.0*im)
	NFFT.nfft_adjoint!(T.FFT, Tdiag, Z)

	return Z
end

function mulT!(T::Freq2wave{2}, v::Vector{Complex{Float64}}, z::Vector{Complex{Float64}})
	@assert size(T,1) == length(v)
	@assert size(T,2) == length(z)

	N = wsize(T)
	Z = reshape_view(z, N)
	mulT!(T, v, Z)

	return z
end

@doc """
	Ac_mul_B(T::Freq2wave, x::vector)
	'*(T::Freq2wave{1}, x::vector)

Compute `T'*x`.
"""->
function Base.Ac_mul_B(T::Freq2wave, x::Vector)
	@assert size(T,1) == length(x)
	y = Array(Complex{Float64}, size(T,2))
	mulT!(T, complex(x), y)
	return y
end

function Base.(:(\))(T::Freq2wave, y::Vector{Complex{Float64}})
	# Non-uniform samples: Scale observations
	if !isuniform(T)
		y = scale(y, T.weights)
	end

	# TODO: Exact solution for 2D uniform samples?

	print("Solution via conjugate gradients... ")
	# TODO: Better initial guess?
	x0 = zeros(Complex{Float64}, wsize(T))
	x = cgnr(T, y, x0)

	return x
end


@doc """
	collect(Freq2wave)
	
Return the full change of basis matrix.

In 2D, the reconstruction grid is sorted by the `y` coordinate, i.e., the order is
(1,1),
(2,1),
(3,1),
(1,2),
(2,2),
(3,2)
etc
"""->
function Base.collect(T::Freq2wave{2})
	M = size(T,1)
	N = wsize(T)
	# TODO: Check if the matrix fits in memory
	F = Array(Complex{Float64}, M, prod(N))

	xsample = T.samples[:,1]/N[1]
	ysample = T.samples[:,2]/N[2]

	translx = [0:N[1]-1;] - N[1]/2
	transly = [0:N[2]-1;] - N[2]/2
	idx = 0
	for ny = 1:N[2]
		yidx = transly[ny]
		for nx = 1:N[1]
			xidx = translx[nx]
			idx += 1
			for m = 1:M
				@inbounds F[m,idx] = T.diag[m]*cis( -2*pi*(xidx*xsample[m] + yidx*ysample[m]) )
			end
		end
	end

	return F
end

