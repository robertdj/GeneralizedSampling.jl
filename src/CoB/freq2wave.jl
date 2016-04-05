# ------------------------------------------------------------
# 1D change of basis matrix

@doc """
	freq2wave(samples, wavename::String, J::Int; B=0)

Make change of basis for switching between frequency responses and wavelets.

- `samples` are the sampling locations as a vector for 1D and 2-by-M / M-by-2 matrix for 2D.
- `wave` is the name of the wavelet; see documentation for possibilities.
- `J` is the scale of the wavelet transform to reconstruct.
- `B` is the bandwidth of the samples; only needed if `samples` are non-uniform.
"""->
function freq2wave(samples::Vector, wavename::AbstractString, J::Int; B::Float64=NaN)
	# TODO: Warning if J is too big

	# Fourier transform of the internal scaling function
	FT = FourScalingFunc( samples, wavename, J )

	# NFFTPlans: Frequencies must be in the torus [-1/2, 1/2)
	# TODO: Should window width m and oversampling factor sigma be changed for higher precision?
	N = 2^J
	if hasboundary(wavename)
		F = waveparse( wavename, true )
		N -= 2*van_moment(F)
		if N <= 0
			error("Too few wavelets: Boundary functions overlap")
		end
	end
	xi = scale(samples, 1/N)
	frac!(xi)
	p = NFFTPlan(xi, N)

	# Weights for non-uniform samples
	if isuniform(samples)
		W = Nullable{ Vector{Complex{Float64}} }()
	else
		if isnan(B)
			error("Samples are not uniform; supply bandwidth")
		end

		W = Nullable(complex( sqrt(weights(samples, B)) ))
		#= W = sqrt(weights(samples, B)) =#
		#= had!( FT, complex(W) ) =#
		#= W = Nullable(W) =#
	end

	# Diagonal matrix for multiplication with internal scaling function
	diag = FT .* cis(-pi*samples)

	# Wavelets w/o boundary
	if !hasboundary(wavename)
		return Freq2NoBoundaryWave(samples, FT, W, J, wavename, diag, p)
	else
		left = FourDaubScaling(samples, F)
		right = FourDaubScaling(samples, F)
		return Freq2BoundaryWave(samples, FT, W, J, wavename, diag, p, left, right)
	end
end

@doc """
	isuniform(Freq2Wave)

Is the change of basis matrix based on uniform samples.
"""->
function isuniform(T::Freq2Wave)
	isnull(T.weights)
end

@doc """
	hasboundary(T::Freq2Wave) -> Bool

Does the wavelet in the change of basis matrix `T` have boundary correction.
"""->
function hasboundary(T::Freq2Wave)
	isdefined(T, :left)
end

@doc """
	hasboundary(wavename::String) -> Bool

Does the `wavename` scaling function have boundary correction.
"""->
function hasboundary(wavename::AbstractString)
	lowername = lowercase(wavename)

	if isdaubechies(lowername)
		if lowername == "haar" || lowername == "db1"
			return false
		else
			return true
		end
	else
		error("Only Daubechies scaling functions are valid")
	end
end


# ------------------------------------------------------------
# Basic operations for 1D Freq2Wave

function Base.collect(T::Freq2NoBoundaryWave{1})
	M, N = size(T)
	F = Array(Complex{Float64}, M, N)
	for n = 1:N
		for m = 1:M
			@inbounds F[m,n] = T.FT[m]*cis( -2*pi*T.samples[m]*(n-1)/N )
			#= @inbounds F[m,n] = T.FT[m]*cis( -2*pi*T.NFFT.x[m]*(n-1) ) =#
		end
	end

	return F
end

function Base.collect(T::Freq2BoundaryWave{1})
	M, N = size(T)
	vm = van_moment(T)
	F = Array(Complex{Float64}, M, N)

	F[:,1:vm] = T.left

	# Internal function
	for n = (vm+1):(N-vm)
		for m = 1:M
			@inbounds F[m,n] = T.FT[m]*cis( -2*pi*T.samples[m]*(n-1)/N )
			#= @inbounds F[m,n] = T.FT[m]*cis( -2*pi*T.NFFT.x[m]*(n-1) ) =#
		end
	end

	F[:,(N-vm+1):N] = T.right

	return F
end


# ------------------------------------------------------------
# 2D change of basis matrix

@doc """
	dim(Freq2wave{D})

Return the dimension `D`.
"""->
function dim{D}(::Freq2Wave{D})
	D
end

@doc """
	wscale(Freq2Wave)

Return the scale of the wavelet coefficients.
"""->
function wscale(T::Freq2Wave)
	T.J
end

@doc """
	wsize(T::Freq2Wave{D})

The size of the reconstructed wavelet coefficients.

- When `D` == 1, the output is (Int,)
- When `D` == 2, the output is (Int,Int)
"""->
function wsize(T::Freq2Wave)
	T.NFFT.N
end


function freq2wave(samples::AbstractMatrix, wavename::AbstractString, J::Int; B::Float64=NaN)
	M = size(samples, 1)
	@assert size(samples,2) == 2
	# TODO: Warning if J is too big

	# Fourier transform of the internal scaling function
	FT = FourScalingFunc( samples, wavename, J )

	# NFFTPlans: Frequencies must be in the torus [-1/2, 1/2)^2
	N = 2^J
	if hasboundary(wavename)
		F = waveparse( wavename, true )
		# TODO: How much to subtract?
		N -= 2*van_moment(F)
		if N < 0
			error("Too few wavelets: Boundary functions overlap")
		end
	end
	xi = samples'
	scale!(xi, 1/N)
	frac!(xi)
	p = NFFTPlan(xi, (N,N))

	# Weights for non-uniform samples
	if isuniform(samples)
		W = Nullable{ Vector{Complex{Float64}} }()
	else
		if isnan(B)
			error("Samples are not uniform; supply bandwidth")
		end

		W = sqrt(weights(samples, B))
		W = Nullable(complex( W ))
	end

	# Diagonal matrix for multiplication with internal scaling function
	samplesx = slice(samples, :, 1)
	samplesy = slice(samples, :, 2)
	diag = vec(prod(FT,2)) .* cis( -pi*(samplesx + samplesy) )

	# Wavelets w/o boundary
	if !hasboundary(wavename)
		return Freq2NoBoundaryWave(samples, FT, W, J, wavename, diag, p)
	else
		F = waveparse( wavename, true )
		left = [ FourDaubScaling(samplesx, F) FourDaubScaling(samplesy, F) ]
		right = [ FourDaubScaling(samplesx, F) FourDaubScaling(samplesy, F) ]
		return Freq2BoundaryWave(samples, FT, W, J, wavename, diag, p, left, right)
	end
end


function Base.A_mul_B!{D}(y::Vector{Complex{Float64}}, T::Freq2NoBoundaryWave{D}, X::AbstractArray{Complex{Float64},D})
	@assert size(T,1) == length(y)
	@assert wsize(T) == size(X)

	NFFT.nfft!(T.NFFT, X, y)
	had!(y, T.diag)

	if !isuniform(T)
		had!(y, get(T.weights))
	end

	return y
end


function Base.A_mul_B!(y::Vector{Complex{Float64}}, T::Freq2BoundaryWave{1}, x::AbstractVector{Complex{Float64}})
	@assert size(T,1) == length(y)
	@assert (Nx = size(T,2)) == length(x)

	vm = van_moment(T)

	# Internal scaling function
	NFFT.nfft!(T.NFFT, x[vm+1:Nx-vm], y)
	had!(y, T.diag)

	# Left boundary functions
	Cone = one(Complex{Float64})
	BLAS.gemv!('N', Cone, T.left, x[1:vm], Cone, y)

	# Right boundary functions
	BLAS.gemv!('N', Cone, T.right, x[Nx-vm+1:Nx], Cone, y)

	if !isuniform(T)
		had!(y, get(T.weights))
	end

	return y
end

function mul!(T::Freq2BoundaryWave{2}, X::AbstractMatrix{Complex{Float64}}, y::Vector{Complex{Float64}})
	@assert size(T,1) == length(y)
	@assert wsize(T) == size(X)
	
	error("not implemented")
end

function mul!(T::Freq2Wave{2}, x::Vector{Complex{Float64}}, y::Vector{Complex{Float64}})
	@assert size(T,1) == length(y)
	@assert size(T,2) == length(x)

	X = reshape_view(x, wsize(T))
	mul!(T, X, y)

	return y
end

@doc """
	*(T::Freq2Wave, x::vector)

Compute `T*x`.
"""->
function Base.(:(*))(T::Freq2Wave, x::AbstractVector)
	@assert size(T,2) == length(x)
	y = Array(Complex{Float64}, size(T,1))
	mul!(T, complex(x), y)
	return y
end

#=

@doc """
	mulT!(T::Freq2Wave, v::Vector, z::Vector)
	
Replace `z` with `T'*v`.
"""->
function mulT!{D}(T::Freq2Wave{D}, v::Vector{Complex{Float64}}, Z::DenseArray{Complex{Float64},D})
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

function mulT!(T::Freq2Wave{2}, v::Vector{Complex{Float64}}, z::Vector{Complex{Float64}})
	@assert size(T,1) == length(v)
	@assert size(T,2) == length(z)

	N = wsize(T)
	Z = reshape_view(z, N)
	mulT!(T, v, Z)

	return z
end

@doc """
	Ac_mul_B(T::Freq2Wave, x::vector)
	'*(T::Freq2Wave{1}, x::vector)

Compute `T'*x`.
"""->
function Base.Ac_mul_B(T::Freq2Wave, x::Vector)
	@assert size(T,1) == length(x)
	y = Array(Complex{Float64}, size(T,2))
	mulT!(T, complex(x), y)
	return y
end

function Base.(:(\))(T::Freq2Wave, y::Vector{Complex{Float64}})
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
	collect(Freq2Wave)
	
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
function Base.collect(T::Freq2NoBoundaryWave{2})
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
=#


# ------------------------------------------------------------
# Common

function Base.size(T::Freq2Wave, ::Type{Val{1}})
	size(T.samples, 1)
end

function Base.size{D}(T::Freq2Wave{D}, ::Type{Val{2}})
	2^(D*wscale(T))
end

function van_moment(T::Freq2Wave)
	van_moment(T.wavename)
end

function Base.show{D}(io::IO, T::Freq2Wave{D})
	println(io, D, "D change of basis matrix")

	isuniform(T) ?  U = " " : U = " non-"
	M = size(T,1)
	println(io, "From: ", M, U, "uniform frequency samples")

	D == 1 ? N = size(T,2) : N = wsize(T)
	print(io, "To: ", N, " ", T.wavename, " wavelets")
end

