# ------------------------------------------------------------
# 1D change of basis matrix

function Freq2Wave(samples::StridedVector, wavename::AbstractString, J::Int, B::Float64=NaN; args...)
	vm = van_moment(wavename)
	( Nint = 2^J ) >= 2*vm-1 || throw(AssertionError("Scale it not large enough for this wavelet"))
	
	M = length(samples)
	if Nint >= M
		warn("The scale is high compared to the number of samples")
	end

	# Weights for non-uniform samples
	if isuniform(samples)
		W = Nullable{ Vector{Complex{Float64}} }()
	else
		isnan(B) && error("Samples are not uniform; supply bandwidth")
		Nint <= 2*B || warn("The scale is high compared to the bandwidth")
		W = sqrt(weights(samples, B))
		W = Nullable(complex( W ))
	end

	# Fourier transform of the internal scaling function
	internal = FourScalingFunc( samples, wavename, J; args... )

	# The number of internal wavelets in the reconstruction
	if hasboundary(wavename)
		Nint -= 2*vm
		Nint <= 0 && error("Too few wavelets: Boundary functions overlap")
	end

	# NFFTPlans: Frequencies must be in the torus [-1/2, 1/2)
	# TODO: Should window width m and oversampling factor sigma be changed for higher precision?
	xi = samples*2.0^(-J)
	frac!(xi)
	p = NFFTPlan(xi, Nint)

	# Wavelets w/o boundary
	if !hasboundary(wavename)
		return Freq2NoBoundaryWave1D(internal, W, J, wavename, p)
	else
		left = FourScalingFunc( samples, wavename, 'L', J; args... )
		right = FourScalingFunc( samples, wavename, 'R', J; args... )

		return Freq2BoundaryWave1D(internal, W, J, wavename, p, left, right)
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


# ------------------------------------------------------------
# Basic operations for 1D Freq2Wave

function Base.collect(T::Freq2NoBoundaryWave1D)
	M, N = size(T)

	F = Array{Complex{Float64}}(M, N)
	offset = div(N, 2) + 1
	for n in 1:N
		for m in 1:M
			@inbounds F[m,n] = T.internal[m]*cis( -twoπ*(n-offset)*T.NFFT.x[m] )
		end
	end

	if !isuniform(T)
		broadcast!(*, F, F, get(T.weights))
	end

	return F
end

function Base.collect(T::Freq2BoundaryWave1D)
	M, N = size(T)
	F = Array{Complex{Float64}}(M, N)

	# Left boundary
	p = van_moment(T)
	F[:,1:p] = T.left

	# Internal function
	offset = div(N, 2) + 1
	for n in p+1:N-p
		for m in 1:M
			@inbounds F[m,n] = T.internal[m]*cis( -twoπ*(n-offset)*T.NFFT.x[m] )
		end
	end

	# Right boundary
	F[:,N-p+1:N] = T.right

	if !isuniform(T)
		broadcast!(*, F, F, get(T.weights))
	end

	return F
end


# ------------------------------------------------------------
# 2D change of basis matrix

@doc """
	dim(Freq2wave)

Return the dimension of the `T`.
"""->
dim(::Freq2Wave1D) = 1
dim(::Freq2Wave2D) = 2

@doc """
	wscale(Freq2Wave)

Return the scale of the wavelet coefficients.
"""->
function wscale(T::Freq2Wave)
	T.J
end

@doc """
	wsize(T::Freq2Wave)

The number of reconstructed wavelet coefficients (in each dimension).

- When `D` == 1, the output is (Int,)
- When `D` == 2, the output is (Int,Int)
"""->
function wsize(T::Freq2BoundaryWave2D)
	N = 2^wscale(T)
	return (N, N)
end
wsize(T::Freq2NoBoundaryWave2D) = T.NFFT.N
wsize(T::Freq2BoundaryWave1D) = (2^wscale(T),)
wsize(T::Freq2NoBoundaryWave1D) = T.NFFT.N

function UnifFourScalingFunc(samples::StridedMatrix{Float64}, wavename::AbstractString, J::Integer; args...)
	# Test if samples are on a grid and in the correct order
	usamplesx = unique(slice(samples,:,1))
	usamplesy = unique(slice(samples,:,2))
	Mx = length(usamplesx)
	My = length(usamplesy)
	grid_samples = grid( (Mx, My), usamplesx[2]-usamplesx[1] )

	samples == grid_samples ||
	throw(AssertionError("Samples are not on a grid or not ordered properly"))

	# Internal scaling function
	internalx = FourScalingFunc( usamplesx, wavename, J )
	internalx = kron(internalx, ones(My))
	internaly = FourScalingFunc( usamplesy, wavename, J )
	internaly = repmat(internaly, Mx)
	internal = hcat(internalx, internaly)

	# Boundary scaling functions
	if !hasboundary(wavename)
		return internal
	else
		vm = van_moment(wavename)

		leftx = FourScalingFunc( usamplesx, wavename, 'L', J; args... )
		lefty = FourScalingFunc( usamplesy, wavename, 'L', J; args... )
		left = cell(2,vm+1)
		left[1] = kron(leftx, ones(My))
		left[2] = repmat(lefty, Mx)

		rightx = FourScalingFunc( usamplesx, wavename, 'R', J; args... )
		righty = FourScalingFunc( usamplesy, wavename, 'R', J; args... )
		right = cell(2,vm+1)
		right[1] = kron(rightx, ones(My))
		right[2] = repmat(righty, Mx)

		# Views of the columns in left/right
		for k in 1:vm, d in 1:2
			left[d,k+1]  = slice( left[d], :, k )
			right[d,k+1] = slice( right[d], :, k )
		end

		return internal, left, right
	end
end

function NotUnifFourScalingFunc(samples::StridedMatrix{Float64}, wavename::AbstractString, J::Integer, B::Float64; args...)
	# Fourier transform of the internal scaling function
	internal = FourScalingFunc( samples, wavename, J )

	# Boundary scaling functions
	if !hasboundary(wavename)
		return internal
	else
		samplesx = slice(samples, :, 1)
		samplesy = slice(samples, :, 2)

		vm = van_moment(wavename)
		left = cell(2,vm+1)
		left[1] = FourScalingFunc( samplesx, wavename, 'L', J; args... )
		left[2] = FourScalingFunc( samplesy, wavename, 'L', J; args... )

		right = cell(2,vm+1)
		right[1] = FourScalingFunc( samplesx, wavename, 'R', J; args... )
		right[2] = FourScalingFunc( samplesy, wavename, 'R', J; args... )

		# Views of the columns in left/right
		vm = van_moment(wavename)
		for k in 1:vm, d in 1:2
			left[d,k+1]  = slice( left[d], :, k )
			right[d,k+1] = slice( right[d], :, k )
		end

		return internal, left, right
	end
end

function Freq2Wave(samples::StridedMatrix{Float64}, wavename::AbstractString, J::Int, B::Float64=NaN; args...)
	vm = van_moment(wavename)
	( Nint = 2^J ) >= 2*vm-1 || throw(AssertionError("Scale it not large enough for this wavelet"))
	M = size(samples, 1)
	size(samples,2) == 2 || throw(DimensionMismatch("Samples must have two columns"))

	if Nint >= M
		warn("The scale is high compared to the number of samples")
	end

	if isuniform(samples)
		W = Nullable{ Vector{Complex{Float64}} }()

		scaling_funcs = UnifFourScalingFunc(samples, wavename, J; args...)
	else
		# Weights for non-uniform samples
		isnan(B) && error("Samples are not uniform; supply bandwidth")
		Nint <= 2*B || warn("The scale is high compared to the bandwidth")
		W = sqrt(weights(samples, B))
		W = Nullable(complex( W ))

		scaling_funcs = NotUnifFourScalingFunc(samples, wavename, J, B; args...)
	end

	if hasboundary(wavename)
		Nint -= 2*vm
		Nint <= 0 && error("Too few wavelets: Boundary functions overlap")

		internal = scaling_funcs[1].'
		left = scaling_funcs[2]
		right = scaling_funcs[3]
	else
		internal = scaling_funcs.'
	end

	# NFFTPlans: Frequencies must be in the torus [-1/2, 1/2)^2
	xi = samples'
	scale!(xi, 2.0^(-J))
	frac!(xi)
	p = NFFTPlan(xi, (Nint,Nint))

	# Wavelets w/o boundary
	if hasboundary(wavename)
		return Freq2BoundaryWave2D(internal, W, J, wavename, p, left, right)
	else
		return Freq2NoBoundaryWave2D(internal, W, J, wavename, p)
	end
end

@doc """
	left(T, d, k) -> Vector

Return (a view of) the `k`'th column of `T`'s left boundary Fourier
transform in dimension `d`.
For the x coordinate `d = 1` and for the y coordinate `d = 2`.
"""->
left(T::Freq2BoundaryWave2D, d::Integer, k::Integer) = T.left[d,k+1]
left(T::Freq2BoundaryWave2D, d::Integer) = T.left[d,1]

@doc """
	right(T, d, k) -> Vector

Return (a view of) the `k`'th column of `T`'s right boundary Fourier
transform in dimension `d`.
For the x coordinate `d = 1` and for the y coordinate `d = 2`.

Note that `right(T, d, 1)` is the right boundary functions *furthest*
from the boundary, cf. the ordering from `FourScalingFunc`.
"""->
right(T::Freq2BoundaryWave2D, d::Integer, k::Integer) = T.right[d,k+1]
right(T::Freq2BoundaryWave2D, d::Integer) = T.right[d,1]


function Base.A_mul_B!(y::DenseVector{Complex{Float64}}, T::Freq2NoBoundaryWave1D, x::DenseVector{Complex{Float64}})
	size(T,1) == length(y) || throw(DimensionMismatch())
	size(T,2) == length(x) || throw(DimensionMismatch())

	nfft!(T.NFFT, x, y)
	had!(y, T.internal)

	isuniform(T) || had!(y, get(T.weights))

	return y
end

function Base.A_mul_B!(y::DenseVector{Complex{Float64}}, T::Freq2NoBoundaryWave2D, X::DenseMatrix{Complex{Float64}})
	(M = size(T,1)) == length(y) || throw(DimensionMismatch())
	wsize(T) == size(X) || throw(DimensionMismatch())

	nfft!(T.NFFT, X, y)
	for m in 1:M, d in 1:2
		@inbounds y[m] *= T.internal[d,m]
	end

	isuniform(T) || had!(y, get(T.weights))

	return y
end


function Base.A_mul_B!(y::DenseVector{Complex{Float64}}, T::Freq2BoundaryWave1D, x::DenseVector{Complex{Float64}})
	(M = size(T,1)) == length(y) || throw(DimensionMismatch())
	size(T,2) == length(x) || throw(DimensionMismatch())

	xleft, xint, xright = split(x, van_moment(T))

	# Internal scaling function
	nfft!(T.NFFT, xint, y)
	had!(y, T.internal)

	# Contribution from the boundaries
	BLAS.gemv!('N', ComplexOne, T.left, xleft, ComplexOne, y)
	BLAS.gemv!('N', ComplexOne, T.right, xright, ComplexOne, y)

	isuniform(T) || had!(y, get(T.weights))

	return y
end

function Base.A_mul_B!(y::StridedVector{Complex{Float64}}, T::Freq2BoundaryWave2D, X::DenseMatrix{Complex{Float64}}, d::Integer, k::Integer)
	p = van_moment(T)
	N = size(X,1)
	xint = slice(X, p+1:N-p, k)

	# Internal scaling function
	if d == 1
		nfft!(T.NFFTx, xint, y)
	elseif d == 2
		nfft!(T.NFFTy, xint, y)
	else
		throw(DimensionMismatch())
	end
	for m in 1:length(y)
		@inbounds y[m] *= T.internal[d,m]
	end

	# Contribution from the boundaries
	xleft = slice(X, 1:p, k)
	BLAS.gemv!('N', ComplexOne, T.left[d], xleft, ComplexOne, y)
	xright = slice(X, N-p+1:N, k)
	BLAS.gemv!('N', ComplexOne, T.right[d], xright, ComplexOne, y)

	isuniform(T) || had!(y, get(T.weights))

	return y
end

function Base.A_mul_B!(y::DenseVector{Complex{Float64}}, T::Freq2BoundaryWave2D, X::DenseMatrix{Complex{Float64}})
	(M = size(T,1)) == length(y) || throw(DimensionMismatch())
	(N = wsize(T)) == size(X) || throw(DimensionMismatch())

	# Internal scaling functions
	vm = van_moment(T)
	S = split(X, vm)
	nfft!(T.NFFT, S.internal, y)
	for m in 1:M, d in 1:2
		@inbounds y[m] *= T.internal[d,m]
	end

	# Update y with border contributions
	for k in 1:vm
		# Left
		A_mul_B!( T.tmpMulVec, T, X, 1, k )
		yphad!(y, left(T,2,k), T.tmpMulVec)

		# Right
		A_mul_B!( T.tmpMulVec, T, X, 1, N[2]-vm+k )
		yphad!(y, right(T,2,k), T.tmpMulVec)

		# Upper
		x = slice(X, k, vm+1:N[2]-vm)
		nfft!( T.NFFTy, x, T.tmpMulVec )
		for m in 1:M
			@inbounds T.tmpMulVec[m] *= T.internal[2,m]
		end
		yphad!(y, left(T,1,k), T.tmpMulVec)

		# Lower
		x = slice(X, N[1]-vm+k, vm+1:N[2]-vm)
		nfft!( T.NFFTy, x, T.tmpMulVec )
		for m in 1:M
			@inbounds T.tmpMulVec[m] *= T.internal[2,m]
		end
		yphad!(y, right(T,1,k), T.tmpMulVec)
	end

	isuniform(T) || had!(y, get(T.weights))

	return y
end

function Base.A_mul_B!(y::DenseVector{Complex{Float64}}, T::Freq2Wave2D, x::DenseVector{Complex{Float64}})
	size(T,2) == length(x) || throw(DimensionMismatch())

	X = reshape(x, wsize(T))
	A_mul_B!(y, T, X)

	return y
end

function Base.(:*)(T::Freq2Wave, x::DenseArray)
	if !isa(x, Array{Complex{Float64}})
		x = map(Complex{Float64}, x)
	end

	y = Array{Complex{Float64}}( size(T,1) )
	A_mul_B!(y, T, x)

	return y
end


function Base.Ac_mul_B!(z::DenseVector{Complex{Float64}}, T::Freq2NoBoundaryWave1D, v::DenseVector{Complex{Float64}})
	(M = size(T,1)) == length(v) || throw(DimensionMismatch())
	size(T,2) == length(z) || throw(DimensionMismatch())

	hadc!(T.tmpMulVec, v, T.internal)
	isuniform(T) || had!(T.tmpMulVec, get(T.weights))

	nfft_adjoint!(T.NFFT, T.tmpMulVec, z)

	return z
end

function Base.Ac_mul_B!(Z::DenseMatrix{Complex{Float64}}, T::Freq2NoBoundaryWave2D, v::DenseVector{Complex{Float64}})
	(M = size(T,1)) == length(v) || throw(DimensionMismatch())
	wsize(T) == size(Z) || throw(DimensionMismatch())

	for m in 1:M
		@inbounds T.tmpMulVec[m] = v[m]
		for d in 1:2
			@inbounds T.tmpMulVec[m] *= conj(T.internal[d,m])
		end
	end
	isuniform(T) || had!(T.tmpMulVec, get(T.weights))

	nfft_adjoint!(T.NFFT, T.tmpMulVec, Z)

	return Z
end

function Base.Ac_mul_B!(z::DenseVector{Complex{Float64}}, T::Freq2BoundaryWave1D, v::DenseVector{Complex{Float64}})
	(M = size(T,1)) == length(v) || throw(DimensionMismatch())
	size(T,2) == length(z) || throw(DimensionMismatch())

	zleft, zint, zright = split(z, van_moment(T))

	# Boundary contributions are first as they don't use T.internal
	copy!(T.tmpMulVec, v)
	isuniform(T) || had!(T.tmpMulVec, get(T.weights))
	Ac_mul_B!( zleft, T.left, T.tmpMulVec )
	Ac_mul_B!( zright, T.right, T.tmpMulVec )

	# Internal scaling function
	hadc!(T.tmpMulVec, T.internal)
	nfft_adjoint!(T.NFFT, T.tmpMulVec, zint)

	return z
end

# Ac_mul_B! for the `k`'th column of `Z` using the `d`'th dimension of T
function Base.Ac_mul_B!(Z::StridedMatrix{Complex{Float64}}, T::Freq2BoundaryWave2D, v::DenseVector{Complex{Float64}}, d::Integer, k::Integer)
	copy!(T.tmpMulVec, v)
	isuniform(T) || had!(T.tmpMulVec, get(T.weights))

	vm = van_moment(T)
	N = size(Z,1)

	zleft = slice(Z, 1:vm, k)
	Ac_mul_B!( zleft, left(T,d), T.tmpMulVec )
	zright = slice(Z, N-vm+1:N, k)
	Ac_mul_B!( zright, right(T,d), T.tmpMulVec )

	# Internal scaling function
	for m in 1:length(v)
		@inbounds T.tmpMulVec[m] *= conj(T.internal[d,m])
	end

	zint = slice(Z, vm+1:N-vm, k)
	if d == 1
		nfft_adjoint!(T.NFFTx, T.tmpMulVec, zint)
	elseif d == 2
		nfft_adjoint!(T.NFFTy, T.tmpMulVec, zint)
	else
		throw(DimensionMismatch())
	end

	return Z
end

function Base.Ac_mul_B!(Z::DenseMatrix{Complex{Float64}}, T::Freq2BoundaryWave2D, v::DenseVector{Complex{Float64}})
	(M = size(T,1)) == length(v) || throw(DimensionMismatch())
	(N = wsize(T)) == size(Z) || throw(DimensionMismatch())
	
	# As in Ac_mul_B! for Freq2BoundaryWave1D
	copy!(T.weigthedVec, v)
	isuniform(T) || had!(T.weigthedVec, get(T.weights))

	vm = van_moment(T)
	S = split(Z, vm)

	# Update borders
	for k in 1:vm
		# Left
		hadc!(T.tmpMulcVec, T.weigthedVec, left(T,2,k))
		Ac_mul_B!( Z, T, T.tmpMulcVec, 1, k )

		# Right
		hadc!(T.tmpMulcVec, T.weigthedVec, right(T,2,k))
		Ac_mul_B!( Z, T, T.tmpMulcVec, 1, N[2]-vm+k )

		# Upper
		hadc!( T.tmpMulcVec, T.weigthedVec, left(T,1,k) )
		for m in 1:M
			@inbounds T.tmpMulcVec[m] *= conj(T.internal[2,m])
		end
		z = slice(Z, k, vm+1:N[2]-vm)
		nfft_adjoint!( T.NFFTy, T.tmpMulcVec, z )

		# Lower
		hadc!( T.tmpMulcVec, T.weigthedVec, right(T,1,k) )
		for m in 1:M
			@inbounds T.tmpMulcVec[m] *= conj(T.internal[2,m])
		end
		z = slice(Z, N[1]-vm+k, vm+1:N[2]-vm)
		nfft_adjoint!( T.NFFTy, T.tmpMulcVec, z )
	end

	# Internal coefficients
	for m in 1:M
		@inbounds T.tmpMulVec[m] = v[m]
		for d in 1:2
			@inbounds T.tmpMulVec[m] *= conj(T.internal[d,m])
		end
	end
	isuniform(T) || had!(T.tmpMulVec, get(T.weights))
	nfft_adjoint!(T.NFFT, T.tmpMulVec, S.internal)

	return Z
end

function Base.Ac_mul_B!(z::DenseVector{Complex{Float64}}, T::Freq2Wave2D, v::DenseVector{Complex{Float64}})
	size(T,2) == length(z) || throw(DimensionMismatch())

	Z = reshape(z, wsize(T))
	Ac_mul_B!(Z, T, v)

	return z
end

function Base.Ac_mul_B(T::Freq2Wave, v::AbstractVector)
	if !isa(v, Array{Complex{Float64}})
		v = map(Complex{Float64}, v)
	end

	z = Array{Complex{Float64}}(size(T,2))
	Ac_mul_B!(z, T, v)

	return z
end

function Base.(:(\))(T::Freq2Wave, Y::AbstractMatrix)
	length(Y) == (M = size(T,1)) || throw(DimensionMismatch())

	y = flatten_view(Y)
	x = T \ y
end

function Base.(:(\))(T::Freq2Wave, y::AbstractVector)
	if !isa(y, Array{Complex{Float64}})
		y = map(Complex{Float64}, y)
	end

	# Non-uniform samples: Scale observations
	if !isuniform(T)
		y .*= get(T.weights)
	end

	x0 = zeros(Complex{Float64}, wsize(T))
	x = cgnr(T, y, x0)
	#x, h = lsqr(T, y)
end


@doc """
	collect(Freq2Wave) -> Matrix
	
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
function Base.collect(T::Freq2NoBoundaryWave2D)
	M = size(T,1)
	Nx, Ny = wsize(T)
	F = Array{Complex{Float64}}(M, Nx*Ny)

	phi = prod(T.internal, 1)

	offsetx = div(Nx, 2) + 1
	offsety = div(Ny, 2) + 1
	row_idx = 0
	for ny in 1:Ny, nx in 1:Nx
		row_idx += 1
		for m in 1:M
			@inbounds F[m,row_idx] = phi[m]*cis( -twoπ*((nx-offsetx)*T.NFFT.x[1,m] + (ny-offsety)*T.NFFT.x[2,m]) )
		end
	end

	if !isuniform(T)
		broadcast!( *, F, F, get(T.weights) )
	end

	return F
end

function Base.collect(T::Freq2BoundaryWave2D)
	M = size(T,1)
	Nx, Ny = wsize(T)
	F = Array{Complex{Float64}}(M, Nx*Ny)

	phix = Array{Complex{Float64}}(M)
	phiy = similar(phix)

	row_idx = 0
	for ny in 1:Ny
		unsafe_FourScaling!(phiy, T, ny, 2)
		for nx in 1:Nx
			unsafe_FourScaling!(phix, T, nx, 1)
			had!(phix, phiy)
			F[:,row_idx+=1] = phix
		end
	end

	if !isuniform(T)
		broadcast!( *, F, F, get(T.weights) )
	end

	return F
end

@doc """
	unsafe_FourScaling!(phi, T::Freq2BoundaryWave{2}, n::Int, d::Int)

Replace `phi` with the `n`'th "column" from dimension `d` of `T`.
"""->
function unsafe_FourScaling!(phi::Vector{Complex{Float64}}, T::Freq2BoundaryWave2D, n::Integer, d::Integer)
	M = length(phi)
	N = wsize(T)[d]
	p = van_moment(T)

	if p < n <= N-p
		offset = div(N, 2) + 1
		for m in 1:M
			@inbounds phi[m] = T.internal[d,m]*cis( -twoπ*(n-offset)*T.NFFT.x[d,m] )
		end
	elseif 1 <= n <= p
		copy!( phi, T.left[d,n+1] )
	else
		copy!( phi, T.right[d,n-N+p+1] )
	end
end


# ------------------------------------------------------------
# Common

function Base.size(T::Freq2Wave)
	( size(T,1), size(T,2) )
end

function Base.size(T::Freq2Wave, d::Integer)
	if d == 1
		size(T.internal, dim(T))
	elseif d == 2
		prod( wsize(T) )
	else
		throw(AssertionError())
	end
end

van_moment(T::Freq2Wave1D) = hasboundary(T) ? size(T.left,2) : van_moment(T.wavename)
van_moment(T::Freq2Wave2D) = hasboundary(T) ? size(T.left[1],2)::Int64 : van_moment(T.wavename)

function Base.eltype(::Freq2Wave)
	return Complex{Float64}
end

function Base.show(io::IO, T::Freq2Wave)
	D = dim(T)
	println(io, D, "D change of basis matrix")

	isuniform(T) ?  U = " " : U = " non-"
	M = size(T,1)
	println(io, "From: ", M, U, "uniform frequency samples")

	D == 1 ? N = size(T,2) : N = wsize(T)
	print(io, "To: ", N, " ", T.wavename, " wavelets")
end


