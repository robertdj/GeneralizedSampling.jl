# ------------------------------------------------------------
# 1D change of basis matrix

function Freq2Wave(samples::DenseVector, wavename::AbstractString, J::Int, B::Float64=NaN; args...)
	vm = van_moment(wavename)
	Nint = 2^J
	Nint >= 2*vm-1 || throw(AssertionError("Scale it not large enough for this wavelet"))
	
	# TODO: Good condition?
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

	# Diagonal matrix used in 'multiplication' with NFFT
	diag = Vector{Complex{Float64}}(M)
	for m in 1:M
		diag[m] = internal[m] * cis(-pi*samples[m])
	end

	# The number of internal wavelets in reconstruction
	if hasboundary(wavename)
		Nint -= 2*van_moment(wavename)
		Nint <= 0 && error("Too few wavelets: Boundary functions overlap")
	end

	# NFFTPlans: Frequencies must be in the torus [-1/2, 1/2)
	# TODO: Should window width m and oversampling factor sigma be changed for higher precision?
	xi = samples*2.0^(-J)
	frac!(xi)
	p = NFFTPlan(xi, Nint)

	# Wavelets w/o boundary
	if !hasboundary(wavename)
		return Freq2NoBoundaryWave1D(internal, W, J, wavename, diag, p)
	else
		left = FourScalingFunc( samples, wavename, 'L', J; args... )
		right = FourScalingFunc( samples, wavename, 'R', J; args... )

		return Freq2BoundaryWave1D(internal, W, J, wavename, diag, p, left, right)
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
	for n = 0:N-1
		for m = 1:M
			@inbounds F[m,n+1] = T.internal[m]*cis( -twoπ*n*T.NFFT.x[m] )
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

	p = van_moment(T)

	# Left boundary
	F[:,1:p] = T.left[1]

	# Internal function
	for n in p:N-p-1
		for m in 1:M
			@inbounds F[m,n+1] = T.internal[m]*cis( -twoπ*n*T.NFFT.x[m] )
		end
	end

	# Right boundary
	F[:,N-p+1:N] = T.right[1]

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

The size of the reconstructed wavelet coefficients.

- When `D` == 1, the output is (Int,)
- When `D` == 2, the output is (Int,Int)
"""->
function wsize(T::Freq2BoundaryWave2D)
	N = 2^wscale(T)
	return (N, N)
end
wsize(T::Freq2NoBoundaryWave1D) = T.NFFT.N
wsize(T::Freq2BoundaryWave1D) = (2^wscale(T),)
wsize(T::Freq2NoBoundaryWave2D) = T.NFFT.N

function Freq2Wave(samples::DenseMatrix, wavename::AbstractString, J::Int, B::Float64=NaN; args...)
	vm = van_moment(wavename)
	Nint = 2^J
	Nint >= 2*vm-1 || throw(AssertionError("Scale it not large enough for this wavelet"))
	# TODO: 2 columns *or* 2 rows
	M = size(samples, 1)
	size(samples,2) == 2 || throw(DimensionMismatch("Samples must have two columns"))

	# TODO: Good condition?
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
	internal = FourScalingFunc( samples, wavename, J )

	# Diagonal matrix for multiplication with internal scaling function
	xi = samples'
	diag = internal.'
	for m in 1:M, d in 1:2
		diag[d,m] *= cis( -pi*xi[d,m] )
	end

	# The number of internal wavelets in reconstruction
	if hasboundary(wavename)
		Nint -= 2*van_moment(wavename)
		Nint <= 0 && error("Too few wavelets: Boundary functions overlap")
	end

	# NFFTPlans: Frequencies must be in the torus [-1/2, 1/2)^2
	scale!(xi, 2.0^(-J))
	frac!(xi)
	p = NFFTPlan(xi, (Nint,Nint))

	# Wavelets w/o boundary
	if !hasboundary(wavename)
		return Freq2NoBoundaryWave2D(internal, W, J, wavename, diag, p)
	else
		samplesx = slice(samples, :, 1)
		samplesy = slice(samples, :, 2)

		left = cell(2,vm+1)
		left[1] = FourScalingFunc( samplesx, wavename, 'L', J; args... )
		left[2] = FourScalingFunc( samplesy, wavename, 'L', J; args... )

		right = cell(2,vm+1)
		right[1] = FourScalingFunc( samplesx, wavename, 'R', J; args... )
		right[2] = FourScalingFunc( samplesy, wavename, 'R', J; args... )

		# Views of the columns in left/right
		for k in 1:vm, d in 1:2
			left[d,k+1]  = slice( left[d], :, k )
			right[d,k+1] = slice( right[d], :, k )
		end

		return Freq2BoundaryWave2D(internal, W, J, wavename, diag, p, left, right)
	end
end

left(T::Freq2BoundaryWave2D, d::Integer, k::Integer) = T.left[d,k+1]
right(T::Freq2BoundaryWave2D, d::Integer, k::Integer) = T.right[d,k+1]

function Base.A_mul_B!(y::DenseVector{Complex{Float64}}, T::Freq2NoBoundaryWave1D, x::DenseVector{Complex{Float64}})
	(M = size(T,1)) == length(y) || throw(DimensionMismatch())
	size(T,2) == length(x) || throw(DimensionMismatch())

	NFFT.nfft!(T.NFFT, x, y)
	had!(y, T.diag)

	isuniform(T) || had!(y, get(T.weights))

	return y
end

function Base.A_mul_B!(y::DenseVector{Complex{Float64}}, T::Freq2NoBoundaryWave2D, X::DenseMatrix{Complex{Float64}})
	(M = size(T,1)) == length(y) || throw(DimensionMismatch())
	wsize(T) == size(X) || throw(DimensionMismatch())

	NFFT.nfft!(T.NFFT, X, y)
	for m in 1:M, d in 1:2
		@inbounds y[m] *= T.diag[d,m]
	end

	isuniform(T) || had!(y, get(T.weights))

	return y
end


function Base.A_mul_B!(y::DenseVector{Complex{Float64}}, T::Freq2BoundaryWave1D, x::DenseVector{Complex{Float64}})
	(M = size(T,1)) == length(y) || throw(DimensionMismatch())
	size(T,2) == length(x) || throw(DimensionMismatch())

	xleft, xint, xright = split(x, van_moment(T))

	# Internal scaling function
	NFFT.nfft!(T.NFFT, xint, y)
	had!(y, T.diag)

	# Contribution from the boundaries
	BLAS.gemv!('N', ComplexOne, T.left, xleft, ComplexOne, y)
	BLAS.gemv!('N', ComplexOne, T.right, xright, ComplexOne, y)

	isuniform(T) || had!(y, get(T.weights))

	return y
end

function Base.A_mul_B!(y::DenseVector{Complex{Float64}}, T::Freq2BoundaryWave2D, X::DenseMatrix{Complex{Float64}}, d::Integer, k::Integer)
	p = van_moment(T)
	N = size(X,1)
	xint = slice(X, p+1:N-p, k)

	# Internal scaling function
	if d == 1
		NFFT.nfft!(T.NFFTx, xint, y)
	elseif d == 2
		NFFT.nfft!(T.NFFTy, xint, y)
	else
		throw(DimensionMismatch())
	end
	for m in 1:length(y)
		@inbounds y[m] *= T.diag[d,m]
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
	
	vm = van_moment(T)
	S = split(X, vm)

	# Internal scaling functions
	NFFT.nfft!(T.NFFT, S.internal, y)
	for m in 1:M, d in 1:2
		@inbounds y[m] *= T.diag[d,m]
	end

	# Update y with border contributions
	for k in 1:vm
		# Left
		A_mul_B!( T.innery, T, X, 1, k )
		yphad!(y, left(T,2,k), T.innery)

		# Right
		A_mul_B!( T.innery, T, X, 1, N[2]-vm+k )
		yphad!(y, right(T,2,k), T.innery)

		# Upper
		x = slice(X, k, vm+1:N[2]-vm)
		NFFT.nfft!( T.NFFTy, x, T.innery )
		for m in 1:M
			@inbounds T.innery[m] *= T.diag[2,m]
		end
		yphad!(y, left(T,1,k), T.innery)

		# Lower
		x = slice(X, N[1]-vm+k, vm+1:N[2]-vm)
		NFFT.nfft!( T.NFFTy, x, T.innery )
		for m in 1:M
			@inbounds T.innery[m] *= T.diag[2,m]
		end
		yphad!(y, right(T,1,k), T.innery)
	end

	isuniform(T) || had!(y, get(T.weights))

	return y
end

function Base.A_mul_B!(y::DenseVector{Complex{Float64}}, T::Freq2Wave2D, x::DenseVector{Complex{Float64}})
	size(T,2) == length(x) || throw(DimensionMismatch())

	X = reshape_view(x, wsize(T))
	A_mul_B!(y, T, X)

	return y
end

function Base.(:(*))(T::Freq2Wave, x::DenseArray)
	if !isa(x, Array{Complex{Float64}})
		x = map(Complex{Float64}, x)
	end

	y = Array{Complex{Float64}}(size(T,1))
	A_mul_B!(y, T, x)

	return y
end


#=
function Base.Ac_mul_B!{D}(Z::DenseArray{Complex{Float64},D}, T::Freq2NoBoundaryWave{D}, v::DenseVector{Complex{Float64}})
	size(T,1) == length(v) || throw(DimensionMismatch())
	wsize(T) == size(Z) || throw(DimensionMismatch())

	Tdiag = vec(prod(T.diag, 2))
	conj!(Tdiag)
	isuniform(T) || had!(Tdiag, get(T.weights))
	had!(Tdiag, v)

	NFFT.nfft_adjoint!(T.NFFT, Tdiag, Z)

	return Z
end

function Base.Ac_mul_B!(z::DenseVector{Complex{Float64}}, T::Freq2BoundaryWave{1}, v::DenseVector{Complex{Float64}})
	size(T,1) == length(v) || throw(DimensionMismatch())
	(Nz = size(T,2)) == length(z) || throw(DimensionMismatch())

	zleft, zint, zright = split(z, van_moment(T))

	# TODO: Is this copying better than had!(*, get(T.weights)) for every part of S?
	if isuniform(T) 
		weigthedV = v
	else
		weigthedV = v .* get(T.weights)
	end

	# Internal scaling function
	Tdiag = conj(T.diag)
	had!(Tdiag, weigthedV)
	NFFT.nfft_adjoint!(T.NFFT, Tdiag, zint)

	# Update z with boundary contributions
	Ac_mul_B!( zleft, T.left, weigthedV )
	Ac_mul_B!( zright, T.right, weigthedV )

	return z
end

function Base.Ac_mul_B!(Z::DenseMatrix{Complex{Float64}}, T::Freq2BoundaryWave{2}, v::DenseVector{Complex{Float64}})
	size(T,1) == length(v) || throw(DimensionMismatch())
	wsize(T) == size(Z) || throw(DimensionMismatch())
	
	vm = van_moment(T)
	S = split(Z, vm)

	# TODO: Is this copying better than had!(*, get(T.weights)) for every part of S?
	if isuniform(T) 
		weigthedV = v
	else
		weigthedV = v .* get(T.weights)
	end

	# Internal coefficients
	innerv = copy(weigthedV)
	hadc!(innerv, T.diag[:,1])
	hadc!(innerv, T.diag[:,2])
	NFFT.nfft_adjoint!(T.NFFT, innerv, S.II)

	# Temporary arrays for holding results of the "inner" multiplication.
	# Assuming that T.NFFT.N[1] == T.NFFT.N[2] only one cornercol and
	# sidecol is needed
	cornercol = Array{Complex{Float64}}(vm)
	Nint = size(S.IL,1)
	sidecol = Array{Complex{Float64}}(Nint)

	p1 = NFFTPlan( vec(T.NFFT.x[1,:]), T.NFFT.N[1] )
	p2 = NFFTPlan( vec(T.NFFT.x[2,:]), T.NFFT.N[2] )
	for k in 1:vm
		# LL
		conj!(innerv, T.left[:,k,2])
		had!(innerv, weigthedV)
		Ac_mul_B!( cornercol, T.left[:,:,1], innerv )
		copy!( S.LL, 1+(k-1)*vm, cornercol, 1, vm )

		# RL: Reusing innerv
		Ac_mul_B!( cornercol, T.right[:,:,1], innerv )
		copy!( S.RL, 1+(k-1)*vm, cornercol, 1, vm )

		# IL: Reusing innerv
		hadc!(innerv, T.diag[:,1])
		NFFT.nfft_adjoint!(p1, innerv, sidecol)
		copy!( S.IL, 1+(k-1)*Nint, sidecol, 1, Nint )

		# LI
		conj!(innerv, T.left[:,k,1])
		had!(innerv, weigthedV)
		hadc!(innerv, T.diag[:,2])
		NFFT.nfft_adjoint!(p2, innerv, sidecol)
		# The LI and RI component should be transposed in order to copy
		# sidecol results in-place. 
		# TODO: Can ReshapedArrays help with this (when they are ready)?
		# TODO: Try home-made unsafe_copy!
		#= BLAS.blascopy!( Nint, sidecol, 1, S.LI, vm ) =#
		#= copy!( S.LI, 1+(k-1)*vm, sidecol, 1, vm ) =#
		for n in 1:Nint
			S.LI[k,n] = sidecol[n]
		end

		# RI
		conj!(innerv, T.right[:,k,1])
		had!(innerv, weigthedV)
		hadc!(innerv, T.diag[:,2])
		NFFT.nfft_adjoint!(p2, innerv, sidecol)
		for n in 1:Nint
			S.RI[k,n] = sidecol[n]
		end

		# LR
		conj!(innerv, T.right[:,k,2])
		had!(innerv, weigthedV)
		Ac_mul_B!( cornercol, T.left[:,:,1], innerv )
		copy!( S.LR, 1+(k-1)*vm, cornercol, 1, vm )

		# RR: Reusing innerv
		Ac_mul_B!( cornercol, T.right[:,:,1], innerv )
		copy!( S.RR, 1+(k-1)*vm, cornercol, 1, vm )

		# IR: Reusing innerv
		hadc!(innerv, T.diag[:,1])
		NFFT.nfft_adjoint!(p1, innerv, sidecol)
		copy!( S.IR, 1+(k-1)*Nint, sidecol, 1, Nint )
	end

	return Z
end

function Base.Ac_mul_B!(z::DenseVector{Complex{Float64}}, T::Freq2Wave{2}, v::DenseVector{Complex{Float64}})
	size(T,2) == length(z) || throw(DimensionMismatch())

	Z = reshape_view(z, wsize(T))
	Ac_mul_B!(Z, T, v)

	return z
end
=#

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

	phi = prod(T.internal, 2)

	idx = 0
	for ny in 0:Ny-1, nx in 0:Nx-1
		idx += 1
		for m in 1:M
			@inbounds F[m,idx] = phi[m]*cis( -twoπ*(nx*T.NFFT.x[1,m] + ny*T.NFFT.x[2,m]) )
		end
	end

	if !isuniform(T)
		broadcast!( *, F, F, get(T.weights) )
	end

	return F
end

function ucollect(T::Freq2NoBoundaryWave2D)
	isuniform(T) || throw(AssertionError())

	M1, M2 = T.NFFT.n
	N1, N2 = T.NFFT.N

	F1 = Array{Complex{Float64}}(M1, N1)
	for n in 1:N1, m in 1:M2:M1
		@inbounds F1[m,n] = T.internal[m]*cis( -twoπ*T.NFFT.x[m]*(n-1) )
	end

	F2 = Array{Complex{Float64}}(M2, N2)
	for n in 1:N1, m in 1:M2
		@inbounds F1[m,n] = T.internal[m]*cis( -twoπ*T.NFFT.x[m]*(n-1) )
	end

	#kron(F1, F2)
	return F1, F2
end

function Base.collect(T::Freq2BoundaryWave2D)
	M = size(T,1)
	Nx, Ny = wsize(T)
	F = Array{Complex{Float64}}(M, Nx*Ny)

	phix = Array{Complex{Float64}}(M)
	phiy = similar(phix)

	idx = 0
	for ny in 0:Ny-1
		unsafe_FourScaling!(phiy, T, ny, 2)
		for nx in 0:Nx-1
			unsafe_FourScaling!(phix, T, nx, 1)
			had!(phix, phiy)
			F[:,idx+=1] = phix
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

	if p <= n < N-p
		for m in 1:M
			@inbounds phi[m] = T.internal[m,d]*cis( -twoπ*n*T.NFFT.x[d,m] )
		end
	elseif 0 <= n < p
		so = n*M + 1
		unsafe_copy!( phi, 1, T.left[d], so, M ) 
	else
		so = (n-N+p)*M + 1
		unsafe_copy!( phi, 1, T.right[d], so, M ) 
	end
end


# ------------------------------------------------------------
# Common

function Base.size(T::Freq2Wave)
	( size(T,1), size(T,2) )
end

function Base.size(T::Freq2Wave, d::Integer)
	if d == 1
		size(T.internal, 1)
	elseif d == 2
		D = dim(T)
		2^(D*wscale(T))
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


