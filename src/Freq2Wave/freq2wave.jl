# ------------------------------------------------------------
# 1D change of basis matrix

function Freq2Wave(samples::StridedVector, wavename::AbstractString, J::Int, B::Float64=NaN; uniform::Bool=true, args...)
	vm = van_moment(wavename)
	( Nint = 2^J ) >= 2*vm-1 || throw(AssertionError("Scale it not large enough for this wavelet"))
	
	M = length(samples)
	if Nint >= M
		warn("The scale is high compared to the number of samples")
	end

	# Weights for non-uniform samples
	uni_samples = isuniform(samples)
	if uni_samples
		sample_dist = samples[2] - samples[1]
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

    if uni_samples && isinteger(1/sample_dist) && uniform
	# For suitable uniform samples, the ordinary FFT can be used instead of NFFT
		p = FFTPlan(samples, J, Nint)
	else
		# NFFTPlans: Frequencies must be in the torus [-1/2, 1/2)
		xi = samples*2.0^(-J)
		frac!(xi)
		p = NFFTPlan(xi, Nint)
	end

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
		scale!(get(T.weights), F)
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
	for n in (p+1):(N-p)
		for m in 1:M
			@inbounds F[m,n] = T.internal[m]*cis( -twoπ*(n-offset)*T.NFFT.x[m] )
		end
	end

	# Right boundary
	F[:,N-p+1:N] = T.right

	if !isuniform(T)
		scale!(get(T.weights), F)
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

wsize(T::Freq2BoundaryWave1D) = (2^wscale(T),)
wsize(T::Freq2NoBoundaryWave2D) = T.NFFT.N
wsize(T::Freq2NoBoundaryWave1D) = T.NFFT.N

function UnifFourScalingFunc(samples::StridedMatrix{Float64}, wavename::AbstractString, J::Integer; args...)
	# Test if samples are on a grid and in the correct order
	usamplesx = unique(view(samples,:,1))
	usamplesy = unique(view(samples,:,2))
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
	internal = Dict( :x => internalx, :y => internaly )

	# Boundary scaling functions
	if !hasboundary(wavename)
		return internal
	else
		vm = van_moment(wavename)

		leftx = FourScalingFunc( usamplesx, wavename, 'L', J; args... )
		lefty = FourScalingFunc( usamplesy, wavename, 'L', J; args... )
		left = Dict( :x => repmat(leftx, My), :y => kron(lefty, ones(Mx)) )

		rightx = FourScalingFunc( usamplesx, wavename, 'R', J; args... )
		righty = FourScalingFunc( usamplesy, wavename, 'R', J; args... )
		right = Dict( :x => repmat(rightx, My), :y => kron(righty, ones(Mx)) )

		return internal, left, right
	end
end

function NotUnifFourScalingFunc(samples::StridedMatrix{Float64}, wavename::AbstractString, J::Integer, B::Float64; args...)
	# Fourier transform of the internal scaling function
	int = FourScalingFunc( samples, wavename, J )
	internal = Dict( :x => int[:,1], :y => int[:,2] )

	# Boundary scaling functions
	if !hasboundary(wavename)
		return internal
	else
		samplesx = view(samples, :, 1)
		samplesy = view(samples, :, 2)

		left = Dict{ Symbol, Matrix{Complex{Float64}} }()
		left[:x] = FourScalingFunc( samplesx, wavename, 'L', J; args... )
		left[:y] = FourScalingFunc( samplesy, wavename, 'L', J; args... )

		right = Dict{ Symbol, Matrix{Complex{Float64}} }()
		right[:x] = FourScalingFunc( samplesx, wavename, 'R', J; args... )
		right[:y] = FourScalingFunc( samplesy, wavename, 'R', J; args... )

		return internal, left, right
	end
end

function Freq2Wave(samples::StridedMatrix{Float64}, wavename::AbstractString, J::Int, B::Float64=NaN; uniform::Bool=true, args...)
	vm = van_moment(wavename)
	( Nint = 2^J ) >= 2*vm-1 || throw(AssertionError("Scale it not large enough for this wavelet"))
	M = size(samples, 1)
	size(samples,2) == 2 || throw(DimensionMismatch("Samples must have two columns"))

	if Nint >= M
		warn("The scale is high compared to the number of samples")
	end

	uni_samples = isuniform(samples)
	if uni_samples
		W = Nullable{ Vector{Complex{Float64}} }()

		sample_dist = samples[2,2] - samples[1,2]
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

		internal = scaling_funcs[1]
		left = scaling_funcs[2]
		right = scaling_funcs[3]
	else
		internal = scaling_funcs
	end

    xi = samples'
	if uni_samples && isinteger(1/sample_dist) && uniform
        # For suitable uniform samples, the ordinary FFT can be used instead of NFFT
        p = FFTPlan(xi, J, (Nint,Nint))
	else
        # NFFTPlans: Frequencies must be in the torus [-1/2, 1/2)^2
        scale!(xi, 2.0^(-J))
        frac!(xi)
        p = NFFTPlan(xi, (Nint,Nint))
    end

	# Wavelets w/o boundary
	if hasboundary(wavename)
		return Freq2BoundaryWave2D(internal, W, J, wavename, p, left, right)
	else
		return Freq2NoBoundaryWave2D(internal, W, J, wavename, p)
	end
end


function Base.A_mul_B!(y::DenseVector{Complex{Float64}}, T::Freq2NoBoundaryWave1D, x::DenseVector{Complex{Float64}})
	size(T,1) == length(y) || throw(DimensionMismatch())
	size(T,2) == length(x) || throw(DimensionMismatch())

	myfft!(y, T.NFFT, x)
	had!(y, T.internal)

	isuniform(T) || had!(y, get(T.weights))

	return y
end

function Base.A_mul_B!(y::DenseVector{Complex{Float64}}, T::Freq2NoBoundaryWave2D, X::DenseMatrix{Complex{Float64}})
	(M = size(T,1)) == length(y) || throw(DimensionMismatch())
	wsize(T) == size(X) || throw(DimensionMismatch())

	myfft!(y, T.NFFT, X)
	for m in 1:M
		@inbounds y[m] *= T.internal[:x][m] * T.internal[:y][m]
	end

	isuniform(T) || had!(y, get(T.weights))

	return y
end


function Base.A_mul_B!(y::DenseVector{Complex{Float64}}, T::Freq2BoundaryWave1D, x::DenseVector{Complex{Float64}})
	(M = size(T,1)) == length(y) || throw(DimensionMismatch())
	size(T,2) == length(x) || throw(DimensionMismatch())

	xleft, xint, xright = split(x, van_moment(T))

	# Internal scaling function
	myfft!(y, T.NFFT, xint)
	had!(y, T.internal)

	# Contribution from the boundaries
	BLAS.gemv!('N', ComplexOne, T.left, xleft, ComplexOne, y)
	BLAS.gemv!('N', ComplexOne, T.right, xright, ComplexOne, y)

	isuniform(T) || had!(y, get(T.weights))

	return y
end


function Base.A_mul_B!(y::DenseVector{Complex{Float64}}, T::Freq2BoundaryWave2D, X::DenseMatrix{Complex{Float64}})
	(M = size(T,1)) == length(y) || throw(DimensionMismatch())
	(N = wsize(T)) == size(X) || throw(DimensionMismatch())

	# Internal scaling functions
	vm = van_moment(T)
	S = split(X, vm)
	myfft!(y, T.NFFT, S.II)
	for m in 1:M
		@inbounds y[m] *= T.internal[:x][m] * T.internal[:y][m]
	end

	onesp = ones(Complex{Float64}, vm)

	# ------------------------------------------------------------
	# Corners

	# LL
	A_mul_B!( T.tmpMulVec, T.left[:x], S.LL )
	had!( T.tmpMulVec, T.left[:y] )
	# y += sum(T.tmpMulVec,2) :
	BLAS.gemv!('N', ComplexOne, T.tmpMulVec, onesp, ComplexOne, y)

	# RL
	A_mul_B!( T.tmpMulVec, T.right[:x], S.RL )
	had!( T.tmpMulVec, T.left[:y] )
	BLAS.gemv!('N', ComplexOne, T.tmpMulVec, onesp, ComplexOne, y)

	# RR
	A_mul_B!( T.tmpMulVec, T.right[:x], S.RR )
	had!( T.tmpMulVec, T.right[:y] )
	BLAS.gemv!('N', ComplexOne, T.tmpMulVec, onesp, ComplexOne, y)

	# LR
	A_mul_B!( T.tmpMulVec, T.left[:x], S.LR )
	had!( T.tmpMulVec, T.right[:y] )
	BLAS.gemv!('N', ComplexOne, T.tmpMulVec, onesp, ComplexOne, y)

	# ------------------------------------------------------------
	# Sides

	# IL
	NFFT.nfft!(T.NFFTx, S.IL, T.tmpMulVec)
	scale!(T.internal[:x], T.tmpMulVec)
	had!(T.tmpMulVec, T.left[:y])
	BLAS.gemv!('N', ComplexOne, T.tmpMulVec, onesp, ComplexOne, y)

	# IR
	NFFT.nfft!(T.NFFTx, S.IR, T.tmpMulVec)
	scale!(T.internal[:x], T.tmpMulVec)
	had!(T.tmpMulVec, T.right[:y])
	BLAS.gemv!('N', ComplexOne, T.tmpMulVec, onesp, ComplexOne, y)

	# LI 
	NFFT.nfft!(T.NFFTy, S.LI, T.tmpMulVecT)
	transpose!(T.tmpMulVec, T.tmpMulVecT)
	scale!(T.internal[:y], T.tmpMulVec)
	had!(T.tmpMulVec, T.left[:x])
	BLAS.gemv!('N', ComplexOne, T.tmpMulVec, onesp, ComplexOne, y)

	# RI 
	NFFT.nfft!(T.NFFTy, S.RI, T.tmpMulVecT)
	transpose!(T.tmpMulVec, T.tmpMulVecT)
	scale!(T.internal[:y], T.tmpMulVec)
	had!(T.tmpMulVec, T.right[:x])
	BLAS.gemv!('N', ComplexOne, T.tmpMulVec, onesp, ComplexOne, y)


	isuniform(T) || had!(y, get(T.weights))

	return y
end

function Base.A_mul_B!(y::DenseVector{Complex{Float64}}, T::Freq2Wave2D, x::DenseVector{Complex{Float64}})
	size(T,2) == length(x) || throw(DimensionMismatch())

	X = reshape(x, wsize(T))
	A_mul_B!(y, T, X)

	return y
end

@compat function Base.:*(T::Freq2Wave, x::DenseArray)
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

	myfft_adjoint!(z, T.NFFT, T.tmpMulVec)

	return z
end

function Base.Ac_mul_B!(Z::DenseMatrix{Complex{Float64}}, T::Freq2NoBoundaryWave2D, v::DenseVector{Complex{Float64}})
	(M = size(T,1)) == length(v) || throw(DimensionMismatch())
	wsize(T) == size(Z) || throw(DimensionMismatch())

	for m in 1:M
		@inbounds T.tmpMulVec[m] = v[m] * conj(T.internal[:x][m]) * conj(T.internal[:y][m])
	end
	isuniform(T) || had!(T.tmpMulVec, get(T.weights))

	myfft_adjoint!(Z, T.NFFT, T.tmpMulVec)

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
	myfft_adjoint!(zint, T.NFFT, T.tmpMulVec)

	return z
end

function Base.Ac_mul_B!(Z::DenseMatrix{Complex{Float64}}, T::Freq2BoundaryWave2D, v::DenseVector{Complex{Float64}})
	(M = size(T,1)) == length(v) || throw(DimensionMismatch())
	(N = wsize(T)) == size(Z) || throw(DimensionMismatch())
	
	# As in Ac_mul_B! for Freq2BoundaryWave1D
	if isuniform(T)
		copy!(T.weigthedVec, v)
	else
		had!(T.weigthedVec, v, get(T.weights))
	end

	vm = van_moment(T)
	S = split(Z, vm)

	# Internal coefficients
	for m in 1:M
		@inbounds T.tmpMulcVec[m] = T.weigthedVec[m] * conj(T.internal[:x][m]) * conj(T.internal[:y][m])
	end
	myfft_adjoint!(S.II, T.NFFT, T.tmpMulcVec)

	# ------------------------------------------------------------
	# Left blocks: All use 'L' for the y coordinate

	# LL
	conj!(T.tmpMulVec, T.left[:y])
	scale!(T.weigthedVec, T.tmpMulVec)
	Ac_mul_B!( S.LL, T.left[:x], T.tmpMulVec )

	# RL: Reuse T.tmpMulVec
	Ac_mul_B!( S.RL, T.right[:x], T.tmpMulVec )

	# IL
	conj!(T.tmpMulcVec, T.internal[:x])
	scale!(T.tmpMulcVec, T.tmpMulVec)
	NFFT.nfft_adjoint!( T.NFFTx, T.tmpMulVec, S.IL )

	# ------------------------------------------------------------
	# Middle blocks: All use 'I' for the y coordinate

	# LI
	conj!(T.tmpMulVec, T.left[:x])
	scale!(T.weigthedVec, T.tmpMulVec)

	conj!(T.tmpMulcVec, T.internal[:y])
	scale!(T.tmpMulcVec, T.tmpMulVec)
	transpose!(T.tmpMulVecT, T.tmpMulVec)
	NFFT.nfft_adjoint!( T.NFFTy, T.tmpMulVecT, S.LI )

	# RI: Reuse T.tmpMulcVec
	conj!(T.tmpMulVec, T.right[:x])
	scale!(T.weigthedVec, T.tmpMulVec)
	scale!(T.tmpMulcVec, T.tmpMulVec)
	transpose!(T.tmpMulVecT, T.tmpMulVec)
	NFFT.nfft_adjoint!( T.NFFTy, T.tmpMulVecT, S.RI )

	# ------------------------------------------------------------
	# Right blocks: All use 'R' for the y coordinate

	# LR
	conj!(T.tmpMulVec, T.right[:y])
	scale!(T.weigthedVec, T.tmpMulVec)
	Ac_mul_B!( S.LR, T.left[:x], T.tmpMulVec )

	# RR: Reuse T.tmpMulVec
	Ac_mul_B!( S.RR, T.right[:x], T.tmpMulVec )

	# IR
	conj!(T.tmpMulcVec, T.internal[:x])
	scale!(T.tmpMulcVec, T.tmpMulVec)
	NFFT.nfft_adjoint!( T.NFFTx, T.tmpMulVec, S.IR )


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

@compat function Base.:\(T::Freq2Wave, Y::AbstractMatrix)
	length(Y) == (M = size(T,1)) || throw(DimensionMismatch())

	y = flatten_view(Y)
	x = T \ y
end

@compat function Base.:\(T::Freq2Wave, y::AbstractVector)
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

	phi = T.internal[:x] .* T.internal[:y]

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
		scale!(get(T.weights), F)
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
		unsafe_FourScaling!(phiy, T, ny, :y)
		for nx in 1:Nx
			unsafe_FourScaling!(phix, T, nx, :x)
			had!(phix, phiy)
			F[:,row_idx+=1] = phix
		end
	end

	if !isuniform(T)
		scale!(get(T.weights), F)
	end

	return F
end

@doc """
	unsafe_FourScaling!(phi, T::Freq2BoundaryWave2D, n::Int, d::Symbol)

Replace `phi` with the `n`'th "column" from dimension `d` of `T`.
"""->
function unsafe_FourScaling!(phi::Vector{Complex{Float64}}, T::Freq2BoundaryWave2D, n::Integer, d::Symbol)
	M = length(phi)
	N = 2^wscale(T)
	p = van_moment(T)

	if p < n <= N-p
		offset = div(N, 2) + 1
        D = (d == :x ? 1 : 2)
		for m in 1:M
			@inbounds phi[m] = T.internal[d][m]*cis( -twoπ*(n-offset)*T.NFFT.x[D,m] )
		end
	elseif 1 <= n <= p
		unsafe_copy!( phi, 1, T.left[d], (n-1)*M+1, M )
	else
		unsafe_copy!( phi, 1, T.right[d], (n-N+p-1)*M+1, M )
	end
end


# ------------------------------------------------------------
# Common

function Base.size(T::Freq2Wave)
	( size(T,1), size(T,2) )
end

function Base.size(T::Freq2Wave, d::Integer)
	if d == 1
        prod(T.NFFT.M)
	elseif d == 2
		prod( wsize(T) )
	else
		throw(AssertionError())
	end
end

van_moment(T::Freq2Wave1D) = hasboundary(T) ? size(T.left,2) : van_moment(T.wavename)
van_moment(T::Freq2Wave2D) = hasboundary(T) ? size(T.left[:x],2)::Int64 : van_moment(T.wavename)

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


