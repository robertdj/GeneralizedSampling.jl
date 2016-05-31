# ------------------------------------------------------------
# 1D change of basis matrix

@doc """
	freq2wave(samples, wavename::String, J::Int, B; ...)

Make change of basis for switching between frequency responses and wavelets.

- `samples` are the sampling locations as a vector for 1D and 2-by-M / M-by-2 matrix for 2D.
- `wave` is the name of the wavelet; see documentation for possibilities.
- `J` is the scale of the wavelet transform to reconstruct.
- `B` is the bandwidth of the samples; only needed if `samples` are non-uniform.
- Optional arguments (if needed) are passed to the functions computing Fourier transforms.
"""->
function freq2wave(samples::DenseVector, wavename::AbstractString, J::Int, B::Float64=NaN; args...)
	vm = van_moment(wavename)
	Nint = 2^J
	@assert Nint >= 2*vm-1 "Scale it not large enough for this wavelet"
	
	# TODO: Good condition?
	if Nint >= length(samples)
		warn("The scale is high compared to the number of samples")
	end

	# Weights for non-uniform samples
	if isuniform(samples)
		const W = Nullable{ Vector{Complex{Float64}} }()
	else
		isnan(B) && error("Samples are not uniform; supply bandwidth")
		Nint > 2*B || warn("The scale is high compared to the bandwidth")
		W = sqrt(weights(samples, B))
		const W = Nullable(complex( W ))
	end

	# Fourier transform of the internal scaling function
	const FT = FourScalingFunc( samples, wavename, J; args... )

	# Diagonal matrix used in 'multiplication' with NFFT
	const diag = FT .* cis(-pi*samples)

	# The number of internal wavelets in reconstruction
	if hasboundary(wavename)
		Nint -= 2*van_moment(wavename)
		Nint <= 0 && error("Too few wavelets: Boundary functions overlap")
	end

	# NFFTPlans: Frequencies must be in the torus [-1/2, 1/2)
	# TODO: Should window width m and oversampling factor sigma be changed for higher precision?
	xi = samples*2.0^(-J)
	frac!(xi)
	const p = NFFTPlan(xi, Nint)

	# Wavelets w/o boundary
	if !hasboundary(wavename)
		return Freq2NoBoundaryWave(samples, FT, W, J, wavename, diag, p)
	else
		# TODO: Extend FourScalingFunc
		left = FourDaubScaling(samples, vm, 'L', J; args...)'
		# TODO: The right scaling functions are indexed "opposite": The
		# function closest to the boundary is first and the function
		# furthest from the boundary is last.
		right = FourDaubScaling(samples, vm, 'R', J; args...)'
		phase_shift = cis( -twoπ*samples )
		broadcast!(*, right, right, phase_shift)

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


# ------------------------------------------------------------
# Basic operations for 1D Freq2Wave

function Base.collect(T::Freq2NoBoundaryWave{1})
	M, N = size(T)

	F = Array{Complex{Float64}}(M, N)
	for n = 1:N
		for m = 1:M
			@inbounds F[m,n] = T.internal[m]*cis( -2*pi*T.samples[m]*(n-1)/N )
			#= @inbounds F[m,n] = T.internal[m]*cis( -2*pi*T.NFFT.x[m]*(n-1) ) =#
		end
	end

	if !isuniform(T)
		broadcast!(*, F, F, get(T.weights))
	end

	return F
end

function Base.collect(T::Freq2BoundaryWave{1})
	M, N = size(T)
	F = Array{Complex{Float64}}(M, N)

	vm = van_moment(T)

	# Left boundary
	F[:,1:vm] = T.left

	# Internal function
	for n = vm+1:N-vm
		for m = 1:M
			@inbounds F[m,n] = T.internal[m]*cis( -2*pi*T.samples[m]*(n-1)/N )
			#= @inbounds F[m,n] = W[m]*T.internal[m]*cis( -2*pi*T.NFFT.x[m]*(n-1) ) =#
		end
	end

	# Right boundary
	F[:,N-vm+1:N] = T.right

	if !isuniform(T)
		broadcast!(*, F, F, get(T.weights))
	end

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
function wsize(T::Freq2Wave{1})
	hasboundary(T) ? (2^wscale(T),) : (T.NFFT.N,)
end

function wsize(T::Freq2Wave{2})
	if hasboundary(T)
		N = 2^wscale(T)
		return (N, N)
	else
		return T.NFFT.N
	end
end


function freq2wave(samples::DenseMatrix, wavename::AbstractString, J::Int, B::Float64=NaN)
	vm = van_moment(wavename)
	Nint = 2^J
	@assert Nint >= 2*vm-1 "Scale it not large enough for this wavelet"
	M = size(samples, 1)
	@assert size(samples,2) == 2

	# TODO: Good condition?
	if Nint >= M
		warn("The scale is high compared to the number of samples")
	end

	# Weights for non-uniform samples
	if isuniform(samples)
		const W = Nullable{ Vector{Complex{Float64}} }()
	else
		isnan(B) && error("Samples are not uniform; supply bandwidth")
		Nint > 2*B || warn("The scale is high compared to the bandwidth")
		W = sqrt(weights(samples, B))
		const W = Nullable(complex( W ))
	end

	# Fourier transform of the internal scaling function
	const FT = FourScalingFunc( samples, wavename, J )

	# Diagonal matrix for multiplication with internal scaling function
	const diag = FT .* cis( -pi*samples )

	# The number of internal wavelets in reconstruction
	if hasboundary(wavename)
		Nint -= 2*van_moment(wavename)
		Nint <= 0 && error("Too few wavelets: Boundary functions overlap")
	end

	# NFFTPlans: Frequencies must be in the torus [-1/2, 1/2)^2
	xi = samples'
	scale!(xi, 2.0^(-J))
	frac!(xi)
	const p = NFFTPlan(xi, (Nint,Nint))

	# Wavelets w/o boundary
	if !hasboundary(wavename)
		return Freq2NoBoundaryWave(samples, FT, W, J, wavename, diag, p)
	else
		samplesx = slice(samples, :, 1)
		samplesy = slice(samples, :, 2)

		vm = van_moment(wavename)
		leftX = FourDaubScaling(samplesx, vm, 'L', J)
		leftY = FourDaubScaling(samplesy, vm, 'L', J)
		const left = cat(3, leftX', leftY' )

		# TODO: In loop?
		rightX = FourDaubScaling(samplesx, vm, 'R', J)'
		phase_shift = cis( -twoπ*samplesx )
		broadcast!(*, rightX, rightX, phase_shift)

		rightY = FourDaubScaling(samplesy, vm, 'R', J)'
		phase_shift = cis( -twoπ*samplesy )
		broadcast!(*, rightY, rightY, phase_shift)
		const right = cat(3, rightX, rightY )

		return Freq2BoundaryWave(samples, FT, W, J, wavename, diag, p, left, right)
	end
end


function Base.A_mul_B!{D}(y::DenseVector{Complex{Float64}}, T::Freq2NoBoundaryWave{D}, X::DenseArray{Complex{Float64},D})
	@assert size(T,1) == length(y)
	@assert wsize(T) == size(X)

	NFFT.nfft!(T.NFFT, X, y)
	for d in 1:D
		had!(y, T.diag[:,d])
	end

	isuniform(T) || had!(y, get(T.weights))

	return y
end


function Base.A_mul_B!(y::DenseVector{Complex{Float64}}, T::Freq2BoundaryWave{1}, x::DenseVector{Complex{Float64}})
	@assert size(T,1) == length(y)
	@assert (Nx = size(T,2)) == length(x)

	xleft, xint, xright = split(x, van_moment(T))

	# Internal scaling function
	NFFT.nfft!(T.NFFT, xint, y)
	had!(y, T.diag)

	# Update y with boundary contributions
	BLAS.gemv!('N', ComplexOne, T.left, xleft, ComplexOne, y)
	BLAS.gemv!('N', ComplexOne, T.right, xright, ComplexOne, y)

	isuniform(T) || had!(y, get(T.weights))

	return y
end

function Base.A_mul_B!(y::DenseVector{Complex{Float64}}, T::Freq2BoundaryWave{2}, X::DenseMatrix{Complex{Float64}})
	@assert size(T,1) == length(y)
	@assert wsize(T) == size(X)
	
	const vm = van_moment(T)
	const S = split(X, vm)

	# Internal scaling function
	NFFT.nfft!(T.NFFT, S.II, y)
	had!(y, T.diag[:,1])
	had!(y, T.diag[:,2])

	# The non-internal contribution have a multiplication from each
	# dimension of X
	innery = similar(y)

	# TODO: Am I using the correct order for the right boundary?

	# Update y with border and side contributions
	# Fourier transform of the internal functions using 1D NFFT
	# TODO: Include these "1D" NFFT's in NFFT package?
	const p1 = NFFTPlan( vec(T.NFFT.x[1,:]), T.NFFT.N[1] )
	const p2 = NFFTPlan( vec(T.NFFT.x[2,:]), T.NFFT.N[2] )
	for k in 1:vm
		# LL
		A_mul_B!( innery, T.left[:,:,1], S.LL[:,k] )
		yphad!(y, T.left[:,k,2], innery)
		# LR
		A_mul_B!( innery, T.left[:,:,1], S.LR[:,k] )
		yphad!(y, T.right[:,k,2], innery)
		# RL
		A_mul_B!( innery, T.right[:,:,1], S.RL[:,k] )
		yphad!(y, T.left[:,k,2], innery)
		# RR
		A_mul_B!( innery, T.right[:,:,1], S.RR[:,k] )
		yphad!(y, T.right[:,k,2], innery)

		# LI
		NFFT.nfft!(p2, vec(S.LI[k,:]), innery) # DFT part of internal
		had!(innery, T.diag[:,2]) # Fourier transform of internal
		yphad!(y, T.left[:,k,1], innery) # Contribution from left

		# IL
		NFFT.nfft!(p1, S.IL[:,k], innery)
		had!(innery, T.diag[:,1])
		yphad!(y, T.left[:,k,2], innery)

		# RI
		NFFT.nfft!(p2, vec(S.RI[k,:]), innery)
		had!(innery, T.diag[:,2])
		yphad!(y, T.right[:,k,1], innery)

		# IR
		NFFT.nfft!(p1, S.IR[:,k], innery)
		had!(innery, T.diag[:,1])
		yphad!(y, T.right[:,k,2], innery)
	end

	isuniform(T) || had!(y, get(T.weights))

	return y
end

function Base.A_mul_B!(y::DenseVector{Complex{Float64}}, T::Freq2Wave{2}, x::DenseVector{Complex{Float64}})
	@assert size(T,2) == length(x)

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


function Base.Ac_mul_B!{D}(Z::DenseArray{Complex{Float64},D}, T::Freq2NoBoundaryWave{D}, v::DenseVector{Complex{Float64}})
	@assert size(T,1) == length(v)
	@assert wsize(T) == size(Z)

	Tdiag = vec(prod(T.diag, 2))
	conj!(Tdiag)
	isuniform(T) || had!(Tdiag, get(T.weights))
	had!(Tdiag, v)

	NFFT.nfft_adjoint!(T.NFFT, Tdiag, Z)

	return Z
end

function Base.Ac_mul_B!(z::DenseVector{Complex{Float64}}, T::Freq2BoundaryWave{1}, v::DenseVector{Complex{Float64}})
	@assert size(T,1) == length(v)
	@assert (Nz = size(T,2)) == length(z)

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
	@assert size(T,1) == length(v)
	@assert wsize(T) == size(Z)
	
	const vm = van_moment(T)
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
	const Nint = size(S.IL,1)
	sidecol = Array{Complex{Float64}}(Nint)

	const p1 = NFFTPlan( vec(T.NFFT.x[1,:]), T.NFFT.N[1] )
	const p2 = NFFTPlan( vec(T.NFFT.x[2,:]), T.NFFT.N[2] )
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
	@assert size(T,2) == length(z)

	Z = reshape_view(z, wsize(T))
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
	@assert length(Y) == (M = size(T,1))

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

	#print("Solution via conjugate gradients... ")
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
function Base.collect(T::Freq2NoBoundaryWave{2})
	const M = size(T,1)
	const Nx, Ny = wsize(T)
	# TODO: Check if the matrix fits in memory
	F = Array{Complex{Float64}}(M, size(T,2))

	phi = prod(T.internal, 2)

	idx = 0
	for ny = 0:Ny-1
		for nx = 0:Nx-1
			idx += 1
			for m = 1:M
				@inbounds F[m,idx] = phi[m]*cis( -twoπ*(nx*T.NFFT.x[1,m] + ny*T.NFFT.x[2,m]) )
			end
		end
	end

	if !isuniform(T)
		broadcast!( *, F, F, get(T.weights) )
	end

	return F
end

function Base.collect(T::Freq2BoundaryWave{2})
	const M = size(T,1)
	const Nx, Ny = wsize(T)
	# TODO: Check if the matrix fits in memory
	F = Array{Complex{Float64}}(M, size(T,2))

	phix = Array{Complex{Float64}}(M)
	phiy = similar(phix)

	idx = 0
	for ny = 0:Ny-1
		unsafe_FourScaling!(phiy, T, ny, 2)
		for nx = 0:Nx-1
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
function unsafe_FourScaling!(phi::Vector{Complex{Float64}}, T::Freq2BoundaryWave{2}, n::Integer, d::Integer)
	M = length(phi)
	N = wsize(T)[d]
	p = van_moment(T)

	if p <= n < N-p
		for m in 1:M
			@inbounds phi[m] = T.internal[m,d]*cis( -twoπ*n*T.samples[m,d]/N )
		end
	elseif 0 <= n < p
		# source offset = dimension offset + column offset
		so = (d-1)*M*p + n*M + 1
		unsafe_copy!( phi, 1, T.left, so, M ) 
	else
		so = (d-1)*M*p + (n-N+p)*M + 1
		unsafe_copy!( phi, 1, T.right, so, M ) 
	end
end


# ------------------------------------------------------------
# Common

function Base.size(T::Freq2Wave, ::Type{Val{1}})
	size(T.samples, 1)
end

function Base.size{D}(T::Freq2Wave{D}, ::Type{Val{2}})
	2^(D*wscale(T))
end

function van_moment(T::Freq2Wave)
	if hasboundary(T)
		return size(T.left,2)
	else
		return van_moment(T.wavename)
	end
end

function Base.eltype(::Freq2Wave)
	return Complex{Float64}
end

function Base.show{D}(io::IO, T::Freq2Wave{D})
	println(io, D, "D change of basis matrix")

	isuniform(T) ?  U = " " : U = " non-"
	M = size(T,1)
	println(io, "From: ", M, U, "uniform frequency samples")

	D == 1 ? N = size(T,2) : N = wsize(T)
	print(io, "To: ", N, " ", T.wavename, " wavelets")
end

