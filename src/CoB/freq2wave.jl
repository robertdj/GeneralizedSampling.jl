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
function freq2wave(samples::DenseVector, wavename::AbstractString, J::Int; B::Float64=NaN)
	# TODO: Warning if J is too big
	# TODO: IS it better with B as a non-optional argument?

	# Fourier transform of the internal scaling function
	const FT = FourScalingFunc( samples, wavename, J )

	# Diagonal matrix used in 'multiplication' with NFFT
	# TODO: Multiply with weights (in the 1D case only)?
	const diag = FT .* cis(-pi*samples)

	# Weights for non-uniform samples
	if isuniform(samples)
		const W = Nullable{ Vector{Complex{Float64}} }()
	else
		isnan(B) && error("Samples are not uniform; supply bandwidth")
		W = sqrt(weights(samples, B))
		const W = Nullable(complex( W ))
	end

	# The number of internal wavelets in reconstruction
	Nint = 2^J
	if hasboundary(wavename)
		Nint -= 2*van_moment(wavename)
		Nint <= 0 && error("Too few wavelets: Boundary functions overlap")
	end

	# NFFTPlans: Frequencies must be in the torus [-1/2, 1/2)
	# TODO: Should window width m and oversampling factor sigma be changed for higher precision?
	xi = samples/2^J
	frac!(xi)
	const p = NFFTPlan(xi, Nint)

	# Wavelets w/o boundary
	if !hasboundary(wavename)
		return Freq2NoBoundaryWave(samples, FT, W, J, wavename, diag, p)
	else
		F = waveparse( wavename, true )
		const left = FourDaubScaling(samples, F, 'L')
		const right = FourDaubScaling(samples, F, 'R')

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

	ishaar(lowername) && return false

	if isdaubechies(lowername)
		return true
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

	if !isuniform(T)
		broadcast!(*, F, F, get(T.weights))
	end

	return F
end

function Base.collect(T::Freq2BoundaryWave{1})
	M, N = size(T)
	F = Array(Complex{Float64}, M, N)

	vm = van_moment(T)

	# Left boundary
	F[:,1:vm] = T.left

	# Internal function
	for n = vm+1:N-vm
		for m = 1:M
			@inbounds F[m,n] = T.FT[m]*cis( -2*pi*T.samples[m]*(n-1)/N )
			#= @inbounds F[m,n] = W[m]*T.FT[m]*cis( -2*pi*T.NFFT.x[m]*(n-1) ) =#
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
function wsize(T::Freq2Wave)
	if hasboundary(T)
		N = 2^wscale(T)
		return dim(T) == 1 ? (N,) : (N, N)
	else
		return T.NFFT.N
	end
end


function freq2wave(samples::DenseMatrix, wavename::AbstractString, J::Int; B::Float64=NaN)
	# TODO: Samples as 2-by-M? Faster access to each point
	M = size(samples, 1)
	@assert size(samples,2) == 2
	# TODO: Warning if J is too big

	# Fourier transform of the internal scaling function
	const FT = FourScalingFunc( samples, wavename, J )

	# Diagonal matrix for multiplication with internal scaling function
	samplesx = slice(samples, :, 1)
	samplesy = slice(samples, :, 2)
	#= const diag = vec(prod(FT,2)) .* cis( -pi*(samplesx + samplesy) ) =#
	const diag = FT .* cis( -pi*samples )

	# Weights for non-uniform samples
	if isuniform(samples)
		const W = Nullable{ Vector{Complex{Float64}} }()
	else
		isnan(B) && error("Samples are not uniform; supply bandwidth")
		W = sqrt(weights(samples, B))
		const W = Nullable(complex( W ))
	end

	# The number of internal wavelets in reconstruction
	Nint = 2^J
	if hasboundary(wavename)
		Nint -= 2*van_moment(wavename)
		Nint <= 0 && error("Too few wavelets: Boundary functions overlap")
	end

	# NFFTPlans: Frequencies must be in the torus [-1/2, 1/2)^2
	xi = samples'
	scale!(xi, 1/2^J)
	frac!(xi)
	const p = NFFTPlan(xi, (Nint,Nint))

	# Wavelets w/o boundary
	if !hasboundary(wavename)
		return Freq2NoBoundaryWave(samples, FT, W, J, wavename, diag, p)
	else
		F = waveparse( wavename, true )
		const left = cat(3, FourDaubScaling(samplesx, F, 'L'), FourDaubScaling(samplesy, F, 'L') )
		const right = cat(3, FourDaubScaling(samplesx, F, 'R'), FourDaubScaling(samplesy, F, 'R') )

		return Freq2BoundaryWave(samples, FT, W, J, wavename, diag, p, left, right)
	end
end


function Base.A_mul_B!{D}(y::DenseVector{Complex{Float64}}, T::Freq2NoBoundaryWave{D}, X::DenseArray{Complex{Float64},D})
	@assert size(T,1) == length(y)
	@assert wsize(T) == size(X)

	NFFT.nfft!(T.NFFT, X, y)
	for d = 1:D
		had!(y, T.diag[:,d])
	end

	!isuniform(T) && had!(y, get(T.weights))

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
	const Cone = one(Complex{Float64})
	BLAS.gemv!('N', Cone, T.left, xleft, Cone, y)
	BLAS.gemv!('N', Cone, T.right, xright, Cone, y)

	!isuniform(T) && had!(y, get(T.weights))

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

	!isuniform(T) && had!(y, get(T.weights))

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

	y = Array(Complex{Float64}, size(T,1))
	A_mul_B!(y, T, x)

	return y
end


function Base.Ac_mul_B!{D}(Z::DenseArray{Complex{Float64},D}, T::Freq2NoBoundaryWave{D}, v::DenseVector{Complex{Float64}})
	@assert size(T,1) == length(v)
	@assert wsize(T) == size(Z)

	Tdiag = vec(prod(T.diag, 2))
	conj!(Tdiag)
	!isuniform(T) && had!(Tdiag, get(T.weights))
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
	innerv = vec(prod(T.diag, 2))
	conj!(innerv)
	had!(innerv, weigthedV)
	NFFT.nfft_adjoint!(T.NFFT, innerv, S.II)

	# Temporary arrays for holding results of the "inner" multiplication.
	# Assuming that T.NFFT.N[1] == T.NFFT.N[2] only one cornercol and
	# sidecol is needed
	cornercol = Array(Complex{Float64}, vm)
	const Nint = size(S.IL,1)
	sidecol = Array(Complex{Float64}, Nint)

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

	z = Array(Complex{Float64}, size(T,2))
	Ac_mul_B!(z, T, v)

	return z
end

function Base.(:(\))(T::Freq2Wave, y::AbstractVector)
	if !isa(y, Array{Complex{Float64}})
		y = map(Complex{Float64}, y)
	end

	# Non-uniform samples: Scale observations
	if !isuniform(T)
		y .*= get(T.weights)
	end

	# TODO: Exact solution for 2D uniform samples?

	println("Solution via conjugate gradients... ")
	# TODO: Better initial guess?
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
	F = Array(Complex{Float64}, M, size(T,2))

	#= xsample = T.samples[:,1]/Nx =#
	#= ysample = T.samples[:,2]/Ny =#

	#= translx = [0:Nx-1;] - Nx/2 =#
	#= transly = [0:Ny-1;] - Ny/2 =#
	idx = 0
	for ny = 1:Ny
		#= yidx = transly[ny] =#
		for nx = 1:Nx
			#= xidx = translx[nx] =#
			idx += 1
			for m = 1:M
				#= @inbounds F[m,idx] = T.diag[m]*cis( -2*pi*(xidx*xsample[m] + yidx*ysample[m]) ) =#
				#= @inbounds F[m,idx] = T.FT[m,1]*T.FT[m,2]*cis( -2*pi*((nx-1)*T.samples[m,1]/Nx + (ny-1)*T.samples[m,2]/Ny) ) =#
				# TODO: T.FT -> T.FT'
				@inbounds F[m,idx] = T.FT[m,1]*T.FT[m,2]*cis( -2*pi*((nx-1)*T.NFFT.x[1,m] + (ny-1)*T.NFFT.x[2,m]) )
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
	F = Array(Complex{Float64}, M, size(T,2))

	idx = 0
	for ny = 1:Ny
		for nx = 1:Nx
			idx += 1
			for m = 1:M
				@inbounds F[m,idx] = FourScaling(T, m, (nx,ny))
			end
		end
	end

	if !isuniform(T)
		broadcast!( *, F, F, get(T.weights) )
	end

	return F
end

@doc """
	FourScaling(T::Freq2Wave{2}, d::Integer, m::Integer, n::Integer) -> Number

Fourier transform of the `n`'th basis function at the `d`'th coordinate at the `m`'th frequency.
"""->
function FourScaling(T::Freq2BoundaryWave{2}, m::Integer, n::Integer, d::Integer)
	# TODO: Prefix this function with unsafe_ and skip all checks
	# TODO: The majority of the time is spent in the following 4 lines
	#= const M = size(T,1) =#
	#= @assert 1 <= m <= M =#
	const N = wsize(T)[d]
	#= @assert 1 <= n <= N =#
	const vm = van_moment(T)

	if vm < n <= N-vm
		@inbounds y = T.FT[m,d]*cis( -2*pi*(n-1)*T.samples[m,d]/N )
	elseif 1 <= n <= vm
		@inbounds y = T.left[m, n, d]
	elseif N-vm < n <= N
		@inbounds y = T.right[m, n-N+vm, d]
	end

	return y
end

function FourScaling(T::Freq2BoundaryWave{2}, m::Integer, n::Tuple{Integer,Integer})
	FourScaling(T,m,n[1],1) * FourScaling(T,m,n[2],2)
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

