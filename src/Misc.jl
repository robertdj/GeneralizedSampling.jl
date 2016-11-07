@doc """
	had!(A, B) -> A

In-place Hadamard product: Replace `A` with `A.*B`.
"""->
function had!{T<:Number}(A::StridedArray{T}, B::AbstractArray{T})
	size(A) == size(B) || throw(DimensionMismatch())
	for idx in eachindex(A,B)
		@inbounds A[idx] *= B[idx]
	end
end

@doc """
	had!(A, B, C) -> A

In-place Hadamard product: Replace `A` with `B.*C`.
"""->
function had!{T<:Number}(A::StridedArray{T}, B::AbstractArray{T}, C::AbstractArray{T})
	size(A) == size(B) == size(C) || throw(DimensionMismatch())
	for idx in eachindex(A,B)
		@inbounds A[idx] = B[idx] * C[idx]
	end
end

@doc """
	hadc!(A, B) -> A

In-place Hadamard product with complex conjugation: Replace `A` with
`A.*conj(B)`.
"""->
function hadc!{T<:Number}(A::StridedArray{T}, B::AbstractArray{T})
	size(A) == size(B) || throw(DimensionMismatch())
	for idx in eachindex(A,B)
		@inbounds A[idx] *= conj(B[idx])
	end
end

@doc """
	hadc!(A, B, C) -> A

In-place Hadamard product with complex conjugation: Replace `A` with
`B.*conj(C)`.
"""->
function hadc!{T<:Number}(A::StridedArray{T}, B::AbstractArray{T}, C::AbstractArray{T})
	size(A) == size(B) == size(C) || throw(DimensionMismatch())
	for idx in eachindex(A,B)
		@inbounds A[idx] = B[idx] * conj(C[idx])
	end
end

@doc """
	yphad!(y, a, b) -> y

Replace `y` with `y + a.*b`.
"""->
function yphad!{T<:Number}(y::StridedVector{T}, a::AbstractVector{T}, b::AbstractVector{T})
	(Ny = length(y)) == length(a) == length(b) || throw(DimensionMismatch())
	for ny in 1:Ny
		@inbounds y[ny] += a[ny] * b[ny]
	end
end

@doc """
	conj!(A, B) -> A

Replace `A` with `conj(B)`.
"""->
function Base.conj!(A::StridedArray{Complex{Float64}}, B::AbstractArray{Complex{Float64}})
	size(A) == size(B) || throw(DimensionMismatch())
	for idx in eachindex(A,B)
		@inbounds A[idx] = conj(B[idx])
	end
	return A
end

@doc """
	isuniform( points::Matrix ) -> Bool

Test if the rows in the `M`-by-2 matrix `points` are on a uniform 2D grid.
"""->
function isuniform{T<:Real}( points::AbstractMatrix{T} )
	M, D = size(points)
	D == 2 || throw(DimensionMismatch())

	M <= 2 && return true

	x = sort( points[:,1] )
	uniquex = unique(x)
	Mx = length(uniquex)

	y = sort( points[:,2] )
	uniquey = unique(y)
	My = length(uniquey)

	if isuniform(uniquex) && isuniform(uniquey) && Mx*My == M
		return true
	else
		return false
	end
end


@doc """
	grid(M, scale) -> Vector

Equidistant 1D points centered around the origin.
With even `M` the grid has one extra point on the negative values.

The points are scaled by `scale` which by default is 1.
"""->
function grid(M::Int, grid_dist::Real=1.0)
	M >= 2 || throw(AssertionError())
	grid_dist > 0 || throw(DomainError())

	startx = -div(M,2)
	endx = (isodd(M) ? -startx : -startx-1)

	x = collect( float(startx):float(endx) )

	if grid_dist != one(grid_dist)
		scale!(grid_dist, x)
	end

	return x
end

@doc """
	grid( (Mx,My), scale) -> Matrix

2D points on an `Mx`-by-`My` grid centered around the origin.

The points are scaled by `scale` which by default is 1.
"""->
function grid( M::Tuple{Integer,Integer}, grid_dist::Real=1.0)
	minimum(M) >= 2 || throw(AssertionError())
	grid_dist > 0 || throw(DomainError())

	# The points are sorted by the x coordinate
	gridx = grid(M[1], grid_dist)
    x = repmat(gridx, M[2])

	gridy = grid(M[2], grid_dist)
    y = kron(gridy, ones(M[1]))

	return hcat(x, y)
end


@doc """
	weights(xi, bandwidth) -> Vector

Compute weights for sampling points `xi`.

- For 1D points `xi` must be a vector. The function is faster if `xi` is sorted.
- For 2D points `xi` must be a matrix with 2 columns.
"""->
function weights(xi::AbstractVector{Float64}, bandwidth::Real)
	bandwidth > 0 || throw(DomainError())
	maxabs(xi) <= bandwidth || throw(AssertionError())
	(Nx = length(xi)) >= 2 || throw(AssertionError())

	is_xi_sorted = issorted(xi)
	if !is_xi_sorted
		sort_idx = sortperm(xi)
		inv_sort_idx = invperm(sort_idx)
		xi = xi[sort_idx]
	end

	mu = Array{Float64}(Nx)

	# Boundary cases
	L = xi[Nx] - 2*bandwidth
	mu[1] = 0.5*(xi[2] - L)

	U = xi[1] + 2*bandwidth
	mu[Nx] = 0.5*(U - xi[Nx])

	# Non-boundary cases
	for n in 2:Nx-1
		@inbounds mu[n] = 0.5*(xi[n+1] - xi[n-1])
	end

	!is_xi_sorted && permute!(mu, inv_sort_idx)

	return mu 
end

function weights(xi::DenseMatrix{Float64}, bandwidth::Real)
	bandwidth > 0 || throw(DomainError())
	size(xi,2) == 2 || throw(DimensionMismatch())
	size(xi,1) >= 2 || throw(DimensionMismatch())
	maxabs(xi) <= bandwidth || throw(AssertionError())

	x = view(xi, :, 1)
	y = view(xi, :, 2)
	rw = [-bandwidth; bandwidth; -bandwidth; bandwidth]
	voronoiarea(x, y, rw)
end


@doc """
	density(xi, K) -> δ

Compute the density of the sampling points `xi` with bandwidth `K`:

- In 1D, `xi` must be a vector and δ is the largest difference between two numerically consecutive points in `xi`.
- In 2D, `xi` must be an `M-by-2` matrix and the density is the smallest number δ that allows a covering of the "bandwidth area" with circles centered in the points of `xi` and radii δ.
Currently, the "bandwidth area" is a square centered at the origin and with sidelength `2*K`. 

To ensure numerical stability in computations with the associated change of basis matrix, `xi` must contain both negative and positive elements. 
"""->
function density(xi::AbstractVector{Float64}, bandwidth::Real)
	(minx = minimum(xi)) < 0 || throw(AssertionError())
	0 < (maxx = maximum(xi)) <= bandwidth || throw(AssertionError())
	(N = length(xi)) >= 2 || throw(DimensionMismatch())

	# Boundary cases
	lower_boundary = maxx - 2*bandwidth
	upper_boundary = minx + 2*bandwidth
	density = max(minx-lower_boundary, upper_boundary-maxx)

	# Non-boundary cases
	xi_sorted = sort(xi)
	for n in 2:N
		@inbounds diff = xi_sorted[n] - xi_sorted[n-1]
		density = max(density, diff)
	end

	return density
end

function density(xi::AbstractMatrix{Float64}, bandwidth::Real)
	bandwidth > 0 || throw(DomainError())
	maxabs(xi) <= bandwidth || throw(AssertionError())
	M, dim = size(xi)
	M >= 2 || throw(DimensionMismatch())
	dim == 2 || throw(DimensionMismatch())

	x = view(xi, :, 1)
	y = view(xi, :, 2)
	rw = [-bandwidth; bandwidth; -bandwidth; bandwidth]
	density(x, y, rw)
end

@doc """
	frac(x)
	frac!(x)

The fractional part of `x` as a number in [-0.5, 0.5).
"""->
function frac(x::DenseArray{Float64})
	y = copy(x)
	frac!(y)

	return y
end

function frac!(x::DenseArray{Float64})
	for idx in eachindex(x)
		@inbounds x[idx] -= round( x[idx] )
	end
end


@doc """
	isdaubechies(wavename::AbstractString) -> Bool

Return `true` if `wavename` is of the form `dbN`, where `N` is an integer.
"""->
function isdaubechies(wavename::AbstractString)
	ishaar(wavename) && return true

	prefix = lowercase(wavename[1:2])
	return prefix == "db" && isnumber(wavename[3:end])
end

@doc """
	ishaar(wavename::AbstractString) -> Bool

Return `true` if `wavename` is "haar" or "db1".
"""->
function ishaar(wavename::AbstractString)
	lowername = lowercase(wavename)

	if lowername == "haar" || lowername == "db1"
		return true
	else
		return false
	end
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


@doc """
	split(x::Vector, border) -> L, I, R

Split `x` into 3 slices:
Left, internal and right, where left and right are the `border` outmost entries.
"""->
function Base.split(x::DenseVector, border::Integer)
	border >= 1 || throw(DomainError())
	(N = length(x)) > 2*border || throw(AssertionError())

	L = view(x, 1:border)
	I = view(x, border+1:N-border)
	R = view(x, N-border+1:N)

	return L, I, R
end


type SplitMatrix{T, A<:AbstractMatrix}
	LL::A
	IL::A
	RL::A
	LI::A
	II::A
	RI::A
	LR::A
	IR::A
	RR::A

	SplitMatrix(LL::AbstractMatrix{T}, IL::AbstractMatrix{T},
	RL::AbstractMatrix{T}, LI::AbstractMatrix{T}, II::AbstractMatrix{T},
	RI::AbstractMatrix{T}, LR::AbstractMatrix{T}, IR::AbstractMatrix{T},
	RR::AbstractMatrix{T}) = new(LL, IL, RL, LI, II, RI, LR, IR, RR)
end
SplitMatrix(LL, IL, RL, LI, II, RI, LR, IR, RR) =
SplitMatrix{eltype(LL), typeof(LL)}(LL, IL, RL, LI, II, RI, LR, IR, RR)

@doc """
	split(A::Matrix, border) -> SplitMatrix

Split `A` into 9 parts, where the outer parts are `border` wide/high.
`LL` is e.g. `border`-by-`border`.

	 ________________ 
	|    |      |    |
	| LL |  LI  | LR |
	|____|______|____|
	|    |      |    |
	| IL |  II  | IR |
	|____|______|____|
	|    |      |    |
	| RL |  RI  | RR |
	|____|______|____|
"""->
function Base.split{T}(A::DenseMatrix{T}, border::Integer)
	border >= 2 || throw(DomainError())
	N = size(A)
	minimum(N) > 2*border || throw(AssertionError())

	Lidx = 1:border
	I1idx = border+1:N[1]-border
	R1idx = N[1]-border+1:N[1]
	I2idx = border+1:N[2]-border
	R2idx = N[2]-border+1:N[2]

	LL = view(A, Lidx, Lidx)
	IL = view(A, I1idx, Lidx)
	RL = view(A, R1idx, Lidx)

	LI = view(A, Lidx, I2idx)
	II = view(A, I1idx, I2idx)
	RI = view(A, R1idx, I2idx)

	LR = view(A, Lidx, R2idx)
	IR = view(A, I1idx, R2idx)
	RR = view(A, R1idx, R2idx)

	SplitMatrix{T, typeof(LL)}( LL, IL, RL, LI, II, RI, LR, IR, RR )
end

