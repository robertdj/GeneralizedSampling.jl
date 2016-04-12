@doc """
	had!(A, B)

In-place Hadamard product/element-wise matrix multiplication; the first argument is modified.
"""->
function had!{T<:Number}(A::Array{T}, B::Array{T})
	@assert size(A) == size(B)
	for idx in eachindex(A)
		@inbounds A[idx] *= B[idx]
	end
end


function isuniform( points::Matrix )
	M, D = size(points)
	@assert D == 2

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
function grid(M::Int, grid_dist=1)
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
function grid( M::Tuple{Integer,Integer}, grid_dist=1)
	# The points are sorted by the x coordinate
	gridx = grid(M[1], grid_dist)
	x = kron(gridx, ones(M[2]))

	gridy = grid(M[2], grid_dist)
	y = repmat(gridy, M[1])

	return hcat(x, y)
end


@doc """
	weights(xi, bandwidth)

Compute weights for sampling points `xi`.

- For 1D points `xi` must be a *sorted* vector.
- For 2D points `xi` must be a matrix with 2 columns.
"""->
function weights(xi::Vector, bandwidth::Real)
	# TODO: Remove this assumption?
	@assert issorted(xi)

	# TODO: Check that bandwidth is sufficiently large

	N = length(xi)
	mu = Array(Float64, N)

	# Boundary cases
	L = xi[N] - 2*bandwidth
	mu[1] = 0.5*(xi[2] - L)

	U = xi[1] + 2*bandwidth
	mu[N] = 0.5*(U - xi[N])

	# Non-boundary cases
	for n = 2:N-1
		@inbounds mu[n] = 0.5*(xi[n+1] - xi[n-1])
	end

	return mu
end

function weights(xi::Matrix, bandwidth::Real)
	# TODO: Which dim is 2?
	@assert size(xi,2) == 2
	@assert maximum(abs(xi)) <= bandwidth

	voronoiarea(xi[:,1], xi[:,2]; rw=[-bandwidth; bandwidth; -bandwidth; bandwidth])
end


@doc """
	density(xi, K)

Compute the density of the sampling points `xi` with bandwidth `K`:

- In 1D, `xi` must be a vector of *sorted* points and the density is the largest difference between two consecutive points in `xi`.
- In 2D, `xi` must be an `M-by-2` matrix and the density is the smallest number δ that allows a covering of the "bandwidth area" with circles centered in the points of `xi` and radii δ.
Currently, the "bandwidth area" is a square centered at the origin and with sidelength `2*K`. 
"""->
function density(xi::Vector, bandwidth::Number)
	@assert issorted(xi)

	N = length(xi)
	if xi[1] < -bandwidth || xi[N] > bandwidth
		error("Bandwidth is not sufficiently big")
	end

	# Boundary cases
	lower_boundary = xi[N] - 2*bandwidth
	upper_boundary = xi[1] + 2*bandwidth
	density = max(xi[1] - lower_boundary, upper_boundary - xi[N])

	# Non-boundary cases
	for n = 2:N
		@inbounds diff = xi[n] - xi[n-1]
		if diff > density
			density = diff
		end
	end

	return density
end

function density(xi::Matrix{Float64}, bandwidth::Real)
	M, dim = size(xi)
	@assert dim == 2

	# Compute Voronoi tesselation
	D = deldir(xi[:,1], xi[:,2]; rw=[-bandwidth; bandwidth; -bandwidth; bandwidth])

	# Corners of Voronoi cells
	x1 = D.vorsgs[:x1]
	y1 = D.vorsgs[:y1]
	x2 = D.vorsgs[:x2]
	y2 = D.vorsgs[:y2]

	# Edge-sampling point relation
	ind = D.vorsgs[:ind1]

	# Compute the distance from each xi to the corners of its Voronoi cell
	density = 0.0
	for n = 1:length(ind)
		idx = ind[n]
		# Distance^2 from end points of Voronoi edge to one of the xi's 
		# with this edge
		diff1 = (x1[n] - xi[idx,1])^2 + (y1[n] - xi[idx,2])^2
		diff2 = (x2[n] - xi[idx,1])^2 + (y2[n] - xi[idx,2])^2

		diff = max(diff1, diff2)
		if diff > density
			density = diff
		end
	end

	return sqrt(density)
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
	lowername = lowercase(wavename)

	if ishaar(lowername)
		return true
	end

	prefix = lowername[1:2]
	N = parse(lowername[3:end])
	return prefix == "db" && isa( N, Integer )
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
	split(x::Vector, border::Int) -> SubVector, SubVector, SubVector

Split `x` into 3 parts:
Left, internal and right, where left and right are `border` outmost entries.
"""->
function split(x::DenseVector, border::Int)
	N = length(x)

	L = slice(x, 1:border)
	I = slice(x, border+1:N-border)
	R = slice(x, N-border+1:N)

	return L, I, R
end

type SplitMatrix
	LL::AbstractMatrix
	LI::AbstractMatrix
	LR::AbstractMatrix
	IL::AbstractMatrix
	II::AbstractMatrix
	IR::AbstractMatrix
	RL::AbstractMatrix
	RI::AbstractMatrix
	RR::AbstractMatrix
end

@doc """
	split(A::Matrix, border) -> SplitMatrix

Split `A` into 9 parts:
4 corners, 4 sides and the internal part. `border` is the width of the
boundary.

With both the horizontal and the vertical part divided in `L`eft,
`I`nternal and `R`ight, the parts are

	 ________________ 
	| LL |  IL  | RL |
	|____|______|____|
	|    |      |    |
	| LI |  II  | RI |
	|____|______|____|
	| LR |  IR  | RR |
	|____|______|____|
"""->
function split(A::DenseMatrix, border::Int)
	@assert border >= 1

	N = size(A)
	@assert minimum(N) > 2*border

	Lidx = 1:border
	I1idx = border+1:N[1]-border
	R1idx = N[1]-border+1:N[1]
	I2idx = border+1:N[2]-border
	R2idx = N[2]-border+1:N[2]

	LL = slice(A, Lidx, Lidx)
	LI = slice(A, I1idx, Lidx)
	LR = slice(A, R1idx, Lidx)

	IL = slice(A, Lidx, I2idx)
	II = slice(A, I1idx, I2idx)
	IR = slice(A, R1idx, I2idx)

	RL = slice(A, Lidx, R2idx)
	RI = slice(A, I1idx, R2idx)
	RR = slice(A, R1idx, R2idx)

	SplitMatrix( LL, LI, LR, IL, II, IR, RL, RI, RR )
end

