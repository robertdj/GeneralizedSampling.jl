@doc """
	had!(A, B)

In-place Hadamard product/element-wise matrix multiplication; the first argument is modified.
"""->
function had!{T<:Number}(A::Matrix{T}, B::Matrix{T})
	m,n = size(A)
	@assert (m,n) == size(B)
	for j in 1:n
		for i in 1:m
			@inbounds A[i,j] *= B[i,j]
		end
	end
end

function had!{T<:Number}(A::Vector{T}, B::Vector{T})
	m = length(A)
	@assert m == length(B)
	for i in 1:m
		@inbounds A[i] *= B[i]
	end
end


@doc """
	isuniform(x; prec)

Test if the sampling points in `x` are on a uniform grid with precision `prec`.

`x` is a vector for 1D points and an `M`-by-2 matrix for 2D points.
""" ->
function isuniform( x::Vector; prec::Float64=sqrt(eps()) )
	M = length(x)

	if M <= 2
		return true
	end

	diff = abs(x[1] - x[2])

	for n = 3:M
		d = abs(x[n-1] - x[n])
		if abs(d - diff) > prec
			return false
		end
	end

	return true
end

function isuniform( points::Matrix; prec::Float64=sqrt(eps()) )
	M, D = size(points)
	@assert D == 2

	if M <= 2
		return true
	end

	x = sort( points[:,1] )
	uniquex = unique(x)
	Mx = length(uniquex)

	y = sort( points[:,2] )
	uniquey = unique(y)
	My = length(uniquey)

	if !isuniform(uniquex; prec=prec) || !isuniform(uniquey; prec=prec) || Mx*My != M
		return false
	end

	return true
end


@doc """
	grid(Mx, My, scale)

2D points on an `Mx`-by-`My` grid centered around the origin.
With even `M`'s the grid has one extra point on the negative values.

By default, `My` = `Mx`.
The points are scaled by `scale` which by default is 1.
"""->
function grid(Mx::Int, My::Int=Mx, scale::Float64=1.0)
	startx = -div(Mx,2)
	endx = (isodd(Mx) ? -startx : -startx-1)
	# The points are sorted by the x coordinate
	x = kron([startx:endx;], ones(My))

	starty = -div(My,2)
	endy = (isodd(My) ? -starty : -starty-1)
	y = repmat([starty:endy;], Mx)

	points = scale*[x y]
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
	@assert size(xi,2) == 2

	# Compute the area of the Voronoi cells corresponding to xi
	# This is done using R and its deldir package
	# TODO: Find a Julia way
	reval("library(deldir)")
	globalEnv[:x] = RObject( xi[:,1] )
	globalEnv[:y] = RObject( xi[:,2] )
	globalEnv[:K] = RObject( bandwidth )
	reval("V = deldir(x, y, rw=c(-K,K,-K,K))")

	area = rcopy("V\$summary\$dir.area")
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

function density(xi::Matrix{Float64}, bandwidth::Number)
	M, dim = size(xi)
	@assert dim == 2

	# Compute Voronoi tesselation
	reval("library(deldir)")
	globalEnv[:x] = RObject( xi[:,1] )
	globalEnv[:y] = RObject( xi[:,2] )
	globalEnv[:K] = RObject( bandwidth )
	reval("V = deldir(x, y, rw=c(-K,K,-K,K))")

	# Corners of Voronoi cells
	x1 = rcopy("V\$dirsgs\$x1")
	y1 = rcopy("V\$dirsgs\$y1")
	x2 = rcopy("V\$dirsgs\$x2")
	y2 = rcopy("V\$dirsgs\$y2")

	# Edge-sampling point relation
	ind1 = rcopy("V\$dirsgs\$ind1")
	ind = round(Int64, ind1)

	# Compute the distance from each xi to the corners of its Voronoi cell
	density = 0.0
	Ncorner = length(ind)
	for n = 1:Ncorner
		idx = ind[n]
		# Distance^2 from corners to one of the xi's that has this edge
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

The fractional part of x as a number in [-0.5, 0.5).
"""->
function frac(x::Array{Float64})
	y = copy(x)
	frac!(y)

	return y
end

function frac!(x::Array{Float64})
	N = length(x)
	for n = 1:N
		@inbounds x[n] -= round(x[n])
	end
end


@doc """
	wavename(name)

Return the characteristics of a wavelet needed for computations.
"""->
function wavename(name::AbstractString)
	lname = lowercase(name)

	if lname == "haar"
		return "Haar"
	elseif lname[1:2] == "db" && typeof(parse(lname[3:end])) <: Integer
		return ("Daubechies", parse(lname[3:end]))
	else
		error("Uknown wavelet name")
	end
end

@doc """
	wavefilter(name)

Return the low pass filter coefficients of the wavelet `name`.

Uses the `Wavelets` package.
"""->
function wavefilter(name::AbstractString)
	parsed_name = wavename(name)

	if parsed_name == "Haar"
		C = wavelet( WT.Haar() )
	elseif parsed_name[1] == "Daubechies"
		C = wavelet( WT.Daubechies{parsed_name[2]}() )
	else
		error("Uknown wavelet name")
	end

	return C.qmf
end

