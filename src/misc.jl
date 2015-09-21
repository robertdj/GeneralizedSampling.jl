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
	upscale(x)

When `x` is a vector, each element is repeated twice.
""" ->
function upscale(x::Vector)
	N = length(x)
	kron(x, ones(2))
end


@doc """
	isuniform(x; prec)

Test if the sampling points in `x` are on a uniform grid with precision `prec`.
""" ->
function isuniform(x::Vector; prec::Float64=eps())
	N = length(x)

	if N == 1
		error("The vector must have at least two elements")
	elseif N == 2
		return true
	end

	diff = abs(x[1] - x[2])

	for n = 3:N
		d = abs(x[n-1] - x[n])
		if abs(d - diff) > prec
			return false
		end
	end

	return true
end

@doc """
	weights(xi, bandwidth)

Compute weights for sampling points `xi`.
When `xi` is a vector it is assumed to be *sorted*.
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

# TODO: Move package call to module
using RCall
function weights(xi::Matrix, bandwidth::Real)
	@assert size(xi,2) == 2

	# Compute the area of the Voronoi cells corresponding to xi
	# This is done using R and its deldir package
	# TODO: Find a Julia way
	reval("library(deldir)")
	globalEnv[:x] = sexp( xi[:,1] )
	globalEnv[:y] = sexp( xi[:,2] )
	globalEnv[:K] = sexp( bandwidth )
	reval("V = deldir(x, y, rw=c(-K,K,-K,K))")

	area = rcopy(reval("V\$summary\$dir.area"))
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

function density(xi::Matrix, bandwidth::Number)
	M, dim = size(xi)
	@assert dim == 2


	error("Not implemented yet")
end

@doc """
	frac(x)
	frac!(x)

The fractional part of x as a number in [-0.5, 0.5).
"""->
function frac(x::Array{Float64})
	y = deepcopy(x)
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
function wavename(name::String)
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
function wavefilter(name::String)
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

