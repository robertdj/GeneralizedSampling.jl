@doc """
In-place Hadamard product/element-wise matrix multiplication.
"""->
function had!{T<:Number}(A::Matrix{T},B::Matrix{T})
	m,n = size(A)
	@assert (m,n) == size(B)
	for j in 1:n
		for i in 1:m
			@inbounds A[i,j] *= B[i,j]
		end
	end
	return A
end

function had!{T<:Number}(A::Vector{T}, B::Vector{T})
	m = length(A)
	@assert m == length(B)
	for i in 1:m
		@inbounds A[i] *= B[i]
	end

	return A
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
	weights(xi, bandwidth)

Compute weights for sampling points `xi`.
When `xi` is a vector it is assumed to be *sorted*.
"""->
function weights(xi::Vector, bandwidth::Real)
	@assert issorted(xi)

	# TODO: Check that bandwidth is sufficiently large

	N = length(xi)
	mu = Array(Float64, N)

	# Boundary cases
	lower_boundary = xi[N] - 2*bandwidth
	mu[1] = 0.5*(xi[2] - lower_boundary)

	upper_boundary = xi[1] + 2*bandwidth
	mu[N] = 0.5*(upper_boundary - xi[N])

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
	density(xi, bandwidth)

Compute the density of the sampling points `xi` with bandwidth `bandwidth`:

- In 1D, `xi` must be a vector of *sorted* points and the density is the largest difference between two consecutive points in `xi`.
- In 2D, `xi` must be an `M-by-2` matrix and the density is the smallest number δ that allows a covering of the "bandwidth area" with circles centered in the points of `xi` and radii δ.
Currently, the "bandwidth area" is a square centered at the origin and with sidelength `2*bandwidth`. 
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
		diff = xi[n] - xi[n-1]
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
	x - round(x)
end

function frac!(x::Array{Float64})
	N = length(x)
	for n = 1:N
		@inbounds x[n] -= round(x[n])
	end
end

#=
@doc """
Compute `A*b` and `A'*c` using the NFFT where A is a Fourier to wavelet change of basis matrix.

- `xi` is the (non-uniform) nodes.
- `b` is the vector to be multiplied.
"""->
=#
function mul{T<:Real}(xi::Vector{T}, b::Vector{Complex{T}}, c::Vector{Complex{T}})
	@assert length(xi) == length(c)

	N = length(b)
	# TODO: Write this nicer: FourHaarScaling needs the scaled xi and NFFT needs the frac part of the scaled xi
	xx = scale(xi, 1/N)
	x = frac(xx)
	p = NFFTPlan(x, N)
	y = nfft(p, b)

	J = Int(log2(N))
	D = 2.0^(-J/2).*FourHaarScaling(xx).*exp(-pi*im*x*N)
	had!(y, D)

	z = conj(D).*c
	z = nfft_adjoint(p, z)

	return y, z
end

function mul2(xi, J, k, b::Vector, c::Vector)
	A = FourHaarScaling(xi, J, k)
	#=
	M = length(xi)
	N = length(k)

	A = Array(Complex{Float64}, M, N)
	for n = 1:N
		@inbounds A[:,n] = exp( -2.0*pi*im*2.0^(-J)*k[n]*xi )
	end
	=#

	return A*b, A'*c
end

