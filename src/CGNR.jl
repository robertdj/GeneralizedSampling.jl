@doc """
	cgnr(A::Matrix, b::Vector, x0::Vector; prec, maxiter)

Conjugate gradient for normal equations residual method for solving the least squares problem `min norm( Ax - b )`.

- `A` is the coefficient matrix (`M-by-N`)
- `b` is the observation vector (`M`).
- `x0` is the starting vector (`N`).

The iteration stops when `norm(xnew - xold) < prec` or after at most `maxiter` iterations.
"""->
function cgnr{T<:Number}(A::Matrix{T}, b::Vector{T}, x0::Vector{T}; prec=sqrt(eps()), maxiter=length(x0))
	@assert size(A,1) == length(b)
	@assert size(A,2) == length(x0)

	# Initialize
	x = copy(x0)
	r = b - A*x
	z = A'*r
	p = copy(z)
	tone = one(T)

	for iter = 1:maxiter
		ztz = norm(z)^2
		mu = ztz / norm(A*p)^2
		BLAS.axpy!(mu, p, x) # x = x + mu*p
		BLAS.gemv!('N', convert(T,-mu), A, p, tone, r) # r = r - mu*A*p
		z = A'*r
		tau = norm(z)^2 / ztz
		# p = z + tau*p
		scale!(p, tau)
		BLAS.axpy!(tone, z, p)

		# Check for convergence: |xnew - xold|
		if mu*norm(p) < prec
			println("Number of iterations: ", iter)
			break
		end
	end

	return x
end

@doc """
	cgnr(T::Freq2wave, b, x0; ...)

Conjugate gradient for normal equations residual method for `Freq2wave`.
The initial point `x0` must be of the same dimension as `T`.
"""->
function cgnr{D}(T::Freq2wave{D}, b::Vector{Complex{Float64}}, x0::DenseArray{Complex{Float64},D}; prec=sqrt(eps()), maxiter=length(x0))
	@assert size(T,1) == length(b)
	@assert wsize(T) == size(x0)

	# Initialize
	x = copy(x0)

	y = similar(b)
	mul!(T, x, y) # y = T*x
	r = b - y

	z = similar(x0)
	mulT!(T, r, z) # z = T'*r
	p = copy(z)

	cone = one(Complex{Float64})

	for iter = 1:maxiter
		ztz = vecnorm(z)^2
		mul!(T, p, y) # y = T*x
		mu = ztz / vecnorm(y)^2
		BLAS.axpy!(mu, p, x) # x = x + mu*p
		BLAS.axpy!(-mu, y, r) # r = r - mu*y
		mulT!(T, r, z) # z = T'*r
		tau = vecnorm(z)^2 / ztz
		# p = z + tau*p
		scale!(p, tau)
		BLAS.axpy!(cone, z, p)

		# Check for convergence: |xnew - xold|
		if mu*vecnorm(p) < prec
			println("Number of iterations: ", iter)
			break
		end
	end

	return x
end

