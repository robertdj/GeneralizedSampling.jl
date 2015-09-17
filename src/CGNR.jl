@doc """
	cgnr(A::Matrix, b::Vector, x0::Vector; prec, maxiter)

Conjugate gradient normal equation residual method for solving the least squares problem `min norm( Ax - b )`.

- `A` is the coefficient matrix (`M-by-N`)
- `b` is the observation vector (`M`).
- `x0` is the starting vector (`N`).

The iteration stops when `norm(xnew - xold) < prec` or after at most `maxiter` iterations.
"""->
function cgnr{T<:Number}(A::Matrix{T}, b::Vector{T}, x0::Vector{T}; prec=sqrt(eps()), maxiter=length(x))
	# Initialize
	x = deepcopy(x0)
	r = b - A*x
	z = A'*r
	p = deepcopy(z)

	for iter = 1:maxiter
		ztz = norm(z)^2
		mu = ztz / norm(A*p)^2
		BLAS.axpy!(mu, p, x) # x = x + mu*p
		BLAS.gemv!('N', convert(T,-mu), A, p, one(T), r) # r = r - mu*A*p
		z = A'*r
		tau = norm(z)^2 / ztz
		# p = z + tau*p
		scale!(p, tau)
		BLAS.axpy!(1.0, z, p)

		# Check for convergence: |xnew - xold|
		if mu*norm(p) < prec
			println("Number of iterations: ", iter)
			break
		end
	end

	return x
end


