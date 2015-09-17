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
	x = copy(x0)
	r = b - A*x
	z = A'*r
	p = copy(z)

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


function cgnr(T::Freq2wave1D, b::Vector{Complex{Float64}}, x0::Vector{Complex{Float64}}; prec=sqrt(eps()), maxiter=length(x))
	# Initialize
	x = copy(x0)

	y = similar(b)
	mul!(T, x, y) # y = T*x
	r = b - y

	z = similar(x0)
	mulT!(T, r, z) # z = T'*r
	p = copy(z)

	for iter = 1:maxiter
		ztz = norm(z)^2
		mul!(T, p, y) # y = T*x
		mu = ztz / norm(y)^2
		x = x + mu*p
		#BLAS.axpy!(mu, p, x) # x = x + mu*p
		r = r - mu*y
		mulT!(T, r, z) # z = T'*r
		tau = norm(z)^2 / ztz
		p = z + tau*p
		#scale!(p, tau)
		#BLAS.axpy!(1.0, z, p)

		# Check for convergence: |xnew - xold|
		#=
		if mu*norm(p) < prec
			println("Number of iterations: ", iter)
			break
		end
		=#
	end

	return x
end

