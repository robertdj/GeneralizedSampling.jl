@doc """
	cgnr(A::Matrix, b::Vector, x0::Vector; prec, maxiter)

Conjugate gradient for normal equations residual method for solving the least squares problem `min norm( Ax - b )`.

- `A` is the coefficient matrix (`M-by-N`)
- `b` is the observation vector (`M`).
- `x0` is the starting vector (`N`).

The iteration stops when `norm(xnew - xold) < prec` or after at most `maxiter` iterations.
"""->
function cgnr{T<:Number}(A::AbstractMatrix{T}, b::AbstractVector{T}, x0::AbstractVector{T}; prec=LARGE_EPS, maxiter=length(x0))
	@assert size(A,1) == length(b)
	@assert size(A,2) == length(x0)
	@assert prec >= SMALL_EPS
	@assert maxiter >= 1

	# Initialize
	x = copy(x0)
	y = A*x
	r = b - y
	z = A'*r
	p = copy(z)
	oneT = one(T)

	for iter = 1:maxiter
		ztz = norm(z)^2
		A_mul_B!(y, A, p) # y = A*p
		mu = ztz / norm(y)^2
		BLAS.axpy!(mu, p, x) # x = x + mu*p
		xdiff = mu*norm(p)

		BLAS.axpy!(-mu, y, r) # r = r - mu*y
		Ac_mul_B!(z, A, r) # z = A'*r
		tau = norm(z)^2 / ztz

		# p = z + tau*p
		scale!(p, tau)
		BLAS.axpy!(oneT, z, p)

		# Check for convergence: |xnew - xold|
		if xdiff < prec
			break
		end
	end

	return x
end

function cgnr{D}(T::Freq2Wave{D}, B::AbstractArray{Complex{Float64},D}, x0::AbstractArray{Complex{Float64},D}; args...)
	# TODO: flatten_view?
	b = reshape_view(B, size(T,1))
	cgnr(T, b, x0; args...)
end

@doc """
	cgnr(T::Freq2wave, b, x0; ...) -> x

Conjugate gradient for normal equations residual method for `Freq2Wave`.
The initial point `x0` must be of the same dimension as `T`.
"""->
function cgnr{D}(T::Freq2Wave{D}, b::AbstractVector{Complex{Float64}}, x0::AbstractArray{Complex{Float64},D}; prec=LARGE_EPS, maxiter=max(length(x0),50))
	# TODO: Assertions in macro?
	@assert size(T,1) == length(b)
	@assert wsize(T) == size(x0)
	@assert prec >= SMALL_EPS
	@assert maxiter >= 1

	# Initialize
	x = copy(x0)
	y = T*x
	r = b - y
	z = T'*r
	p = copy(z)
	Cone = one(Complex{Float64})

	for iter = 1:maxiter
		ztz = vecnorm(z)^2
		A_mul_B!(y, T, p) # y = T*p
		mu = ztz / vecnorm(y)^2
		BLAS.axpy!(mu, p, x) # x = x + mu*p
		xdiff = mu*norm(p)

		BLAS.axpy!(-mu, y, r) # r = r - mu*y
		Ac_mul_B!(z, T, r) # z = T'*r
		tau = vecnorm(z)^2 / ztz

		# p = z + tau*p
		scale!(p, tau)
		BLAS.axpy!(Cone, z, p)

		# Check for convergence: |xnew - xold|
		if xdiff < prec
			break
		end
	end

	return x
end

