#=
@doc """
Conjugate gradient normal equation residual method.
"""->
=#
function cgnr{T<:Number}(A::Matrix{T}, b::Vector{T}, x::Vector{T}; prec=eps(), maxiter=100)
	# Initialize
	r = b - A*x
	At = A'
	z = At*r
	p = deepcopy(z)

	for iter = 1:maxiter
		ztz = norm(z)^2
		mu = ztz / norm(A*p)^2
		xdiff = mu*p
		x += xdiff
		r -= mu*A*p
		z = At*r
		tau = norm(z)^2 / ztz
		p = z + tau*p

		# Check for convergence
		if norm(xdiff) < prec
			println("Number of iterations: ", iter)
			break
		end
	end

	return x
end

