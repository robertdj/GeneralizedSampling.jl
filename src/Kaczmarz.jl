# General Kaczmarz algorithm for solving Ax = b
#
# Input:
# A : Coefficient matrix
# b : Vector
# x0 : Initial estimate
#
# K : Max number of iterations

#@doc """
#""" ->
@debug function Kaczmarz(A::Matrix, b::Vector; prec::Float64=1e-2)
	M, N = size(A)

	# Create alias tables for sampling row and column indices
	A_squared = abs(A).^2
	row_norm = vec( sum(A_squared,2) )
	col_norm = vec( sum(A_squared,1) )
	matrix_norm = sum(col_norm)

	row_distribution = Categorical( row_norm/matrix_norm )
	row_sampler = sampler(row_distribution)

	col_distribution = Categorical( col_norm/matrix_norm )
	col_sampler = sampler(col_distribution)

	# The Kaczmarz solver
	x = fill(0.0 + 0.0*im, N)
	xold = zeros(Float64, N)
	z = deepcopy(b)
	zold = deepcopy(b)

	for iter = 1:N^2
		col_index = rand(col_sampler)
		col = A[:,col_index]
		z -= dot(col, z)/col_norm[col_index] * col

		row_index = rand(row_sampler)
		row = vec( A[row_index, :] )
		x += (b[row_index] - z[row_index] - mydot(x,row))/row_norm[row_index] * conj(row)

		@printf("%u : %f\n", iter, norm(x))
	end

	return x
end

#function mydot{T<:Number}(a::Vector{T}, b::Vector{T})
function mydot(a::Vector, b::Vector)
	N = length(a)

	d = 0.0
	for n = 1:N
		d += a[n]*b[n]
	end

	return d
end

