@doc """
	Kaczmarz(A, b; maxiter)

General Kaczmarz algorithm for solving `Ax = b` for complex `A` and `b`.
By default, `maxiter=100`.
""" ->
function Kaczmarz(A::Matrix, b::Vector; maxiter=100)
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

	for iter = 1:maxiter
		col_index = rand(col_sampler)
		col = A[:,col_index]
		z -= dot(col, z)/col_norm[col_index] * col

		row_index = rand(row_sampler)
		row = vec( A[row_index, :] )
		x += (b[row_index] - z[row_index] - BLAS.dotu(N,x,1,row,1))/row_norm[row_index] * conj(row)
	end

	return x
end

