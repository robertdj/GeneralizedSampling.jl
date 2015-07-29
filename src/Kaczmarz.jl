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
	z = deepcopy(b)

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

function Kaczmarz(xi::Vector, b::Vector, x0::Vector, J::Integer; prec::Float64=1e-2, maxiter=100)
	# TODO: Test input
	M = length(xi)
	N = length(x0)

	# Squared norms of rows, columns and matrix
	wavelet_fourier = FourHaarScaling(xi, J, 0)
	row_norm = N*abs(wavelet_fourier).^2
	matrix_norm = sum(row_norm)
	col_norm = matrix_norm / N

	# Create alias tables for sampling row and column indices
	row_dist = Categorical( row_norm/matrix_norm )
	row_sampler = sampler(row_dist)

	#col_sampler = [1:N;]
	col_sampler = sampler(Categorical( ones(Float64,N)/N ))

	# The Kaczmarz solver
	x = complex(deepcopy(x0))
	z = complex(deepcopy(b))

	for iter = 1:maxiter
		# Update z
		col_index = rand(col_sampler)
		col = exp(-2.0*pi*im*2.0^(-J)*(col_index-1)*xi)
		broadcast!(*, col, col, wavelet_fourier)
		z -= dot(col, z)/col_norm * col

		# Update x
		row_index = rand(row_sampler)
		row = exp(-2.0*pi*im*2.0^(-J)*[0:N-1;]*xi[row_index])
		scale!(row, wavelet_fourier[row_index])
		x += (b[row_index] - z[row_index] - BLAS.dotu(N,x,1,row,1))/row_norm[row_index] * conj(row)
	end

	return x
end
