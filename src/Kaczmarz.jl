#@doc """
#	Kaczmarz(A, b; prec)
#
#General Kaczmarz algorithm for solving `Ax = b` for complex `A` and `b`.
#By default, `prec=1e-4`.
#""" ->
function Kaczmarz(A::Matrix, b::Vector; prec=1e-4)
	M, N = size(A)
	maxiter = N^2

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
		# Update z
		col_index = rand(col_sampler)
		col = A[:,col_index]
		zdiff = dot(col, z)/col_norm[col_index] * col
		z -= zdiff

		# Update x
		row_index = rand(row_sampler)
		row = vec( A[row_index, :] )
		xdiff = (b[row_index] - z[row_index] - BLAS.dotu(N,x,1,row,1))/row_norm[row_index] * conj(row)
		x += xdiff

		# Check for convergence
		if mod(iter,8*N) == 0
			xnorm = norm(x)
			test1 = norm(A*x - b + z)/(matrix_norm*xnorm)
			test2 = norm(A'*z)/(matrix_norm^2*xnorm)

			if test1 < prec && test2 < prec
				println("Number of iterations: ", iter)
				break
			end
		end
	end

	return x
end

function Kaczmarz(xi::Vector, b::Vector, x0::Vector; prec::Float64=1e-2, maxiter=100)
	M = length(xi)
	N = length(x0)

	# Wavelet scale
	J = Int(log2(N))

	# Squared norms of rows, columns and matrix
	# TODO: This is only for Haar
	wavelet_fourier = FourHaarScaling(xi, J, 0)
	row_norm = N*abs(wavelet_fourier).^2
	matrix_norm = sum(row_norm)
	col_norm = matrix_norm / N

	# Create alias tables for sampling row and column indices
	row_dist = Categorical( row_norm/matrix_norm )
	row_sampler = sampler(row_dist)

	col_sampler = [1:N;]

	# The Kaczmarz solver
	x = complex(deepcopy(x0))
	z = complex(deepcopy(b))

	for iter = 1:maxiter
		# Update z
		col_index = rand(col_sampler)
		col = exp(-2.0*pi*im*2.0^(-J)*(col_index-1)*xi)
		broadcast!(*, col, col, wavelet_fourier)
		#had!(col, wavelet_fourier)
		z -= dot(col, z)/col_norm * col

		# Update x
		row_index = rand(row_sampler)
		row = exp(-2.0*pi*im*2.0^(-J)*[0:N-1;]*xi[row_index])
		scale!(row, wavelet_fourier[row_index])
		x += (b[row_index] - z[row_index] - BLAS.dotu(N,x,1,row,1))/row_norm[row_index] * conj(row)
	end

	return x
end


function Kaczmarz(xi::Vector, weights::Vector, Jstart::Integer, Jend::Integer; prec::Float64=1e-3, maxiter=100)
end

