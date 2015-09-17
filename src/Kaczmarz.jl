@doc """
	REK(A::Matrix, b, x0; prec, maxiter)

General Random Extended Kaczmarz algorithm for solving the least squares problem `min norm(Ax - b)` for complex `A` and `b` with starting point `x0`..
""" ->
function REK{T<:Number}(A::Matrix{T}, b::Vector{T}, x0::Vector{T}; prec=1e-6, maxiter::Int=size(A,2)^2)
	M, N = size(A)
	@assert M == length(b)
	@assert N == length(x0)

	# Create alias tables for sampling row and column indices
	# TODO: Try NumericExtensions package
	A_squared = abs2(A)
	row_norm = vec( sum(A_squared,2) )
	col_norm = vec( sum(A_squared,1) )
	matrix_norm = sum(col_norm)

	row_distribution = Categorical( row_norm/matrix_norm )
	row_sampler = sampler(row_distribution)

	col_distribution = Categorical( col_norm/matrix_norm )
	col_sampler = sampler(col_distribution)

	# The Kaczmarz solver
	z = deepcopy(b)
	x = deepcopy(x0)
	AH = A'

	for iter = 1:maxiter
		# Update z
		col_index = rand(col_sampler)
		col = view(A, :, col_index)
		col_val = -BLAS.dotc(M,col,1,z,1)/col_norm[col_index]
		BLAS.axpy!(col_val, col, z) # z = z + col_val*col

		# Update x
		row_index = rand(row_sampler)
		row = view(AH, :, row_index)
		row_val = (b[row_index] - z[row_index] - BLAS.dotc(N,row,1,x,1))/row_norm[row_index] 
		BLAS.axpy!(row_val, row, x) # x = x + row_val*row

		# Check for convergence
		if mod(iter,8*N) == 0
			xnorm = norm(x)
			Ax = A*x
			Atz = A'*z
			test1 = norm(Ax - b + z)/(matrix_norm*xnorm)
			test2 = norm(Atz)/(matrix_norm^2*xnorm)

			if test1 < prec && test2 < prec
				println("Number of iterations: ", iter)
				break
			end
		end
	end

	return x
end


@doc """
	REK(T::Freq2wave1D, b, x0; prec, maxiter)

Random Extended Kaczmarz algorithm for solving the least squares problem `min norm(Tx - b)` with starting point `x0`.

By default, `prec=1e-6` and `maxiter=length(x0)`.
""" ->
function REK(T::Freq2wave1D, b::Vector{Complex{Float64}}, x0::Vector{Complex{Float64}}; prec::Float64=1e-4, maxiter::Int=length(x0)^2)
	M, N = size(T)
	# TODO: Scale b with weights??

	# Wavelet scale
	J = Int(log2(N))

	# Squared norms of rows, columns and matrix
	row_norm = N*abs2(T.column1)
	matrix_norm = sum(row_norm)
	col_norm = matrix_norm / N

	# Create alias tables for sampling row and column indices
	row_dist = Categorical( row_norm/matrix_norm )
	row_sampler = sampler(row_dist)

	col_sampler = [1:N;]

	# The Kaczmarz solver
	z = deepcopy(b)
	x = deepcopy(x0)

	column = cis(-2.0*pi*2.0^(-J)*T.samples)
	col = similar(column)
	row_base = 2.0*pi*2.0^(-J)*[0:N-1;]
	row = Array(Complex{Float64}, N)

	for iter = 1:maxiter
		# Update z
		col_index = rand(col_sampler)-1
		for m = 1:M
			@inbounds col[m] = T.column1[m]*column[m]^col_index
		end
		col_val = -BLAS.dotc(M,col,1,z,1)/col_norm
		BLAS.axpy!(col_val, col, z) # z = z + col_val*col
		#z += col_val*col

		# Update x
		row_index = rand(row_sampler)
		for n = 1:N
			row[n] = cis( row_base[n]*T.samples[row_index] )*conj(T.column1[row_index])
		end
		row_val = (b[row_index] - z[row_index] - BLAS.dotc(N,row,1,x,1))/row_norm[row_index]
		BLAS.axpy!(row_val, row, x) # x = x + row_val*row
		#x += row_val*row

		# Check for convergence
		if mod(iter,8*N) == 0
			xnorm = norm(x)
			Ax = T*x
			Atz = H(T,z)
			test1 = norm(Ax - b + z)/(matrix_norm*xnorm)
			test2 = norm(Atz)/(matrix_norm^2*xnorm)

			if test1 < prec && test2 < prec
				println("Number of iterations: ", iter)
				break
			end
		end
	end

	return x
end

