# ------------------------------------------------------------
# Exact Fourier transforms of Haar wavelet/scaling function

@doc """
	FourHaarScaling(xi[, J, k])

The Fourier transform of the Haar scaling function on scale `J` and translation `k` evaluated at `xi`.
If not supplied, `J = 0` and `k = 0`.
"""->
function FourHaarScaling{T<:Real}(xi::T)
	if xi == 0
		return one(Complex{Float64})
	else
		return (1.0 - exp(-2.0*pi*im*xi)) / (2.0*pi*im*xi)
	end
end

@doc """
	FourHaarWavelet(xi[, J, k])

The Fourier transform of the Haar wavelet on scale `J` and translation `k` evaluated at `xi`.
If not supplied, `J = 0` and `k = 0`.
"""->
function FourHaarWavelet{T<:Real}(xi::T)
	if xi == 0
		return zero(Complex{Float64})
	else
		return (1.0 - exp(-pi*im*xi))^2 / (2.0*pi*im*xi)
	end
end


# Vectorize and add dilation/scale and translation. 
# Identical for scaling and wavelet
for name in [:FourHaarScaling, :FourHaarWavelet]
	@eval begin
		@vectorize_1arg Real $name

		function $name{T<:Real}(xi::DenseArray{T}, J::Int)
			M = length(xi)
			C2 = 2.0^(-J)
			C = 2.0^(-J/2)
			y = Array(Complex{Float64}, M)
			for m = 1:M
				@inbounds y[m] = C*$name( C2*xi[m] )
			end
			return y
		end

		function $name{T<:Real}(xi::DenseArray{T}, J::Int, k::Int)
			y = $name(xi, J)
			D = exp( -2.0*pi*im*2.0^(-J)*k*xi )
			had!(y, D)
			return y
		end
	end
end


# ------------------------------------------------------------
# Fourier transform of Daubechies wavelets

@doc """
	FourDaubScaling(xi, N; ...)

The Fourier transform of the Daubechies `N` scaling function evaluated at `xi`.
`N` is the number of zeros at -1.

The function is computed as an 'infinite' product;
to control this there are optional arguments:

- `prec`: Include factors that are numerically smaller than 1-prec.
- `maxCount`: The maximum number of factors.
"""->
function FourDaubScaling{T<:Real}( XI::DenseArray{T}, N::Int; prec=eps(), maxCount=100)
	Z = Array(Complex{Float64}, size(XI))

	# Filter coefficients
	C = wavefilter( string("db", N) )
	scale!(C, 1/sum(C))

	M = length(XI)
	almost1 = 1.0 - prec

	# Fixed factor in the complex exponential
	Four_idx = 2.0*pi*im*[0:2*N-1;]
	Y = y = zero(Complex{Float64})

	# The infinite product of low pass filters for each xi
	for m = 1:M
		xi = XI[m] / 2.0
		Y = y = dot(C, exp(xi*Four_idx))

		# When xi is (almost) an even integer y is (approx) 1.
		# To ensure that the while loop below does not exit 
		# prematurely, a minimum number of iterations is set,
		# which is the number of iterations needed for abs(xi) < 1
		minCount = (abs(xi) > 1.0 ? ceil(Int, log2(abs(xi))) : 1)
		for count = 1:minCount
			xi /= 2.0
			Y *= dot(C, exp(xi*Four_idx))
		end

		count = 1
		y_almost1 = false
		while !y_almost1 && count <= maxCount
			xi /= 2.0
			y = dot(C, exp(xi*Four_idx))
			Y *= y

			y_almost1 = abs(y) <= almost1
			count += 1
		end

		Z[m] = Y
	end

	return Z
end


@doc """
### FourDaubWavelet(xi, N[, J, k]; ...)

The Fourier transform of the Daubechies `N` wavelet function evaluated at `xi`.
`N` is the number of zeros at -1.

The optional arguments are passed to `FourDaubScaling`.
"""->
function FourDaubWavelet{T<:Real}( XI::DenseArray{T}, N::Int; args... )
	# Fixed factor in the complex exponential
	Four_idx = 2*pi*im*[0:2*N-1;]

	# Filter coefficients
	C = wavefilter( string("db", N) )
	C = scale!(C, 1/sum(C))

	# Fourier transform of low pass filter
	XI /= 2
	Z = FourDaubScaling(XI, N; args...)

	M = length(XI)

	# The infinite product for each xi
	# TODO: Vectorize this?
	for m = 1:M
		xi = XI[m]
		# High pass filter
		H = conj(dot( C, exp((xi+0.5)*Four_idx) ))

		Z[m] *= exp( 2*pi*im*xi )*H
	end

	return Z
end

# Vectorize and add dilation/scale and translation. 
# Identical for scaling and wavelet
for name in [:FourDaubScaling, :FourDaubWavelet]
	@eval begin
		function $name{T<:Real}(xi::DenseArray{T}, N::Int, J::Int; args...)
			scale_xi = scale( 2.0^(-J), xi )
			y = $name( scale_xi, N; args... )
			scale!( 2.0^(-J/2), y )
			return y
		end

		function $name{T<:Real}(xi::DenseArray{T}, N::Int, J::Int, k::Int; args...)
			y = $name(xi, N, J; args... )
			D = exp( -2.0*pi*im*2.0^(-J)*k*xi )
			had!(y, D)
			return y
		end
	end
end

