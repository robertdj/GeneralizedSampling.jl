# ------------------------------------------------------------
# Exact Fourier transforms of Haar wavelet/scaling function

# Fourier transform of Haar scaling function

@doc """
	FourHaarScaling(xi)

The Fourier transform of the Haar scaling function evaluated at `xi`.
"""->
function FourHaarScaling(xi::Number)
	if xi == 0
		return complex(1.0)
	else
		return (1.0 - exp(-2.0*pi*im*xi)) / (2.0*pi*im*xi)
	end
end

@vectorize_1arg Number FourHaarScaling



@doc """
	FourHaarScaling(xi, J, k::Int)

The Fourier transform of the Haar scaling function on scale `J` and translation `k` evaluated at `xi`.
"""->
function FourHaarScaling(xi::Vector, J::Int, k::Int)
	y = exp( -2.0*pi*im*2.0^(-J)*k*xi ) .* 2.0^(-J/2) .* FourHaarScaling(2.0^(-J)*xi)
end

@doc """
	FourHaarScaling(xi, J, k::Vector) -> F

`F[n,m]` is the Fourier transform at translation `k[m]` evaluated at `xi[n]`.
"""->
function FourHaarScaling(xi::Vector, J::Int, k::Vector{Int})
	M = length(xi)
	N = length(k)

	y = Array(Complex{Float64}, M, N)

	for n = 1:N
		@inbounds y[:,n] = FourHaarScaling(xi, J, k[n])
	end

	return y
end


# Fourier transform of Haar wavelet

@doc """
	FourHaarWavelet(xi)

The Fourier transform of the Haar wavelet evaluated at `xi`.
"""->
function FourHaarWavelet(xi::Number)
	if xi == 0
		return complex(0.0)
	else
		return (1.0 - exp(-pi*im*xi))^2.0 / (2.0*pi*im*xi)
	end
end

@vectorize_1arg Number FourHaarWavelet


# Fourier transform of Haar scaling function with translation (J) and dilation (k)
# Entry (n,m) is the Fourier transform at translation k[m] in xi[n]

@doc """
	FourHaarWavelet(xi, J, k)

The Fourier transform of the Haar wavelet on scale `J` and translation `k` evaluated at `xi`.
"""->
function FourHaarWavelet(xi::Vector, J::Int, k::Int)
	y = exp( -2*pi*im*2.0^(-J)*k*xi ) .* 2.0^(-J/2) .* FourHaarWavelet(2.0^(-J)*xi);
end

@doc """
	FourHaarWavelet(xi, J, k::Vector) -> F

`F[n,m]` is the Fourier transform at translation `k[m]` evaluated at `xi[n]`.
"""->
function FourHaarWavelet(xi::Vector, J::Int, k::Vector)
	M = length(xi)
	N = length(k)

	y = Array(Complex{Float64}, M, N)

	for n = 1:N
		@inbounds y[:,n] = FourHaarWavelet(xi, J, k[n])
	end

	return y
end


# ------------------------------------------------------------

# Fourier transform of Daubechies scaling function

@doc """
### FourDaubScaling(xi, N)

The Fourier transform of the Daubechies `N` scaling function evaluated at `xi`.
`N` is the number of zeros at -1.

The function is computed as an 'infinite' product;
to control this there are optional arguments:

- `prec`: Include factors that are numerically smaller than 1-prec.
- `M`: The maximum number of factors.
"""->
function FourDaubScaling{T<:Real}( XI::AbstractArray{T}, N::Int; prec=eps(), M::Int=20 )
	Z = Array(Complex{Float64}, size(XI))
	N_xi = length(XI)

	# Fixed factor in the complex exponential
	Four_idx = -2*pi*im*[0:2*N-1;]

	# Filter coefficients
	C = qmf(wavelet( WT.Daubechies{N}() ))
	C = scale!(C, 1/sum(C))

	# The infinite product for each xi
	for n = 1:N_xi
		xi = XI[n] / 2.0
		Y = y = dot(C, exp(xi*Four_idx))

		m = 1
		while abs(y) <= 1-prec && m <= M
			xi /= 2.0
			y = dot(C, exp(xi*Four_idx))
			Y *= y

			m += 1
		end

		Z[n] = Y
	end

	return Z
end

@doc """
	DaubLowPass(xi, N::Int)

Low-pass filter for Daubechies `N` wavelet.
Uses the Wavelets package to compute the filter coefficients.
"""->
function DaubLowPass(xi::Real, N::Int)
	C = qmf(wavelet( WT.Daubechies{N}() ))
	C = scale!(C, 1/sum(C))

	Y = dot(C, exp(-2*pi*im*xi*[0:2*N-1;]))

	return Y
end

