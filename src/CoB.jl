# ------------------------------------------------------------
# 1D functions

# Change of basis matrix from frequencies to Haar wavelets 

@doc """
	freq2Haar(xi::Vec, J::Int) -> T

Change of basis matrix from frequency responses at `xi` to coefficients of
Haar scaling functions at scale `J`.
"""->
function freq2Haar( xi::Vector, J::Int )
	k = [0:2^J-1;]

	T = FourHaarScaling( xi, J, k )
end

@doc """
	freq2Haar(xi::Vec, J::Vec{Int}) -> T

Change of basis matrix from frequency responses at `xi` to coefficients of
Haar scaling functions at scale `J[1]` and Haar wavelet functions at scales 
`J[2:end]`.
"""->
function freq2Haar( xi::Vector, J::Vector{Int} )
	M = length(xi)
	J2 = 2.^J
	N = cumsum( J2 )

	T = Array(Complex{Float64}, M, N[end])

	# Scaling functions
	k = [0:N[1]-1;]
	@inbounds T[:, 1:N[1]] = FourHaarScaling( xi, J[1], k )

	# Wavelet functions
	for n = 2:length(J)
		j = J[n]
		k = [0:J2[n]-1;]
		@inbounds T[:, N[n-1]+1:N[n]] = FourHaarWavelet( xi, j, k )
	end

	return T
end


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
# (Row) Weigthed matrices

@doc """ 
	freq2Haar(xi, J, weights) -> T

The rows of `freq2Haar(xi, J)` weighted by the vector `weights`.
"""->
function freq2Haar( T::Array{Complex{Float64},2}, weights::AbstractVector )
	@assert size(T,1) == length(weights)

	broadcast!(*, T, weights, T)

	return T
end

function freq2Haar( xi::Vector, J::Int, weights::AbstractVector )
	@assert length(xi) == length(weights)

	T = freq2Haar( xi, J )
	freq2Haar( T, weights )

	return T
end


function freq2Haar( xi::Vector, J::Vector{Int}, weights::AbstractVector )
	@assert length(xi) == length(weights)

	T = freq2Haar( xi, J )
	freq2Haar( T, weights )

	return T
end


# ------------------------------------------------------------
# 2D functions

@doc """
	freq2Haar(xi1::Vector, xi2::Vector, J) -> T

Change of basis matrix from 2D frequency responses in rectangular area `xi1`-by-`xi2`.
"""->
function freq2Haar( xi_x::Vector, xi_y::Vector, J::Int )
	Tx = freq2Haar(xi_x, J)
	Ty = freq2Haar(xi_y, J)

	return Tx, Ty
end


# General points in 2D

@doc """
	freq2Haar(xi::Matrix, J) -> T

Change of basis matrix from 2D frequency responses at the rows of `xi`.
"""->
function freq2Haar( xi::AbstractMatrix, J::Int )
	M = size(xi,1)
	T = Array(Complex{Float64}, M, 4^J)

	k = [0:2^J-1;]
	Tx = FourHaarScaling( xi[:,1], J, k)
	Ty = FourHaarScaling( xi[:,2], J, k)

	for m = 1:M
		@inbounds T[m,:] = kron( Tx[m,:], Ty[m,:] )
	end

	return T
end

# With row weighing

@doc """ 
	freq2Haar(xi::Matrix, J, weights) -> T

The rows of `freq2Haar(xi, J)` weighted by the vector `weights`.
"""->
function freq2Haar( xi::AbstractMatrix, J::Int, weights::AbstractVector )
	@assert size(xi,1) == length(weights)

	T = freq2Haar( xi, J )
	freq2Haar( T, weights )

	return T
end


