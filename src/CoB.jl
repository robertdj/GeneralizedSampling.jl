# ------------------------------------------------------------
# 1D functions

# Change of basis matrix from frequencies to Haar wavelets 

# Input:
# xi: The frequencies 
# J: Scale of the wavelet transform
# 
# Output:
# T: Change of basis matrix from frequencies to dilated Haar scaling functions

function freq2Haar( xi::Vector, J::Int )
	k = [0:2^J-1;]

	T = FourHaarScaling( xi, J, k )
end

# Input:
# xi: The frequencies 
# J: The scales of the wavelet transform
# 
# Output:
# T: Change of basis matrix from frequencies to Haar scaling and wavelet functions
# Scaling functions are at scale J[1] and wavelet functions are at scale J[2:end]

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

function FourHaarScaling(xi::Number)
	if xi == 0
		return complex(1.0)
	else
		return (1.0 - exp(-2.0*pi*im*xi)) / (2.0*pi*im*xi)
	end
end

@vectorize_1arg Number FourHaarScaling


# Fourier transform of Haar scaling function with translation (J) and dilation (k)
# Entry (n,m) is the Fourier transform at translation k[m] in xi[n]

function FourHaarScaling(xi::Vector, J::Int, k::Int)
	y = exp( -2.0*pi*im*2.0^(-J)*k*xi ) .* 2.0^(-J/2) .* FourHaarScaling(2.0^(-J)*xi)
end

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

function FourHaarWavelet(xi::Vector, J::Int, k::Int)
	y = exp( -2*pi*im*2.0^(-J)*k*xi ) .* 2.0^(-J/2) .* FourHaarWavelet(2.0^(-J)*xi);
end

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

# Rectangular domain xi_x by xi_y in the frequency plane

function freq2Haar( xi_x::Vector, xi_y::Vector, J::Int )
	Tx = freq2Haar(xi_x, J)
	Ty = freq2Haar(xi_y, J)

	return Tx, Ty
end


# General points in 2D
#
# Input:
# xi: M-by-2 matrix, where M is the number of observations

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

function freq2Haar( xi::AbstractMatrix, J::Int, weights::AbstractVector )
	@assert size(xi,1) == length(weights)

	T = freq2Haar( xi, J )
	freq2Haar( T, weights )

	return T
end


