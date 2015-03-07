# Reconstruction of impulse responses, i.e., a single coefficient in the wavelet basis
# from equidistant frequency samples

# Input: 
# M: The number of frequency samples
# epsislon: The distance between two consecutive samples
# N: The number of wavelet coefficients recovered

M = 8
epsilon = 0.5
J = floor(Int, log2(M))

xi = [-M:epsilon:M;]

# Change of basis matrix
T = freq2Haar(xi, J)

# Fourier observations
b = FourHaarScaling(xi, J, 3)

# Wavelet coefficients
y = pinv(T)*b

#=function Haar_impulse(M::Int, epsilon::Float64, N::Int)=#
#=end=#

