# ------------------------------------------------------------
# Reconstruction of impulse responses, i.e., a single coefficient in the Haar 
# basis from equidistant frequency samples


# ------------------------------------------------------------
# Reconstruction of scaling function

# Input: 
# M: The number of frequency samples
# epsislon: The distance between two consecutive samples
# N: The number of wavelet coefficients recovered

M = 8
epsilon = 0.25
J = floor(Int, log2(M))

# Sampling points
xi = [-M:epsilon:M;]

# Change of basis matrix
T = freq2Haar(xi, J)

# Fourier observations
b = FourHaarScaling(xi, J, 3)

# Scaling function coefficients: 1 in the 4th entry, 0 elsewhere
y = pinv(T)*b


# ------------------------------------------------------------
# Reconstruction of wavelet function

J2 = [0:2;]

# Change of basis matrix
T2 = freq2Haar(xi, J2)

# Fourier observations
b2 = FourHaarWavelet(xi, 2, 1)
#b2 = FourHaarScaling(xi, 1, 0)

# Wavelet coefficients: 1 in the 5th entry, 0 elsewhere
y2 = pinv(T2)*b2

