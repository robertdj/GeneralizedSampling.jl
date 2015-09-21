# ------------------------------------------------------------
# Reconstruction of impulse responses, i.e., a single coefficient in the Haar 
# basis from equidistant frequency samples


# Define scaling function and mother wavelet

function Haar_scaling(x::Real)
	if x >= 0 && x <= 1.0
		return 1.0
	else
		return 0.0
	end
end

function Haar_wavelet(x::Real)
	if x >= 0 && x < 0.5
		return 1.0
	elseif x >= 0.5 && x <= 1.0
		return -1.0
	else
		return 0.0
	end
end

# ------------------------------------------------------------
# Reconstruction of scaling function

# Input: 
# M: The number of frequency samples
# epsislon: The distance between two consecutive samples
# N: The number of wavelet coefficients recovered

M = 32
epsilon = 1
#J = floor(Int, log2(M))
J = 2

# Sampling points
xi = [-M:M;]*epsilon

# Change of basis matrix
Tx, Ty = freq2Haar( xi, xi, J )
T = kron(Ty, Tx)

# Fourier observations
bx = FourHaarScaling(xi, J, 3)
by = FourHaarScaling(xi, J, 3)
b = bx * by.'

# Scaling function coefficients: 1 in the last (16th) entry, 0 elsewhere
y = pinv(T)*b[:] # Brute force
yy = pinv(Tx)*b*pinv(Ty).' # With Kronecker products

