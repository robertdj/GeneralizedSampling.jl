#= using GeneralizedSampling =#
#= import WaveletPlot: inner =#
using Base.Test

println("Testing Fourier transforms...")

#=
- The Fourier transform of the Haar wavelet should be equal to that of the
Daubechies 1 wavelet, which is computed by a different function.

- The Fourier transform of the scaling function of an MRA wavelet must 
satisfy the 'nested' relation (see e.g. (2.3) of Hernandez & Weiss)

- Thanks to Parseval we know that the Fourier transforms of the scaling
functions (boundary and internal) are orthonormal in L2.
=#

# ------------------------------------------------------------
# Haar and Daubechies 1
xi = randn(10)
@test_approx_eq_eps FourHaarScaling(xi) FourScalingFunc(xi,"db1") 1e-7


# ------------------------------------------------------------
# Internal scaling functions

EPS = 1e-3

for N = 2:8
	C = ifilter(N)
	scale!(C, 1/sum(C))

	x = linspace(-7,7,1001)
	y = FourDaubScaling(x, C)
	for k = 0:N+1
		z = FourDaubScaling(x, C, 0, k)
		I = inner(x, y, z)
		if k == 0
			@test_approx_eq_eps I 1.0 EPS
		else
			@test_approx_eq_eps I 0.0 EPS
		end
	end
end

# ------------------------------------------------------------

# ------------------------------------------------------------

# ------------------------------------------------------------
# Fourier transform of the boundary scaling functions

