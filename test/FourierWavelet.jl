using GeneralizedSampling
import WaveletPlot: inner, ifilter
using Base.Test

println("Testing Fourier transforms...")

#=
- The Fourier transform of the Haar wavelet should be equal to that of the
Daubechies 1 wavelet, which is computed by a different function.

- Thanks to Parseval we know that the Fourier transforms of the scaling
functions (boundary and internal) are orthonormal in L2.
=#

# ------------------------------------------------------------
# Haar and Daubechies 1
xi = randn(10)
@test_approx_eq_eps FourHaarScaling(xi) FourDaubScaling(xi,1) 1e-7


# ------------------------------------------------------------
# Internal scaling functions

# For N == 2 and k == 1 the inner product is only accurate within 3
# decimals. For all other N/k combinations it is better.
EPS = 1e-3

for N in 2:8
	C = ifilter(N)
	scale!(C, 1/sum(C))

	x = linspace(-7,7,1001)
	y = FourDaubScaling(x, C)

	for k in 0:N+1
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
# Boundary scaling functions

# Large integration intervals are needed; instead a low precision is accepted.
EPS = 1e-1
# TODO: Right side boundary scaling functions seem to require a *huge* interval
#= for side = ['L'; 'R'] =#
	side = 'L'
	for N = 3:8
		if N == 3
			x = linspace(-50,50,2001)
		else 
			x = linspace(-30,30,2001)
		end

		Y = FourDaubScaling(x, N, side)

		for k in 1:N
			# Boundary-boundary
			I = inner(x, Y[:,1], Y[:,k])
			if k == 1
				@test_approx_eq_eps I 1.0 EPS
			else
				@test_approx_eq_eps I 0.0 EPS
			end

			# Boundary-internal
			C = ifilter(N)
			scale!(C, 1/sum(C))
			@test_approx_eq_eps inner(x, Y[:,k], FourDaubScaling(x, C)) 0.0 EPS
		end
	end
#= end =#

