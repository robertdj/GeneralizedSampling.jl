using GeneralizedSampling
import IntervalWavelets: coef, inner, l2norm, ifilter
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
@test_approx_eq_eps FourHaarScaling(xi) FourDaubScaling(xi, 1) 1e-7


# ------------------------------------------------------------
# Internal scaling functions

# For p == 2 and k == 1 the inner product is only accurate within 3
# decimals. For all other p/k combinations it is better (as expected due
# to the low regularity of that scaling function)
EPS = 1e-3

for p in 2:8
	C = coef(ifilter(p))
	scale!(C, 1/sum(C))

	x = linspace(-7,7,1001)
	y = FourDaubScaling(x, C)
	#@show l2norm(x, y)
	@test_approx_eq_eps l2norm(x,y) 1.0 EPS

	for k in 1:p+1
		z = FourDaubScaling(x, C, 0, k)
		#@show inner(x, y, z)
		@test_approx_eq_eps inner(x, y, z) 0.0 EPS
	end
end


# ------------------------------------------------------------
# Boundary scaling functions

# In some cases a low precision is needed
EPS = 1e-1
for side = ['L'; 'R']
	for p = 2:8
		x = linspace(-15,15,3001)

		Y = FourDaubScaling(x, p, side)
		for l in 1:p
			#@show p, l, l2norm(x, vec(Y[l,:]))
			@test_approx_eq_eps l2norm(x, vec(Y[l,:])) 1.0 EPS

			# Boundary-boundary
			for k in l+1:p
				IP = inner(x, vec(Y[l,:]), vec(Y[k,:]))
				#@show p, k, IP
				@test_approx_eq_eps IP 0.0 EPS
			end

			# Boundary-internal for appropriate translation of internal
			C = coef(ifilter(p))
			scale!(C, 1/sum(C))
			transl = (side == 'L' ? p : -p-1 )
			IP = inner(x, vec(Y[l,:]), FourDaubScaling(x, p, 0, transl))
			#@show p, l, IP
			@test_approx_eq_eps IP 0.0 EPS
		end
	end
end

