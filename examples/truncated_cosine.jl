# ------------------------------------------------------------
# Truncated cosine and its Fourier transform

function tcos(x)
	if 0.5 <= x <= 1.0
		return cos(2.0*pi*x)
	else
		return 0.0
	end
end

function ftcos(xi)
	if abs(xi) == 1.0
		return complex(0.25)
	else
		return im*xi*( exp(-2.0*pi*xi*im) + exp(-pi*xi*im) ) / (2.0*pi*(xi^2 - 1.0))
	end
end

# Overload tcos and ftcos to vector input
@vectorize_1arg Float64 tcos
@vectorize_1arg Float64 ftcos


# ------------------------------------------------------------
# Reconstruct truncated cosine in Haar basis

using GeneralizedSampling

J = 6
M = 2^(J+1)
# Both GeneralizedSampling and Winston (below) have a grid function
xi = GeneralizedSampling.grid(M)
f = ftcos(xi)

T = Freq2Wave(xi, "haar", J)
wcoef = T \ f


# ------------------------------------------------------------
# Plot reconstruction

using IntervalWavelets
using Winston

x, yw = weval( real(wcoef), "haar", 10 )
plot(x,yw)

oplot(x, tcos(x), "r-")

