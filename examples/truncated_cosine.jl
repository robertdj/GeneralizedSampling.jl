# ------------------------------------------------------------
# Truncated cosine and its Fourier transform

function tcos(x::Float64)
	@assert 0.0 <= x <= 1.0
	x <= 0.5 ? 0.0 : cos(2.0*pi*x)
end

function ftcos(xi::Float64)
	if abs(xi) == 1.0
		return complex(0.25)
	else
		return im*xi*( exp(-2.0*pi*xi*im) + exp(-pi*xi*im) ) / (2.0*pi*(xi^2-1))
	end
end

@vectorize_1arg Float64 tcos
@vectorize_1arg Float64 ftcos


# ------------------------------------------------------------
# Reconstruct truncated cosine in Haar basis

using GeneralizedSampling

J = 6
M = 2^(J+1)
xi = GeneralizedSampling.grid(M)
f = ftcos(xi)

T = freq2wave(xi, "Haar", J)
wcoef = T \ f


# ------------------------------------------------------------
# Plot reconstruction

using WaveletPlot
using Winston

x, yw = weval( real(wcoef), 10 )
plot(x,yw)

oplot(x, tcos(x), "r-")

