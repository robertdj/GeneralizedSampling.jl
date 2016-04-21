# ------------------------------------------------------------
# Truncated cosine and its Fourier transform

function tcos(x::Float64)
	@assert 0 <= x <= 1
	x <= 0.5 ? 0.0 : cos(2*pi*x)
end

function ftcos(xi::Float64)
	if abs(xi) == 1.0
		return complex(0.25)
	else
		return im*xi*( exp(-2*pi*xi*im) + exp(-pi*xi*im) ) / (2*pi*(xi^2-1))
	end
end

@vectorize_1arg Float64 tcos
@vectorize_1arg Float64 ftcos


# ------------------------------------------------------------
# Reconstruct truncated cosine in Haar basis

using GeneralizedSampling

J = 6
N = 2^J
xi = float( collect(-N+1:N) )
f = ftcos(xi)

T = freq2wave(xi, "Haar", J)
coef = T \ f


# ------------------------------------------------------------
# Plot reconstruction

using PyPlot

x,y = weval( real(coef), 10 )
plot(x,y)

