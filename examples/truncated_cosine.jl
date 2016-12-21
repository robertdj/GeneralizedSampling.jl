# ------------------------------------------------------------
# Truncated cosine and its Fourier transform

# Support = [-1/2,1/2]
function tcos(x)
	if 0.0 <= x <= 0.5
		return cos(2*pi*x)
	else
		return 0.0
	end
end

function ftcos(xi)
	if abs(xi) == 1.0
		return 0.25 + 0.0*im
	else
		return im*xi*( 1 + exp(-pi*xi*im) ) / (2*pi*(1 - xi^2))
	end
end


# Overload ftcos to vector input
@vectorize_1arg Float64 ftcos


# ------------------------------------------------------------
# Reconstruct truncated cosine in Haar basis

using GeneralizedSampling

J = 5
M = 2^(J+2)
# Both GeneralizedSampling and Plots (below) have a grid function
xi = GeneralizedSampling.grid(M, 0.5)
f = ftcos(xi)

T = Freq2Wave(xi, "haar", J)
wcoef = T \ f


# ------------------------------------------------------------
# Plot reconstruction

using IntervalWavelets

x, yw = weval( real(wcoef), "haar", 10 )

using Plots

p = plot(x, tcos, 
         linecolor=:red,
         linewidth=2,
         label="true"
        )

plot!(p, x, yw, 
      linecolor=:black,
      linewidth=2,
      label="approximation"
     )

display(p)

