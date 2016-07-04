# ------------------------------------------------------------
# Truncated cosine and its Fourier transform

# [-1/2,1/2]
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
		return im*xi*( 1.0 + exp(-pi*xi*im) ) / (2.0*pi*(1 - xi^2))
	end
end

#=
# [-1,1]
function tcos(x)
	if 0.0 <= x <= 1.0
		return cos(pi*x)
	else
		return 0.0
	end
end

# [-1,1] wrt exp(-2pi*im*xi*x)
#= function ftcos(xi) =#
#= 	if abs(xi) == 0.5 =#
#= 		return 0.5 + 0.0*im =#
#= 	else =#
#= 		return 2*im*xi*( 1 + exp(-2*im*pi*xi) ) / (pi*(1 - 4*xi^2)) =#
#= 	end =#
#= end =#

# [-1,1] wrt exp(-pi*im*xi*x)
function ftcos(xi)
	if abs(xi) == 1.0
		return 0.5 + 0.0*im
	else
		return im*xi*( 1 + exp(-im*pi*xi) ) / (pi*(1 - xi^2))
	end
end
=#

# Overload tcos and ftcos to vector input
@vectorize_1arg Float64 tcos
@vectorize_1arg Float64 ftcos


# ------------------------------------------------------------
# Reconstruct truncated cosine in Haar basis

using GeneralizedSampling

J = 5
M = 2^(J+2)
# Both GeneralizedSampling and Winston (below) have a grid function
xi = GeneralizedSampling.grid(M, 0.5)
f = ftcos(xi)

T = Freq2Wave(xi, "haar", J)
#= T = Freq2Wave(xi, "db2", J) =#
wcoef = T \ f
#= wcoef, h = lsqr(T, f) =#
A = collect(T)
#= wcoef = A \ f =#


# ------------------------------------------------------------
# Plot reconstruction

using IntervalWavelets
using Winston

x, yw = weval( real(wcoef), "haar", 10 )
#= x, yw = weval( real(wcoef), "db2", 10 ) =#

plot(x, yw, x, tcos(x))

