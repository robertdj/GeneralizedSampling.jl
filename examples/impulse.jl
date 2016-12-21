# ------------------------------------------------------------
# Reconstruction of impulse responses

using GeneralizedSampling

J = 3
M = 2^(J+1)
wname = "db2"

xi = GeneralizedSampling.grid(M, 0.5)
f = FourScalingFunc(xi, wname, J)

T = Freq2Wave(xi, wname, J)
wcoef = T \ f

println("Reconstructing a ", wname, " function from its Fourier measurements.")
println("Reconstructed coefficients is a unit vector:")
display(wcoef)


# ------------------------------------------------------------
# Plot reconstruction

using IntervalWavelets

x, yw = weval( real(wcoef), wname, 10 )

using Plots

p = plot(x, yw, label=wname)
display(p)

