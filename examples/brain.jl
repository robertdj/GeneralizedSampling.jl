#= 
Reconstruct "brain" from Fourier measurements.
The measurements are generated with the Matlab files from
http://bigwww.epfl.ch/algorithms/mriphantom
=#

# ------------------------------------------------------------
# Load data

using JLD
data_file = joinpath(dirname(@__FILE__), "brain512.jld")
F = load( data_file )
f = F["freqs"]

# Sampling locations
M = size( f )
J = Int( log2(M[1]) ) - 1


# ------------------------------------------------------------
# Compute reconstruction

using GeneralizedSampling

xi = GeneralizedSampling.grid(M)
T = Freq2Wave( xi, "db4", J )
w = T \ f


# ------------------------------------------------------------
# Plot reconstruction

using IntervalWavelets

y = weval( abs(w), "db4", 10 )

using Plots

heatmap( y )

