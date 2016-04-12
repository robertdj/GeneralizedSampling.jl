#= using GeneralizedSampling =#
using Base.Test

#=
Test collecting of Freq2Wave:
- The Fourier transform of the non-boundary part should be the same in a Freq2NoBoundaryWave and Freq2BoundaryWave
=#

# ------------------------------------------------------------
# 1D

J = 3
M = 2^(J+1)

samples = grid(M, 0.5)

# No boundary
T = freq2wave(samples, "haar", J)
A = collect(T)

# With boundary (using the same internal Fourier transform as above)
TT = freq2wave(samples, "db2", J)
Tb = Freq2BoundaryWave(TT.samples, T.FT, TT.weights, TT.J, TT.wavename, T.diag, TT.NFFT, TT.left, TT.right)
Ab = collect(Tb)

vm = van_moment(Tb)
# The internal parts are the same
N = size(T,2)
@test_approx_eq A[:,vm+1:N-vm] Ab[:,vm+1:N-vm]

