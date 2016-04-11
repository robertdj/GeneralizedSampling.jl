using GeneralizedSampling
using Base.Test

#=
Test if the (N)FFT part of a Freq2Wave type behaves like a DFT matrix by
removing the parts related to wavelets.
=#


# 1D

N = 16.0
samples = collect( 0:N-1 )

J = Int(log2(N))
T = freq2wave(samples, "haar", J)

TF1 = Freq2NoBoundaryWave(T.samples, ones(T.FT), T.weights, T.J, T.wavename, ones(T.diag), T.NFFT)
F1 = collect(TF1)
B1 = F1'*F1 / N

@test_approx_eq B eye(B)


# 2D

M = 4
samples = grid(M)

J = Int(log2(M))
T = freq2wave(samples, "haar", J)

TF2 = Freq2NoBoundaryWave(T.samples, ones(T.FT), T.weights, T.J, T.wavename, ones(T.diag), T.NFFT)
F2 = collect(TF2)
N = size(TF2,2)
B2 = F2'*F2 / N

@test_approx_eq B eye(B)

