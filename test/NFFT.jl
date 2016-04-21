#= using GeneralizedSampling =#
using Base.Test

println("Testing NFFT...")

#=
Test if the (N)FFT part of a Freq2Wave type behaves like a DFT matrix by
removing the parts related to wavelets.

For integer samples [0:N-1] the DFT matrix is Hermitian.
=#


# ------------------------------------------------------------
# 1D

begin
	N = 16.0
	samples = collect( 0:N-1 )

	J = Int(log2(N))
	T = freq2wave(samples, "haar", J)

	TF = Freq2NoBoundaryWave(T.samples, ones(T.FT), T.weights, T.J, T.wavename, ones(T.diag), T.NFFT)
	F = collect(TF)
	B = F'*F / N

	@test_approx_eq B eye(B)
end


# ------------------------------------------------------------
# 2D

begin
	M = 4
	samples = grid(M)

	J = Int(log2(M))
	T = freq2wave(samples, "haar", J)

	TF = Freq2NoBoundaryWave(T.samples, ones(T.FT), T.weights, T.J, T.wavename, ones(T.diag), T.NFFT)
	F = collect(TF)
	N = size(TF,2)
	B = F'*F / N

	@test_approx_eq B eye(B)
end

