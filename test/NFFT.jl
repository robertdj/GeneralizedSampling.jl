using GeneralizedSampling
using Base.Test

println("Testing NFFT...")

#=
Test if the (N)FFT part of a Freq2Wave type behaves like a DFT matrix by
removing the parts related to wavelets.

For integer samples [0:N-1] the DFT matrix is Hermitian.

The special version of NFFT with uniform sampling locations should agree 
with the general version.
=#


# ------------------------------------------------------------
# 1D

begin
	N = 16.0
	samples = collect( 0:N-1 )

	J = Int(log2(N))
	T = Freq2Wave(samples, "haar", J)
	for idx in eachindex(T.internal)
		T.internal[idx] = 1.0 + 0.0*im
	end

	F = collect(T)
	B = F'*F / N

	@test_approx_eq B eye(B)
end


# ------------------------------------------------------------
# 2D

begin
	M = 4
	samples = grid(M)

	J = Int(log2(M))
	T = Freq2Wave(samples, "haar", J)
	for idx in eachindex(T.internal)
		T.internal[idx] = 1.0 + 0.0*im
	end

	F = collect(T)
	N = size(T,2)
	B = F'*F / N

	@test_approx_eq B eye(B)
end


# ------------------------------------------------------------
# Uniform NFFT

J = rand(2:5)
M = 2*rand( 2^(J-1):2^J )
epsilon = 1/rand(2:5)
xi = GeneralizedSampling.grid(M, epsilon)

using NFFT

N1 = 2*rand( 2:2^(J-1) )
chi = xi * 2.0^-J
P1 = NFFTPlan( frac(chi), N1 )

f = rand(N1) + rand(N1)*im
y1 = nfft(P1, f)

P2 = GeneralizedSampling.FFTPlan( xi, J, N1 )
y2 = similar(y1)
GeneralizedSampling.myfft!(y2, P2, f)

@test_approx_eq_eps y1 y2 1e-5

