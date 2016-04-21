#= using GeneralizedSampling =#
using Base.Test

println("Testing collection of Freq2Wave type...")

#=
Test collecting of Freq2Wave:
- The Fourier transform of the non-boundary part should be the same in a Freq2NoBoundaryWave and Freq2BoundaryWave
- TODO: Test the boundary parts
=#

# ------------------------------------------------------------
# 1D

begin
	J = 3
	M = 2^(J+1)

	samples = grid(M, 0.5)

	# No boundary
	T = freq2wave(samples, "haar", J)
	A = collect(T)

	# With boundary (using the same internal Fourier transform as above)
	TMP = freq2wave(samples, "db2", J)
	Tb = Freq2BoundaryWave(TMP.samples, T.FT, TMP.weights, TMP.J, TMP.wavename, T.diag, TMP.NFFT, TMP.left, TMP.right)
	Ab = collect(Tb)

	vm = van_moment(Tb)
	# The internal parts are the same
	N = size(T,2)
	@test_approx_eq A[:,vm+1:N-vm] Ab[:,vm+1:N-vm]
end


# ------------------------------------------------------------
# 2D

begin 
	J = 3
	M = 2^(J+1)

	samples = grid((M,M), 0.5)

	# No boundary
	T = freq2wave(samples, "haar", J)
	A = collect(T)

	# With boundary (using the same internal Fourier transform as above)
	TMP = freq2wave(samples, "db2", J)
	Tb = Freq2BoundaryWave(TMP.samples, T.FT, TMP.weights, TMP.J, TMP.wavename, T.diag, TMP.NFFT, TMP.left, TMP.right)
	Ab = collect(Tb)

	# Select the indices of the internal scaling functions
	vm = van_moment(Tb)
	I = zeros(Bool, wsize(T))
	S = split(I,vm)
	for idx in eachindex(S.II)
		S.II[idx] = true
	end

	# @show norm( A[:,vec(I)] - Ab[:,vec(I)] )
	@test_approx_eq A[:,vec(I)] Ab[:,vec(I)]
end

