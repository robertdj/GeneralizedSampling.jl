using GeneralizedSampling
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
	T = Freq2Wave(samples, "haar", J)
	A = collect(T)

	# With boundary (using the same internal Fourier transform as above)
	TMP = Freq2Wave(samples, "db2", J)
	Tb = GeneralizedSampling.Freq2BoundaryWave1D(T.internal, TMP.weights, TMP.J, TMP.wavename, TMP.NFFT, TMP.left, TMP.right)
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
	T = Freq2Wave(samples, "haar", J)
	A = collect(T)

	# With boundary (using the same internal Fourier transform as above)
	TMP = Freq2Wave(samples, "db2", J)
	Tb = GeneralizedSampling.Freq2BoundaryWave2D(T.internal, TMP.weights, TMP.J, TMP.wavename, TMP.NFFT, TMP.left, TMP.right)
	Ab = collect(Tb)

	# Select the indices of the internal scaling functions
	vm = van_moment(Tb)
	J = zeros(Bool, wsize(T))
	S = split(J,vm)
	for idx in eachindex(S.II)
		S.II[idx] = true
	end

	#@show norm( A[:,vec(J)] - Ab[:,vec(J)] )
	@test_approx_eq A[:,vec(J)] Ab[:,vec(J)]
end

