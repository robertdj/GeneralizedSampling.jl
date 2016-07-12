using GeneralizedSampling
using Base.Test

println("Testing multiplication with Freq2Wave type...")

#=
Multiplication with a Freq2Wave element should give the same result as
multiplication with the collected matrix.

For both 1D and 2D samples and the 4 combinations of non-uniform samples
and wavelets with and without boundaries are tested with the
multiplications * and '*
=#

EPS = 1e-5

# ------------------------------------------------------------
# 1D

begin
	J = 3
	M = 2^(J+1)
	N = 2^J

	x = rand(N)
	v = rand(M)

	# ----------------------------------------
	# Compute matrices

	# Uniform points, No Boundary
	Usamples = grid(M, 0.5)
	TUNB = Freq2Wave(Usamples, "haar", J)
	AUNB = collect(TUNB)

	y1UNB = TUNB*x
	y2UNB = AUNB*x
	#@show norm(y1UNB - y2UNB)
	@test_approx_eq_eps y1UNB y2UNB EPS

	z1UNB = TUNB'*v
	z2UNB = AUNB'*v
	#@show norm(z1UNB - z2UNB)
	@test_approx_eq_eps z1UNB z2UNB EPS


	# Uniform points, with Boundary
	TUB = Freq2Wave(Usamples, "db2", J)
	AUB = collect(TUB)

	y1UB = TUB*x
	y2UB = AUB*x
	#@show norm(y1UB - y2UB)
	@test_approx_eq_eps y1UB y2UB EPS

	z1UB = TUB'*v
	z2UB = AUB'*v
	#@show norm(z1UB - z2UB)
	@test_approx_eq_eps z1UB z2UB EPS


	# Non-uniform points
	K = N/2 # bandwidth
	NUsamples = N*rand(M) - K

	# Non-Uniform points, No Boundary
	TNUNB = Freq2Wave(NUsamples, "haar", J, K)
	ANUNB = collect(TNUNB)

	y1NUNB = TNUNB*x
	y2NUNB = ANUNB*x
	#@show norm(y1NUNB - y2NUNB)
	@test_approx_eq_eps y1NUNB y2NUNB EPS

	z1NUNB = TNUNB'*v
	z2NUNB = ANUNB'*v
	#@show norm(z1NUNB - z2NUNB)
	@test_approx_eq_eps z1NUNB z2NUNB EPS


	# Non-Uniform points, with Boundary
	TNUB = Freq2Wave(NUsamples, "db2", J, K)
	ANUB = collect(TNUB)

	y1NUB = TNUB*x
	y2NUB = ANUB*x
	#@show norm(y1NUB - y2NUB)
	@test_approx_eq_eps y1NUB y2NUB EPS

	z1NUB = TNUB'*v
	z2NUB = ANUB'*v
	#@show norm(z1NUB - z2NUB)
	@test_approx_eq_eps z1NUB z2NUB EPS
end


# ------------------------------------------------------------
# 2D

begin
	J = 3
	M = 2^(J+1)
	N = 2^J

	x = rand(N^2)
	v = rand(M^2)

	# TODO: Run this in loop
	# ----------------------------------------
	# Compute matrices

	# Uniform points, No Boundary
	Usamples = grid((M,M), 0.5)
	TUNB = Freq2Wave(Usamples, "haar", J)
	AUNB = collect(TUNB)

	y1UNB = TUNB*x
	y2UNB = AUNB*x
	#@show norm(y1UNB - y2UNB)
	@test_approx_eq_eps y1UNB y2UNB EPS

	z1UNB = TUNB'*v
	z2UNB = AUNB'*v
	#@show norm(z1UNB - z2UNB)
	@test_approx_eq_eps z1UNB z2UNB EPS


	# Uniform points, with Boundary
	TUB = Freq2Wave(Usamples, "db2", J)
	AUB = collect(TUB)

	y1UB = TUB*x
	y2UB = AUB*x
	#@show norm(y1UB - y2UB)
	@test_approx_eq_eps y1UB y2UB EPS

	z1UB = TUB'*v
	z2UB = AUB'*v
	#@show norm(z1UB - z2UB)
	@test_approx_eq_eps z1UB z2UB EPS

	# Non-uniform points
	K = N/2 # bandwidth
	NUsamples = N*rand(M^2,2) - K

	# Non-Uniform points, No Boundary
	TNUNB = Freq2Wave(NUsamples, "haar", J, K)
	ANUNB = collect(TNUNB)

	y1NUNB = TNUNB*x
	y2NUNB = ANUNB*x
	#@show norm(y1NUNB - y2NUNB)
	@test_approx_eq_eps y1NUNB y2NUNB EPS

	z1NUNB = TNUNB'*v
	z2NUNB = ANUNB'*v
	#@show norm(z1NUNB - z2NUNB)
	@test_approx_eq_eps z1NUNB z2NUNB EPS


	# Non-Uniform points, with Boundary
	TNUB = Freq2Wave(NUsamples, "db2", J, K)
	ANUB = collect(TNUB)

	y1NUB = TNUB*x
	y2NUB = ANUB*x
	#@show norm(y1NUB - y2NUB)
	@test_approx_eq_eps y1NUB y2NUB EPS

	z1NUB = TNUB'*v
	z2NUB = ANUB'*v
	#@show norm(z1NUB - z2NUB)
	@test_approx_eq_eps z1NUB z2NUB EPS
end

