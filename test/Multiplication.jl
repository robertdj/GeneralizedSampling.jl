#= using GeneralizedSampling =#
using Base.Test

#=
Multiplication with a Freq2Wave element should give the same result as
multiplication with the collected matrix.

For both 1D and 2D samples and the 4 combinations of non-uniform samples
and wavelets with and without boundaries are tested with the
multiplications * and '*
=#

#= const EPS = 1e-5 =#
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
	TUNB = freq2wave(Usamples, "haar", J)
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
	TUB = freq2wave(Usamples, "db2", J)
	AUB = collect(TUB)

	y1UB = TUB*x
	y2UB = AUB*x
	#@show norm(y1UB - y2UB)
	@test_approx_eq_eps y1UB y2UB EPS

	z1UB = TUNB'*v
	z2UB = AUNB'*v
	#@show norm(z1UB - z2UB)
	@test_approx_eq_eps z1UB z2UB EPS


	# Non-uniform points
	K = N/2 # bandwidth
	NUsamples = N*rand(M) - K
	#= sort!(NUsamples) =#

	# Non-Uniform points, No Boundary
	TNUNB = freq2wave(NUsamples, "haar", J; B=K)
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
	TNUB = freq2wave(NUsamples, "db2", J; B=K)
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

#=
J = 3
M = 2^(J+1)

# ----------------------------------------
# Compute matrices

# Uniform points
samples1 = grid((M,M), 0.5)

# No boundary
T1 = freq2wave(samples1, "haar", J)
A1 = collect(T1)

# With boundary
T1b = freq2wave(samples1, "db2", J)
#= T = freq2wave(samples1, "db2", J) =#
#= T1b = Freq2BoundaryWave(T.samples, ones(T.FT), T.weights, T.J, T.wavename, ones(T.diag), T.NFFT, zeros(T.left), zeros(T.right)) =#
A1b = collect(T1b)


# Non-uniform points
N = size(T1,2)
K = N/2 # bandwidth
samples2 = N*rand(size(T1,1),2) - K

# No boundary
# TODO: This one is not being tested
T2 = freq2wave(samples2, "haar", J; B=K)
A2 = collect(T2)

# With boundary
T2b = freq2wave(samples2, "db2", J; B=K)
A2b = collect(T2b)


# ----------------------------------------
# Perform computations

x = rand(N)
y1 = T1*x
y2 = A1*x
#@show norm(y1 - y2)
@test_approx_eq_eps y1 y2 EPS

X = rand(wsize(T1b))
y1b = T1b*vec(X)
y2b = A1b*vec(X)
#@show norm(y1b - y2b)
@test_approx_eq_eps y1b y2b EPS

v = rand(size(T1,1))
z1 = T1'*v
z2 = A1'*v
@show norm(z1 - z2)
#@test_approx_eq_eps z1 z2 EPS

z1b = T2b'*v
z2b = A2b'*v
#@show norm(z1b - z2b)
@test_approx_eq_eps z1b z2b EPS
=#

