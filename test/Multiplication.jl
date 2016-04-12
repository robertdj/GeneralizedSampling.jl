#= using GeneralizedSampling =#
using Base.Test

#=
Multiplication with a Freq2Wave element should give the same result as
multiplication with the collected matrix.
=#

#= const EPS = 1e-5 =#
EPS = 1e-5

# ------------------------------------------------------------
# 1D

begin
	J = 3
	M = 2^(J+1)

	# ----------------------------------------
	# Compute matrices

	# Uniform points
	samples1 = grid(M, 0.5)

	# No boundary
	T1 = freq2wave(samples1, "haar", J)
	A1 = collect(T1)

	# With boundary
	T1b = freq2wave(samples1, "db2", J)
	A1b = collect(T1b)


	# Non-uniform points
	N = size(T1,2)
	K = N/2 # bandwidth
	samples2 = N*rand(M) - K
	sort!(samples2)

	# No boundary
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

	y1b = T1b*x
	y2b = A1b*x
	#@show norm(y1b - y2b)
	@test_approx_eq_eps y1b y2b EPS

	v = rand(M)
	z1 = T1'*v
	z2 = A1'*v
	#@show norm(z1 - z2)
	@test_approx_eq_eps z1 z2 EPS

	z1b = T2b'*v
	z2b = A2b'*v
	#@show norm(z1b - z2b)
	@test_approx_eq_eps z1b z2b EPS

end


# ------------------------------------------------------------
# 2D

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
A1b = collect(T1b)


# Non-uniform points
N = size(T1,2)
K = N/2 # bandwidth
samples2 = N*rand(size(T1,1),2) - K

# No boundary
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

y1b = T1b*x
y2b = A1b*x
@show norm(y1b - y2b)
#@test_approx_eq_eps y1b y2b EPS

v = rand(size(T1,1))
z1 = T1'*v
z2 = A1'*v
#@show norm(z1 - z2)
@test_approx_eq_eps z1 z2 EPS

#=
z1b = T2b'*v
z2b = A2b'*v
#@show norm(z1b - z2b)
@test_approx_eq_eps z1b z2b EPS
=#

