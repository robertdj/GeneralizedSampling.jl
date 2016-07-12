using GeneralizedSampling
import IntervalWavelets: unit
using Base.Test

println("Testing least squares solver...")

#=
Test the least squares solvers:
- Solution with Freq2Wave and the collected version should be identical.
- When sampling from a scaling function used in the Freq2Wave type, the solution should be a unit vector.
=#


# ------------------------------------------------------------
# Reconstruct unit vector

J = 3
M = 2^(J+1)
N = 2^J
xi = grid(M, 0.5)
wavename = "db2"
vm = van_moment(wavename)

T = Freq2Wave(xi, wavename, J)

for k = 1:N
	if k <= vm
		F = FourScalingFunc(xi, wavename, 'L', J)
		f = F[:,k]
	elseif k > N-vm
		F = FourScalingFunc(xi, wavename, 'R', J)
		f = F[:,k-N+vm]
	else
		trans = k - div(N,2) - 1
		f = FourDaubScaling(xi, vm, J, trans)
	end
	y = T \ f

	u = unit(Complex{Float64}, length(y), k)
	#@show k, norm( y - u )
	@test_approx_eq_eps y u 1e-4
end


# ------------------------------------------------------------

begin
	J = 10
	M = 2^(J+2)
	N = 2^J

	# Uniform samples
	xi = grid(M, 0.5)
	TU = Freq2Wave(xi, "Haar", J)
	AU = collect(TU)

	b = rand(M)
	x0 = rand(N)

	x1 = AU \ b
	x2 = TU \ b

	#@show norm(x1 - x2, Inf)
	@test_approx_eq_eps x1 x2 1e-5


	# Non-uniform samples: Requires high M/N ratio 
	K = N/2
	xi = 2*K*rand(M) - K
	sort!(xi)
	TNU = Freq2Wave(xi, "Haar", J, K)
	ANU = collect(TNU)
	b = rand(M) + rand(M)*im

	bb = b .* get(TNU.weights)
	z1 = ANU \ bb
	z2 = TNU \ b

	#@show norm(z1 - z2, Inf)
	@test_approx_eq_eps z1 z2 1e-4
end

# 2D
begin
	J = 5
	M = 2^(J+1)
	N = 2^J

	# Uniform samples
	xi = grid((M,M), 0.5)
	TU = Freq2Wave(xi, "Haar", J)
	AU = collect(TU)
	b = rand(size(TU,1))

	x1 = AU \ b
	x2 = TU \ b

	#@show norm(x1 - vec(x2), Inf)
	@test_approx_eq_eps x1 x2 1e-5

	# Non-uniform samples: Requires high M/N ratio 
	K = N/2
	xi = 2*K*rand(M^2,2) - K
	TNU = Freq2Wave(xi, "Haar", J, K)
	ANU = collect(TNU)
	b = rand(M^2) + rand(M^2)*im

	bb = b .* get(TNU.weights)
	z1 = ANU \ bb
	z2 = TNU \ b

	#@show norm(z1 - z2, Inf)
	@test_approx_eq_eps z1 z2 1e-4
end

