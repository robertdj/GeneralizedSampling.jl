#= using GeneralizedSampling =#
using Base.Test

println("Testing least squares solver...")

#=
Test the least squares solvers for Freq2wave types and their collected matrices.
=#

begin
	J = 10
	M = 2^(J+2)
	N = 2^J

	# Uniform samples
	xi = float( [-M/2:M/2-1;] )
	TU = freq2wave(xi, "Haar", J)
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
	TNU = freq2wave(xi, "Haar", J; B=K)
	ANU = collect(TNU)
	b = rand(M) + rand(M)*im

	# Avoid weighing b for ANU with "\"
	z0 = zeros(Complex{Float64}, N)
	z1 = cgnr( ANU, b, z0 )
	z2 = cgnr( TNU, b, z0 )

	#@show norm(z1 - z2, Inf)
	@test_approx_eq_eps z1 z2 1e-4
end

# 2D
begin
	J = 5
	M = 2^(J+1)

	# Uniform samples
	xi = grid((M,M), 0.5)
	TU = freq2wave(xi, "Haar", J; B=K)
	AU = collect(TU)
	b = rand(size(TU,1))

	x1 = AU \ b
	x2 = TU \ b

	#@show norm(x1 - vec(x2), Inf)
	@test_approx_eq_eps x1 x2 1e-5


	# Non-uniform samples: Solutions are far from identical. Collected
	# matrix has *high* condition number
end

