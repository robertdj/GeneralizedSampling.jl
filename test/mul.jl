using GeneralizedSampling
using Base.Test

for J = 3:10
	#println(J)
	M = 2^(J+1)
	N = 2^J

	# Test uniform points
	xi = float( [-M/2:M/2-1;] )
	T = freq2wave(xi, "Haar", J)
	A = collect(T)

	x = complex( rand(N) )
	y1 = A*x
	y2 = T*x
	#@show norm(y1-y2)
	@test_approx_eq_eps( y1, y2, 1e-5 )

	v = complex( rand(M) )
	z1 = A'*v
	z2 = T'*v
	#@show norm(z1-z2)
	@test_approx_eq_eps( z1, z2, 1e-5 )


	# Test non-uniform points
	K = N/2
	xi = 2*K*rand(M) - K
	sort!(xi)
	T = freq2wave(xi, "Haar", J; B=K)
	A = collect(T)

	y1 = A*x
	y2 = T*x
	#@show norm(y1-y2)
	@test_approx_eq_eps( y1, y2, 1e-5 )
end
