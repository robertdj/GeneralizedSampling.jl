using GeneralizedSampling
using Base.Test

J = 10
M = 2^(J+2)
N = 2^J

# Uniform samples
xi = float( [-M/2:M/2-1;] )
T = freq2wave(xi, "Haar", J)
A = collect(T)

b = complex(rand(M))
x0 = complex(rand(N))

x = A \ b

x1 = cgnr(A, b, x0)
x2 = cgnr(T, b, x0)

#@show norm(x - x1, Inf)
#@show norm(x - x2, Inf)
@test_approx_eq_eps( x, x1, 1e-5 )
@test_approx_eq_eps( x, x2, 1e-5 )


# Non-uniform samples: Requires high M/N ratio 
K = N/2
xi = 2*K*rand(M) - K
sort!(xi)
T = freq2wave(xi, "Haar", J; B=K)
A = collect(T)
b = rand(M) + rand(M)*im

z = A \ b;

z1 = cgnr(A, b, x0)
z2 = cgnr(T, b, x0)

#@show norm(z - z1, Inf)
#@show norm(z - z2, Inf)
@test_approx_eq_eps( z, z1, 1e-5 )
@test_approx_eq_eps( z, z2, 1e-4 )

