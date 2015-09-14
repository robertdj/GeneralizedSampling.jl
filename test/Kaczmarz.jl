using GeneralizedSampling
using Base.Test

J = 3
M = 2^(J+3)
N = 2^J

# Uniform Kaczmarz
xi = float( [-M/2:M/2-1;] )
T = freq2wave(xi, "Haar", J)
A = collect(T)
b = rand(M) + rand(M)*im
x = A \ b

x0 = complex(zeros(N))
prec = 1e-4
x1 = REK(T, b, x0; prec=prec)

#@show norm(x - x1, Inf)
@test_approx_eq_eps( x, x1, 0.01 )


# Non-uniform Kaczmarz
K = N/2
xi = 2*K*rand(M) - K
sort!(xi)
T = freq2wave(xi, "Haar", J; B=K)
A = collect(T)
b = rand(length(xi)) + rand(length(xi))*im
x = A \ b;

x1 = REK(T, b, x0; prec=prec, maxiter=200)
x2 = REK(A, b, x0; prec=prec, maxiter=200)

#@show norm(x - x1, Inf)
#@show norm(x1 - x2, Inf)
@test_approx_eq_eps( x, x1, 0.01 )
@test_approx_eq_eps( x1, x2, 0.01 )

