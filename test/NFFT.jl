using GeneralizedSampling
using Base.Test

M = 2^8
xi = collect(linspace(-0.4, 0.4, M))
xi = rand(M)*10

J = 6

# Collect change of basis matrix
T = freq2wave(xi, "Haar", J)

A = collect(T)
B = freq2Haar(xi, J)

@test_approx_eq_eps( A, B, sqrt(eps()) )

# Multiplication with change of basis matrix
N = 2^J
b = rand(N) + rand(N)*im
y1 = T*b
y2 = A*b
@test_approx_eq_eps(y1, y2, 1e-4) 

b = rand(N)
y1 = T*b
y2 = A*b
@test_approx_eq_eps(y1, y2, 1e-4) 

# Multiplication with adjoint change of basis matrix
c = rand(M) + rand(M)*im
z1 = H(T, c)
z2 = A'*c
@test_approx_eq_eps(z1, z2, 1e-4) 

c = rand(M)
z1 = H(T, c)
z2 = A'*c
@test_approx_eq_eps(z1, z2, 1e-4) 

