# M-by-M sampling points
M = 16
xi = grid(M, M, 0.2)

# Reconstruction coefficients on N-by-N  grid
N = 8
x = complex(rand(N^2))
X = reshape_view(x, (N,N))

# NFFT package: 'direct' NDFT and fast NDFT
p = NFFTPlan(xi', (N,N))
a1 = ndft(p, X)
#= a2 = nfft(p, X) =#

# NDFT matrix
F = NDFT( xi[:,1], xi[:,2], N );
b1 = F*x

#= @show norm(a1 - b1, Inf) =#
@test_approx_eq_eps( a1, b1, 1e-10 )
#= @show norm(a2 - b1, Inf) =#
#= @test_approx_eq_eps( a2, b1, 1e-5 ) =#

