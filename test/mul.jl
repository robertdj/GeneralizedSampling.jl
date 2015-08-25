include("../src/FourierTransforms.jl")
include("../src/misc.jl")
using NFFT 

M = 2^4
xi = collect(linspace(-0.4, 0.4, M)) # nodes at which the NFFT is evaluated
xi = rand(M)*10

J = 3
N = 2^J
k = [0:2^J-1;]

b = rand(N) + rand(N)*im
c = rand(M) + rand(M)*im

y1, z1 = mul(xi, b, c)
y2, z2 = mul2(xi, J, k, b, c)

@show maximum(abs(y1 - y2))
@show maximum(abs(z1 - z2))

