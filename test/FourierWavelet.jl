using GeneralizedSampling
using Base.Test

# The Fourier transform of the scaling function of an MRA wavelet must 
# satisfy 'nested' relation (see e.g. (2.3) of Hernandez & Weiss)

# Filter
N = 2
C = wavefilter( string("db", N) )
scale!(C, 1/sum(C))

y1 = xi -> FourDaubScaling(xi, C) * GeneralizedSampling.feval(xi, C)
y2 = xi -> FourDaubScaling(2*xi, C)
y = xi -> abs( y1(xi) - y2(xi) )

# TODO: Would be better to measure in L2-norm, but quadgk takes forever
xi = collect( linspace(0,5,100) )
Y = map(xi -> y(xi), xi)
@test_approx_eq_eps( norm(Y), 0.0, sqrt(eps()) )

