using GeneralizedSampling
using Base.Test

xi = collect(-5:4)
dens = density( xi, max(abs(xi)) )
@test dens == 1

