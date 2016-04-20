#= using GeneralizedSampling =#
using Base.Test

#=
The Fourier transform of the scaling function of an MRA wavelet must 
satisfy 'nested' relation (see e.g. (2.3) of Hernandez & Weiss)

Thanks to Parseval we know that the Fourier transforms of the scaling
functions (boundary and internal) are orthogonal in L2.
=#


# ------------------------------------------------------------
# Fourier transform of the boundary scaling functions

N = 2
F = scalingfilters(N)
xi = linspace(-10,10,1001)
y = Array(Complex{Float64}, length(xi), N)
count = 0
for x in xi
	y[count+=1,:] = FourDaubScaling(x, F)
end

