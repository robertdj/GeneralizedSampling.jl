#= using GeneralizedSampling =#
using Base.Test

#=
Testing the misc functions:
- Assembling the spilt parts should yield the original array.
=#

x = rand(rand(5:10))
L, I, R = split(x,2)
@test x == vcat(L, I, R)

A = rand(rand(5:10), rand(5:10))
S = split(A,2)

@test A == [ S.LL S.LI S.LR ; S.IL S.II S.IR ; S.RL S.RI S.RR ]
