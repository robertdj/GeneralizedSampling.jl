using GeneralizedSampling
using Base.Test

#=
Various small functions are defined to interact with a Freq2Wave type.
=#


# ------------------------------------------------------------
# Creating matrices

J = 3
M = 2^(J+1)
N = 2^J

# Uniform samples, No Boundary
Usamples = grid(M, 0.5)
TUNB = Freq2Wave(Usamples, "haar", J)

# Uniform samples, with Boundary
wavename = "db3"
TUB = Freq2Wave(Usamples, wavename, J)

# Non-Uniform samples, No Boundary
K = N/2
NUsamples = N*rand(M) - K
TNUNB = Freq2Wave(NUsamples, "haar", J, K)

# Non-Uniform samples, with Boundary
TNUB = Freq2Wave(NUsamples, wavename, J, K)


# ------------------------------------------------------------
# Testing

@test isuniform(TUNB)
@test isuniform(TUB)
@test !isuniform(TNUNB)
@test !isuniform(TNUB)

@test !hasboundary(TUNB)
@test hasboundary(TUB)
@test !hasboundary(TNUNB)
@test hasboundary(TNUB)

@test dim(TUNB) == dim(TUB) == dim(TNUNB) == dim(TNUB) == 1

@test wscale(TUNB) == wscale(TUB) == wscale(TNUNB) == wscale(TNUB) == J

@test wsize(TUNB) == wsize(TUB) == wsize(TNUNB) == wsize(TNUB) == (N,)

@test size(TUNB) == size(TUB) == size(TNUNB) == size(TNUB) == (M,N)

@test van_moment(TUNB) == van_moment(TNUNB) == 1
@test van_moment(TUB) == van_moment(TNUB) == van_moment(wavename)

