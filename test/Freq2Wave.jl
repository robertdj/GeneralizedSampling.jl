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
TUNB = freq2wave(Usamples, "haar", J)

# Uniform samples, with Boundary
wavename = "db3"
TUB = freq2wave(Usamples, wavename, J)

# Non-Uniform samples, No Boundary
K = N/2
NUsamples = N*rand(M) - K
TNUNB = freq2wave(NUsamples, "haar", J, K)

# Non-Uniform samples, with Boundary
TNUB = freq2wave(NUsamples, wavename, J, K)


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

