using GeneralizedSampling
using Base.Test

println("Testing miscellaneous functions...")

#=
Testing the miscellaneous functions:
- had!, hadc!, yphad!, conj! should the same result as the standard not in-place functions.
- isuniform for matrices should return true for output from grid and false otherwise.
- In 2D the weights/Voronoi areas are known beforehand for some point configurations.
- For some point configurations, the density is known beforehand.
- Testing the relevant wavelet names for isdaubechies and ishaar.
- Assembling the spilt parts should yield the original array.
=#

# ------------------------------------------------------------
# had!, hadc!, yphad!, conj!

M = rand(2:10)
N = rand(2:10)

# had!
A = rand(M, N)
B = rand(M, N)
hadtmp = A.*B
had!(A, B)
@test hadtmp == A

# hadc!
A = rand(M, N) + im*rand(M,N)
B = rand(M, N) + im*rand(M,N)
hadctmp = A.*conj(B)
hadc!(A, B)
@test hadctmp == A

# yphad!
a = rand(M)
b = rand(M)
y = rand(M)
yphadtmp = y + a.*b
yphad!(y, a, b)
@test yphadtmp == y

# conj!
conj!(A, B)
@test A == conj(B)


# ------------------------------------------------------------
# isuniform and grid

for M in 15:16, scaling in rand(2)
	x = grid(M, scaling)
	@test isuniform(x)

	x[1] += 1
	@test !isuniform(x)

	@test !isuniform(rand(M))

	P = grid( (M,M), scaling )
	@test isuniform(P)

	P[1] += 1
	@test !isuniform(P)

	@test !isuniform(rand(M,2))
end
 

# ------------------------------------------------------------
# weights

# A corner in the "middle" of each quadrant
xi = [0.5 0.5 ; -0.5 0.5 ; -0.5 -0.5 ; 0.5 -0.5]
@test_approx_eq weights(xi, 1.0) ones(4)


# ------------------------------------------------------------
# Density

# 1D: Equidistant points
for N in 3:10
	x = linspace(-1,1,N)
	@test_approx_eq density(x,1) 2/(N-1)
end

@test density( [-1.0; 1.0], 1 ) == 2.0

# A corner in the "middle" of each quadrant
xi = [0.5 0.5 ; -0.5 0.5 ; -0.5 -0.5 ; 0.5 -0.5]
@test_approx_eq density(xi, 1.0) sqrt(2)/2


# ------------------------------------------------------------
# ishaar and isdaubechies

for N in 2:10, prefix in ["db", "DB", "Db", "dB"]
	name = string(prefix, N)
	@test isdaubechies(name)
	@test !ishaar(name)
end

@test ishaar("haar")
@test ishaar("db1")

@test !isdaubechies("ddb2")


# ------------------------------------------------------------
# Split 

# Vector
x = rand( rand(5:10) )
L, I, R = split(x,2)
@test x == vcat(L, I, R)

# Matrix
A = rand(rand(5:10), rand(5:10))
S = split(A,2)
@test A == [ S.LL S.LI S.LR ; S.IL S.II S.IR ; S.RL S.RI S.RR ]

