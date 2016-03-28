#= using GeneralizedSampling =#
using Base.Test

# Dyadic rationals
#=
upper = 3
res = 4
dy_res = dyadic_rationals(upper, res)

for level = 1:res
	idx = dyadic_rationals(upper, res, level)
	dy_level = dy_res[idx]
	@test isinteger( 2^level*dy_level )
end
=#

# ------------------------------------------------------------
# Translation between index and values for dyadic rationals

vm = 2
F = scalingfilters(vm)
IS, BS = support(F)

for support_type in [:IS, :BS]
	@eval begin
		# Level 0 = integers
		supportx = dyadic_rationals( $support_type, 0 )
		for x in round(Int, supportx)
			@test supportx[ x2index(x,$support_type) ] == x
		end

		for idx = 1:length(supportx)
			@test supportx[idx] == index2x(idx,$support_type)
		end

		# Higher level
		R = 2
		supportx = dyadic_rationals( $support_type, R )
		for x in supportx
			@test supportx[ x2index(x,$support_type,R) ] == x
		end

		for idx = 1:length(supportx)
			@test supportx[idx] == index2x(idx,$support_type,R)
		end
	end
end


# ------------------------------------------------------------
# Boundary wavelets satisfy dilation equation

vm = 4
F = scalingfilters(vm)
IS, BS = support(F)

R = 2

phi = DaubScaling(F.internal, R)[2]

xval, Y = DaubScaling(F, R, Val{'L'})
for x in xval
	doublex = 2*x

	for k = 0:vm-1
		curF = bfilter(F.left,k)
		curF_length = length(curF)

		doubleY = zeros(Float64, curF_length)
		if isinside(doublex, BS)
			doubleY[1:vm] = Y[x2index(doublex,BS,R),:]
		end

		for m = vm:curF_length-1
			if isinside(doublex-m, IS)
				doubleY[m+1] = phi[x2index(doublex-m,IS,R)]
			end
		end

		@test_approx_eq_eps sqrt(2)*dot(doubleY,curF) Y[x2index(x,BS,R),k+1] sqrt(eps())
		#= println("x = ", x, ", ", sqrt(2)*dot( doubleY, curF ) - Y[x2index(x,BS,R),k+1], ", ", sqrt(2)*dot( doubleY, curF ), ", ", Y[x2index(x,BS,R),k+1]) =#
	end
end

#=
function trapezquad(x::Vector, y::Vector)
	@assert 2 <= (Nx = length(x)) == length(y)
	@assert issorted(x)

	integral = 0.0
	for n = 2:Nx
		integral += (x[n]-x[n-1])*0.5*(y[n]+y[n-1])
	end

	return integral
end

# Test trapezquad
#= x = dyadic_rationals( (-1,1), 3 ) =#
#= @test trapezquad( x, ones(length(x)) ) == 2 =#
#= @test trapezquad( x, x ) == 0 =#

# The L2 inner product of y and z over the x's
function inner(x::Vector, y::Vector, z::Vector)
	trapezquad( x, y.*conj(z) )
end

function l2norm(x::Vector, y::Vector)
	trapezquad( x, abs2(y) )
end


# ------------------------------------------------------------
# Orthogonality of integer translates of the internal scaling function

R = 10
N = 3
F = scalingfilters(N)

# Union of supports of all translated scaling functions
supp = support(F.internal)
extended_supp = (supp[1], 2*supp[2]-1)
x = dyadic_rationals( extended_supp, R )

# Base scaling funtion
y1 = zeros(size(x))
phi = DaubScaling(F.internal, R)[2]
supp_index = [1:length(phi);]
y1[supp_index] = phi

# Integer translates
y2 = similar(y1)
for k in 1:supp[2]-1
	fill!(y2, 0.0)
	y2[k*2^R-1+supp_index] = phi
	@show inner(x, y1, y2)
end


# ------------------------------------------------------------
# Orthogonality of internal and boundary scaling functions

x, Lphi = DaubScaling(F, R, Val{'L'})
for k in 1:N
	@show inner( x, phi, Lphi[:,k] )
	for l in 1:N
		@show inner( x, Lphi[:,l], Lphi[:,k] )
	end
end
=#

