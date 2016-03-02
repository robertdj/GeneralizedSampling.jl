# ------------------------------------------------------------
# Scaling functions

@doc """
	HaarScaling(xi[, J, k])

The Haar scaling function evaluated in `xi` at level `J` and translation `k`.
By default, `J=0` and `k=0`.
"""->
function HaarScaling{T<:Real}(xi::T)
	0 <= xi < 1 ? 1 : 0
end

@doc """
	DaubScaling(N, L[, J, k]) -> x,y

A Daubechies `N` scaling function at scale `J` and translation `k` evaluated in `0:2^(-L):2N-1`.
"""->
function DaubScaling(N::Int, L::Int)
	wt = wavelet( WT.Daubechies{N}() )
	phi, x = DaubScaling( wt.qmf, L )
end

@doc """
	DaubScaling(C::Vector) -> y

Compute function values of the scaling function defined by the filter `C` at the integers in the support.
"""->
function DaubScaling(C::Vector{Float64})
	L = dyadic_dil_matrix(C)

	# Eigenvector of eigenvalue 1 (the largest)
	eigenvec = eigs(L; nev=1)[2]
	E = vec(real(eigenvec))

	# eigs is not consistent with the sign
	if E[2] < 0
		scale!(E, -1.0)
	end

	# The first and last entry are both 0
	# TODO: Don't compute them
	E[1] = E[end] = 0.0

	return E
end

@doc """
	DaubScaling(C::Vector, R::Int) -> x, y

Compute function values of the scaling function defined by the filter
`C` at the dyadic rationals of resolution `R` in the support.
"""->
function DaubScaling(C::Vector{Float64}, R::Int)
	max_supp = support(C)
	x = dyadic_rationals(max_supp, R)
	Nx = length(x)
	phi = zeros(Float64, Nx)

	# Base level
	cur_idx = dyadic_rationals(max_supp, R, 0)
	phi[cur_idx] = DaubScaling(C)

	# Recursion: Fill remaining levels
	NC = length(C)
	sqrt2 = sqrt(2)
	for L = 1:R
		cur_idx = dyadic_rationals(max_supp, R, L)

		for phin = cur_idx
			for Ck = 1:NC
				pidx = dyadic_parent(phin, Ck-1, R)
				# TODO: Calculate only the necessary pidx
				if 1 <= pidx <= Nx
					phi[phin] += sqrt2*C[Ck]*phi[pidx]
				end
			end
		end
	end

	return x, phi
end

@doc """
	DaubScaling(N::Int, edge::Char)

Compute function values of all the boundary scaling function with `N` vanishing moments at the integers in its support.

`edge` can be either `L`(eft) og `R`(ight).

The output is
"""->
function DaubScaling(N::Int, edge::Char)
	# TODO: Check input
end

@doc """
	boundary_coef_mat(F::BoundaryFilter) -> Matrix

The boundary coefficients collected in a matrix with the `i`'th row
containing the coefficients of the `i`'th boundary scaling function.
"""->
function boundary_coef_mat(F::BoundaryFilter)
	const vm = van_moment(F)
	coef_mat = zeros(Float64, vm, vm)

	const sqrt2 = sqrt(2)
	for row = 1:vm
		coef = bfilter(F, row-1)
		coef_mat[row,:] = sqrt2*coef[1:vm]
	end

	return coef_mat
end

#=
@doc """
	DaubScaling(F::BoundaryFilter, boundary) -> Matrix 

Compute function values of the boundary scaling function defined by the filters `F` at the integers in the support.

`boundary` must be either `Val{'L'}` or `Val{'R'}`.

The output is a Matrix where the `(i,j)`'th entry is the `j`'th boundary function evaluated at `i-1`.
"""->
=#
@debug function DaubScaling(F::BoundaryFilter, ::Type{Val{'L'}})
	# The number of vanishing moments = the number of boundary functions
	const vm = van_moment(F)

	# Bounds of support of boundary and internal scaling functions
	const boundary_support = (0, 2*vm-1)
	const internal_support = (-vm+1, vm)

	# Convert between function values and indices for boundary and
	# internal scaling functions
	boundary_x2index(x::Int) = x + 1
	boundary_index2x(idx::Int) = idx - 1
	internal_x2index(x::Int) = x + vm

	# Function values of the internal scaling function
	C = wavelet( WT.Daubechies{vm}() ).qmf
	internal = DaubScaling(C)

	# The number of integers to evaluate
	const Nint = internal_support[2] - internal_support[1] + 1
	Y = zeros(Float64, Nint, vm)

	# Boundary scaling functions at 0
	# TODO: Move to separate function?
	boundary_mat = boundary_coef_mat(F)
	# TODO: Are we sure that the eigvector corresponding to 1 is always the first?
	eigenvecs = eigvecs(boundary_mat)
	Y[1,:] = eigenvecs[:,1]

	# Remaining boundary scaling functions values from highest to lowest to use the recursion
	# The order of the loops are intentional: All y values are needed for each x value
	const sqrt2 = sqrt(2)
	for i = (Nint-1):-1:2
		xval = boundary_index2x(i)
		xval2 = 2*xval
		xval2index = boundary_x2index(xval2)

		for j = 1:vm
			cur_filter = bfilter(F, j-1)

			# Boundary contribution
			if boundary_support[1] <= xval2 < boundary_support[2]
				for l = 0:vm-1
					Y[i,j] += sqrt2*cur_filter[l+1]*Y[xval2index,l+1]
				end
			end

			# Internal contribution
			for m = vm:(vm+2*(j-1))
				if internal_support[1] < xval2-m < internal_support[2]
					Y[i,j] += sqrt2*cur_filter[m+1]*internal[ internal_x2index(xval2-m) ]
				end
			end
		end
	end

	return Y
end

function DaubScaling(C::Vector{Float64}, R::Int, ::Val{'L'})
end


@doc """
	dyadic_rationals(upper, res)

The dyadic rationals of resolution `R` in the integer interval `[0, upper]`.
"""->
function dyadic_rationals(upper::Int, res::Int)
	collect( 0:2.0^(-res):upper )
end

@doc """
	dyadic_rationals(upper, res, level)

In a vector of dyadic rationals up to resolution `res` in `[0,upper]`,
return the indices of those at exactly `level`.
"""->
function dyadic_rationals(upper::Int, res::Int, level::Int)
	N = upper*2^res + 1

	if level == 0
		step = 2^res
		return 1:step:N
	else
		step = 2^(res-level)
		return 1+step:2*step:N
	end
end

@doc """
	dyadic_rationals(dya_rat::Vector, level::Int)

In a vector of dyadic rationals return the indices of those up to `level`.
"""->
function dyadic_rationals(dy_rat::Vector{Float64}, level::Int)
	# Check input
	res = -log2(dy_rat[2])
	@assert isuniform(dy_rat) && isinteger( res ) "Input is not a vector of dyadic rationals"
	@assert 0 <= level <= res "Resolution must be greater than level"

	power2 = 2^level
	max_supp = dy_rat[end]
	Nlevel = Int(max_supp*power2+1)
	dyadic_level = Array(Int, Nlevel)

	# An entry in dy_rat (at level res) is also in dyadic_level if
	# and only if it is an integer when multiplied with power2
	count = 1
	Nres = length(dy_rat)
	for n = 1:Nres
		if isinteger( power2*dy_rat[n] )
			dyadic_level[count] = n
			count += 1
		end
	end

	return dyadic_level
end

@doc """
	dyadic_dil_matrix(C::Vector)

The "dyadic dilation matrix" of the filter `C` with `(i,j)`'th entry
`C[2i-j]`.
"""->
function dyadic_dil_matrix(C::Vector{Float64})
	NC = length(C)
	dydil_mat = zeros(Float64, NC, NC)

	sqrt2 = sqrt(2)
	for nj = 1:NC, ni = 1:NC
		Cidx = 2*ni - nj
		# TODO: Avoid this check?
		if 1 <= Cidx <= NC
			dydil_mat[ni, nj] = sqrt2*C[Cidx]
		end
	end

	return dydil_mat
end

@doc """
	support(C) -> length

Return the length of the support of the scaling function defined by the
filter vector `C`.
"""->
function support(C::Vector{Float64})
	length(C) - 1
end

@doc """
	support(phi, J, k) -> lower, upper

From a vector of scaling function values `phi`, return the `lower` and
`upper` bound of the support of the version that is dilated with `J` and
translated with `k`.
"""->
function support(phi::Vector{Float64}, J::Int, k::Int)
	upper, res = factor_support( length(phi) )

	lower = k*2^(res-J) + 1
	upper = (upper+k)*2^(res-J) + 1

	return lower, upper
end

@doc """
	factor_support(L::Int) -> upper, resolution

From the length of a vector of scaling function evaluations from
`DaubScaling`, return the `upper` endpoint of the scaling function's
support and the `resolution` of the the dyadic rationals.
"""->
function factor_support(L::Int)
	# L = upper*2^R + 1, where upper is uneven
	R = 0
	L -= 1
	while iseven(L)
		R += 1
		L = div(L, 2)
	end

	return L, R
end

@doc """
	dyadic_parent(i::Int, k::Int, L::Int)

In the vector `x` where the `i`'th entry is `x_i = (i-1)/2^L`,
`dyadic_parent` returns the index of `2x_i - k`.
"""->
function dyadic_parent(i::Int, k::Int, L::Int)
	2*i - 1 - k*2^L
end


# ------------------------------------------------------------
# Dilation and translation

function HaarScaling{T<:Real}(xi::T, J::Int)
	2.0^(J/2)*HaarScaling(2^J*xi)
end

function HaarScaling{T<:Real}(xi::T, J::Int, k::Int)
	@assert 0 <= k < (dil_fact = 2^J) "Translation is too large"
	2.0^(J/2)*HaarScaling(dil_fact*xi - k)
end

@doc """
	DaubScaling(phi::Vector{Float64}, J::Int, k::Int)

From the base scaling function `phi`, compute `phi_jk` of dilation `J`
and translation `k`.
"""->
function DaubScaling(phi::Vector{Float64}, J::Int, k::Int)
	Nphi = length(phi)
	upper, resolution = factor_support( Nphi )

	@assert 0 <= J <= resolution 
	@assert 0 <= k < upper*2^J

	# If x_{i,L} = (i-1)/2^L, the dilated & translated scaling function
	# phi_{J,k}(x_{i,L}) = 2^(j/2)*phi( x_{i-k*2^{L-J}, L-J} )
	# TODO: Pass this as an argument?
	dy_ras = dyadic_rationals(upper, resolution)
	xidx = dyadic_rationals(dy_ras, resolution-J)

	lower, upper = support(phi, J, k)
	supp = [lower:upper;]
	@assert (Nsupp = length(supp)) == length(xidx)

	phijk = zeros(Float64, Nphi)
	norm_fact = 2.0^(J/2)
	for nsupp = 1:Nsupp
		phijk[ supp[nsupp] ] = norm_fact*phi[ xidx[nsupp] ]
	end

	return phijk
end

# ------------------------------------------------------------
# Vectorization

@vectorize_1arg Real HaarScaling

function HaarScaling{T<:Real}(xi::AbstractArray{T}, J::Int, k::Int)
	N = length(xi)
	y = Array(Float64, size(xi))
	for n = 1:N
		@inbounds y[n] = HaarScaling( xi[n], J, k )
	end
	return y
end


# ------------------------------------------------------------

# Reconstruction in Haar basis
function weval(coef::AbstractArray{Float64}, res::Int)
	x = dyadic_rationals(1, res)
	Nx = length(x)
	Ncoef = length(coef)
	# Scale of wavelet transform.
	# TODO: Error check
	J = Int( log2(Ncoef) )

	y = zeros(Float64, Nx)
	for n = 1:Nx
		for m = 1:Ncoef
			# TODO: Only include the functions that have x[nx] in their support
			@inbounds y[n] += coef[m]*HaarScaling( x[n], J, m-1 )
		end
	end

	return x, y
end

@doc """
	weval(coef::Vector, N::Int, J::Int, res::Int)

Evaluate `coef` vector in the Daubechies `N` basis at scale `J` in the dyadic rationals of resolution `res`.

`J` must be so large that the support of the wavelet is contained in `[0,1]`.
"""->
function weval(coef::AbstractArray{Float64}, N::Int, J::Int, res::Int)
	# TODO: J must be sufficiently large 

	x = dyadic_rationals(N, res)
	Nx = length(x)
	Ncoef = length(coef)
	# TODO: Compute J from coef?

	y = zeros(Float64, Nx)
	for n = 1:Nx
		for m = 1:Ncoef
			# TODO: Only include the functions that have x[nx] in their support
			@inbounds y[n] += coef[m]*HaarScaling( x[n], J, m-1 )
		end
	end

	return x, y
end

