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
	DaubScaling(N, L[, J, k])

A Daubechies `N` scaling function at scale `J` and translation `k` evaluated in `0:2^(-L):2N-1`.
"""->
function DaubScaling(N::Int, L::Int)
	wt = wavelet( WT.Daubechies{N}() )
	phi, x = DaubScaling( wt.qmf, L )
end

function DaubScaling(C::Vector{Float64}, R::Int)
	max_supp = support(C)*2^R
	x = collect(0:max_supp) / 2^R
	Nx = max_supp + 1
	phi = zeros(Float64, Nx)

	# Base level
	cur_idx = level2idx(0, R, Nx)
	phi[cur_idx] = DaubScaling(C)

	# Recursion: Fill remaining levels
	NC = length(C)
	sqrt2 = sqrt(2)
	for L = 1:R
		cur_idx = level2idx(L, R, Nx)

		for phin = cur_idx
			for Ck = 1:NC
				pidx = dyadic_parent(phin, Ck-1, R)
				if 1 <= pidx <= Nx
					phi[phin] += sqrt2*C[Ck]*phi[pidx]
				end
			end
		end
	end

	return x, phi
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
	DaubScaling(C::Vector)

Compute function values of the scaling function defined by the filter
`C` at the integers in the support.
"""->
function DaubScaling(C::Vector{Float64})
	L = dyadic_dil_matrix(C)

	# Eigenvector of eigenvalue 1 (the largest)
	# TODO: The first and last entry are both 0
	# TODO: eigs is not consistent with sign
	eigenvec = eigs(L; nev=1)[2]
	return real(vec(eigenvec))
end

@doc """
	level2idx(L::Int, R::Int, N::Int) -> idx

Return the indices at level `L` in a scaling function vector of length `N` at resolution `R`.
"""->
function level2idx(L::Int, R::Int, N::Int)
	@assert 0 <= L <= R
	@assert N >= 2^R

	if L == 0
		return collect( 1:2^R:N )
	else
		return 1 + collect( 2^(R-L):2^(R-L+1):N )
	end
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
	dyadic_parent(i::Int, k::Int, L::Int)

In the vector `x` where the `i`'th entry is `x_i = (i-1)/2^L`,
`dyadic_parent` returns the index of `2x_i - k`.
"""->
function dyadic_parent(i::Int, k::Int, L::Int)
	2*i - 1 - k*2^L
end


# ------------------------------------------------------------
# Dilation, translation

function HaarScaling{T<:Real}(xi::T, J::Int)
	2.0^(J/2)*HaarScaling(2^J*xi)
end

function HaarScaling{T<:Real}(xi::T, J::Int, k::Int)
	# TODO: Is this error check slowing things down?
	@assert 0 <= k < 2^J
	2.0^(J/2)*HaarScaling(2^J*xi - k)
end

function DaubScaling(C::Vector{Float64}, L::Int, J::Int, k::Int=0)
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

@doc """
	weval(coef, wave, L::Int)

Evaluate `coef` vector in the `wave` basis in the points `0:2^-L:1-2^-L`.
"""->
function weval(coef::AbstractArray{Float64}, J::Int, wave::AbstractString, L::Int)
end

# Reconstruction in Haar basis
function weval(coef::AbstractArray{Float64}, L::Int)
	inc = 2.0^(-L)
	x = collect( 0:inc:1-inc )
	Nx = length(x)
	Ncoef = length(coef)
	# Scale of wavelet transform.
	# TODO: Error check
	J = Int( log2(Ncoef) )

	y = zeros(Float64, Nx)
	for n = 1:Nx
		for m = 1:Ncoef
			# TODO: Only include the functions that have x[nx] in their
			# support
			@inbounds y[n] += coef[m]*HaarScaling( x[n], J, m-1 )
		end
	end

	return x, y
end

# Reconstruction in general Daubechies basis
function weval(coef::AbstractArray{Float64}, J::Int, N::Int, L::Int)
	# TODO: J must be so large that the support of the wavelet is
	# contained in [0,1]
end

