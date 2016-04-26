# ------------------------------------------------------------
# Exact Fourier transforms of Haar scaling function

@doc """
	FourHaarScaling(xi[, J, k])

The Fourier transform of the Haar scaling function on scale `J` and translation `k` evaluated at `xi`.
If not supplied, `J = 0` and `k = 0`.
"""->
function FourHaarScaling(xi::Real)
	if xi == zero(xi)
		return ComplexOne
	else
		return (1.0 - cis(-2.0*pi*xi)) / (2.0*pi*im*xi)
	end
end


# ------------------------------------------------------------
# Fourier transform of Daubechies scaling functions

@doc """
	feval(xi::Float64, C::Vector{Float64})

*F*ilter *eval*uation at `xi` of the filter `C`, i.e., compute

	sum( C[n]*exp(-2*pi*n*xi) )
"""->
function feval(xi::Float64, C::Vector{Float64})
	y = zero(Complex{Float64})
	for n in 1:length(C)
		@inbounds y += C[n]*cis(-2.0*pi*(n-1)*xi)
	end
	return y
end

@doc """
	FourDaubScaling(xi, C; ...)

The Fourier transform of the Daubechies scaling function with filter `C` evaluated at `xi`. 

The function is computed as an 'infinite' product;
to control this there are optional arguments:

- `prec`: Include factors that are numerically smaller than 1-prec.
- `maxcount`: The maximum number of factors.
"""->
function FourDaubScaling( xi::Real, C::Vector{Float64}; prec=SMALL_PREC, maxcount=100)
	# TODO: Move these checks to vector version?
	# TODO: Assertions in macro?
	@assert isapprox(sum(C), 1.0)
	@assert prec >= SMALL_PREC
	@assert maxcount >= 1

	const almost1 = 1.0 - prec
	const half = 0.5*one(xi)
	Y = ComplexOne
	count = 1
	while count <= maxcount
		xi *= 0.5
		y = feval(xi, C)
		Y *= y

		# Convergence check: |y(xi) - 1| is small for small xi. But y is
		# exactly 1 in all even integers, so prevent premature exit
		if abs(xi) <= half
			abs(y) >= almost1 && break
			count += 1
		end
	end

	return Y
end

function FourDaubScaling( xi::AbstractArray{Float64}, C::Vector{Float64}; args... )
	Y = Array(Complex{Float64}, size(xi))
	for idx in eachindex(xi)
		@inbounds Y[idx] = FourDaubScaling( xi[idx], C; args... )
	end

	return Y
end

# TODO: method w/o J and k?
function FourDaubScaling( xi, N::Int, J::Integer=0, k::Integer=0; args... )
	C = ifilter(N)
	scale!(C, 1/sum(C))
	FourDaubScaling( xi, C, J, k; args... )
end


@doc """
	FourScalingFunc( xi, wavename, J=0, k=0; ... )

Compute the Fourier transform at `xi` of the scaling function `wavename` with
scale `J` and translation `k`.

Optional arguments are passed to the Fourier transform of `wavename` (if relevant).
"""->
function FourScalingFunc( xi, wavename::AbstractString, J::Integer=0, k::Integer=0; args... )
	@assert J >= 0 "Scale must be a non-negative integer"

	if ishaar(wavename)
		return FourHaarScaling(xi, J, k)
	elseif isdaubechies(wavename)
		return FourDaubScaling(xi, van_moment(wavename), J, k; args...)
	else
		error("Fourier transform for this wavelet is not implemented")
	end

end


# ------------------------------------------------------------
# Fourier transform of Daubechies boundary wavelets

@doc """
	UVmat(F::BoundaryFilter) -> Matrix, Matrix

Return the matrices `U` and `V` used in `FourDaubScaling` for boundary scaling functions.
"""->
function UVmat(B::BoundaryFilter)
	const vm = length(B)

	UV = zeros(Complex{Float64}, vm, 3*vm-1)

	for i in 1:vm
		@inbounds UV[i, 1:vm+2*i-1] = bfilter(B, i-1) / sqrt2
	end

	return UV[:,1:vm], UV[:,vm+1:end]
end

@doc """
	FourDaubScaling( xi, N[, J], side::Char; args... ) -> Matrix

The Fourier transform at `xi` of the Daubechies boundary scaling function with `N` vanishing moments at scale `J` at the left or right `side` ('L' or 'R').

The output is an `N`-by-`length(xi)` matrix where column `j` holds the Fourier transform of the `N` boundary functions at `xi[j]`.
"""->
function FourDaubScaling{T<:Real}( xi::AbstractVector{T}, N::Integer, side::Char; maxcount=40 )
	@assert 1 <= maxcount <= 100
	#=
	- F is a vector with the Fourier transforms of the boundary scaling
	function at the current xi
	- F0 is F at the zero vector
	- phihat is the Fourier transform of the internal scaling function
	at the current xi (with different phaseshifts)

	For large j,
	F(xi) \approx U^j*F0 + sum_{l=0}^{j-1} U^l*V*phihat( xi/2^{l+1} )
	where U & V are from UVmat
	=#

	const U, V = UVmat(bfilter(N, side))
	const C = ifilter(N) / sum(ifilter(N))

	const F0 = complex( (eye(U) - U) \ vec(sum(V,2)) )

	# ----------------------------------------
	# Precompute U^l and U^l*V

	Upower = Dict{Int, Matrix{Complex{Float64}}}()
	sizehint!(Upower, maxcount)
	Upower[0] = eye(U)

	UpowerV = Dict{Int, Matrix{Complex{Float64}}}()
	sizehint!(UpowerV, maxcount)
	UpowerV[0] = V

	for l in 1:maxcount
		Upower[l] = U*Upower[l-1]
		UpowerV[l] = U*UpowerV[l-1]
	end

	# ----------------------------------------

	const Nxi = length(xi)
	Y = Array(Complex{Float64}, N, Nxi)

	# Arrays with temporary results
	F = Array(Complex{Float64}, N)
	const Vcol = size(V,2) - 1
	phihat = Array(Complex{Float64}, Vcol+1)

	for nxi in 1:Nxi
		fill!(F, zero(Complex{Float64}))

		current_xi = xi[nxi]
		for l in 0:maxcount-1
			current_xi *= 0.5
			internal = FourDaubScaling(current_xi, C)
			for k in 0:Vcol
				# TODO: 2pi as a const?
				phihat[k+1] = internal*cis( -2.0*pi*(N+k)*current_xi )
			end

			# F = F + U^l*V*z
			BLAS.gemv!('N', ComplexOne, UpowerV[l], phihat, ComplexOne, F)

			# TODO: Test for convergence after each pass through the loop?
		end
		# F = F + U^maxcount*F0
		BLAS.gemv!('N', ComplexOne, Upower[maxcount], F0, ComplexOne, F)

		Y[:,nxi] = F
	end

	return Y
end


# ------------------------------------------------------------
# Dilation and translation. 

@vectorize_1arg Real FourHaarScaling

function FourHaarScaling(xi::Real, J::Integer)
	sqrt(2.0^(-J)) * FourHaarScaling(2.0^(-J)*xi)
end

function FourHaarScaling(xi::Real, J::Integer, k::Integer)
	cis( -2.0*pi*2.0^(-J)*k*xi ) * FourHaarScaling(xi, J)
end

function FourHaarScaling{T<:Real}( xi::AbstractArray{T}, J::Integer)
	# TODO: Save 2^-J and do above calculations explicitly?
	Y = Array(Complex{Float64}, size(xi))
	for idx in eachindex(xi)
		@inbounds Y[idx] = FourHaarScaling( xi[idx], J )
	end

	return Y
end

function FourHaarScaling{T<:Real}( xi::AbstractArray{T}, J::Integer, k::Integer)
	Y = Array(Complex{Float64}, size(xi))
	for idx in eachindex(xi)
		@inbounds Y[idx] = FourHaarScaling( xi[idx], J, k )
	end

	return Y
end


function FourDaubScaling(xi::Real, C::Vector{Float64}, J::Integer; args...)
	sqrt(2.0^(-J)) * FourDaubScaling(2.0^(-J)*xi, C; args...)
end

function FourDaubScaling(xi::Real, C::Vector{Float64}, J::Integer, k::Integer; args...)
	cis( -2.0*pi*2.0^(-J)*k*xi ) * FourDaubScaling(xi, C, J; args...)
end

function FourDaubScaling{T<:Real}( xi::AbstractArray{T}, C::Vector{Float64}, J::Integer; args... )
	# TODO: Save 2^-J and do above calculations explicitly?
	Y = Array(Complex{Float64}, size(xi))
	for idx in eachindex(xi)
		@inbounds Y[idx] = FourDaubScaling( xi[idx], C, J; args... )
	end

	return Y
end

function FourDaubScaling{T<:Real}( xi::AbstractArray{T}, C::Vector{Float64}, J::Integer, k::Integer; args... )
	Y = Array(Complex{Float64}, size(xi))
	for idx in eachindex(xi)
		@inbounds Y[idx] = FourDaubScaling( xi[idx], C, J, k; args... )
	end

	return Y
end


function FourDaubScaling( xi, N::Integer, J::Integer, side::Char; args... )
	@assert J >= 0

	xi *= 2.0^(-J)
	y = FourDaubScaling(xi, N, side; args...)
	scale!( y, sqrt(2.0^(-J)) )

	return y
end

