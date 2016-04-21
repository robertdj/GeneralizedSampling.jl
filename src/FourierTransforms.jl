# ------------------------------------------------------------
# Exact Fourier transforms of Haar wavelet/scaling function

@doc """
	FourHaarScaling(xi[, J, k])

The Fourier transform of the Haar scaling function on scale `J` and translation `k` evaluated at `xi`.
If not supplied, `J = 0` and `k = 0`.
"""->
function FourHaarScaling{T<:Real}(xi::T)
	if xi == zero(xi)
		return ComplexOne
	else
		return (1.0 - cis(-2.0*pi*xi)) / (2.0*pi*im*xi)
	end
end

@doc """
	FourHaarWavelet(xi[, J, k])

The Fourier transform of the Haar wavelet on scale `J` and translation `k` evaluated at `xi`.
If not supplied, `J = 0` and `k = 0`.
"""->
function FourHaarWavelet{T<:Real}(xi::T)
	if xi == zero(xi)
		return zero(Complex{Float64})
	else
		return (1.0 - cis(-pi*xi))^2 / (2.0*pi*im*xi)
	end
end


# ------------------------------------------------------------
# Fourier transform of Daubechies wavelets

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
	# TODO: Make this check optional?
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

function FourDaubScaling{T<:Real}( xi::AbstractArray{T}, C::Vector{Float64}; args... )
	Y = Array(Complex{Float64}, size(xi))
	for idx in eachindex(xi)
		@inbounds Y[idx] = FourDaubScaling( xi[idx], C; args... )
	end

	return Y
end

# TODO: method w/o J and k?
function FourDaubScaling( xi, N::Int, J::Integer, k::Integer; args... )
	C = ifilter(N)
	scale!(C, 1/sum(C))
	FourDaubScaling( xi, C, J, k; args... )
end

@doc """
	FourDaubWavelet(xi, N[, J, k]; ...)

The Fourier transform of the Daubechies `N` wavelet function evaluated at `xi`.
`N` is the number of zeros at -1.

The optional arguments are passed to `FourDaubScaling`.
"""->
function FourDaubWavelet{T<:Real}( xi::T, C::Vector{Float64}; args... )
	@assert isapprox(sum(C), 1.0)

	xi /= 2.0
	Y = FourDaubScaling(xi, C; args...)

	# High pass filter
	Y *= feval( -xi-0.5, C )
	Y *= cis( 2*pi*xi )

	return Y
end

function FourDaubWavelet{T<:Real}( xi::AbstractArray{T}, N::Integer, J::Integer=0, k::Integer=0; args... )
	# Filter coefficients
	C = ifilter(N)
	scale!(C, 1/sum(C))

	Y = Array(Complex{Float64}, size(xi))
	for idx in eachindex(xi)
		@inbounds Y[idx] = FourDaubWavelet( xi[idx], C, J, k; args... )
	end

	return Y
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

Return the matrices `U` and `V` used in `FourDaubScaling`.
This notation is copied directly from the article of Poon & Gataric (see references in docs).

Only intended for internal use!
"""->
function UVmat(F::BoundaryFilter)
	const vm = van_moment(F)

	UV = zeros(Float64, vm, 3*vm-1)

	for i in 1:vm
		@inbounds UV[i, 1:vm+2*i-1] = bfilter(F, i-1) / sqrt2
	end

	return UV[:,1:vm], UV[:,vm+1:end]
end


@doc """
	FourDaubScaling( xi, N, side::Char; args... )

The Fourier transform at `xi` of the Daubechies boundary scaling function with `N` vanishing moments at the left or right `side` ('L' or 'R').
"""->
function FourDaubScaling( xi, N::Integer, side::Char; args... )
	const B = bfilter(N, side)
	const U, V = UVmat(B)
	C = ifilter(N)
	scale!(C, 1/sum(C))

	FourDaubScaling( xi, C, U, V; args... )
end

@doc """
	FourDaubScaling( xi::Real, C::Vector, U::Matrix, V::Matrix; ... ) -> Vector

The Fourier transform at `xi` of all the Daubechies boundary scaling functions. 
`C` is the vector of coefficients for the internal scaling function.
`U` and `V` are the matrices calculated from the boundary filter.

The optional arguments are the precision `prec` and maximum number of iterations `maxcount`.
"""->
function FourDaubScaling( xi::Real, C::Vector{Float64}, U::Matrix{Float64}, V::Matrix{Float64}; prec=LARGE_PREC, maxcount=50 )
	# TODO: Assertions in macro?
	@assert prec >= SMALL_PREC
	@assert maxcount >= 1
	# TODO: Skip this assertion?
	@assert isapprox(sum(C), 1.0)

	# The notation is copied directly from the article of Poon & Gataric (see references in docs).
	# v1(xi) is the vector of desired Fourier transforms and v10 = v1(0)
	# TODO: These lines are independent of xi. Move to array method?
	const Vcol = size(V,2)
	const v20 = ones(Float64, Vcol)
	const v10 = (eye(U) - U) \ V*v20

	#= const vm = size(U,1) =#
	const v2index = size(U,1) + [0:(Vcol-1);]
	v2(xi) = FourDaubScaling(xi, C)*cis(-2*pi*v2index*xi)

	# For large j,
	# v1(xi) \approx U^j*v10 + sum_{l=0}^{j-1} U^l*V*v2(xi/2^{l+1}) 
	v1 = U^maxcount*v10
	for l in 0:maxcount-1
		v1 += U^l*V*v2(xi/2^(l+1))
	end

	return v1
end

function FourDaubScaling{T<:Real}( xi::AbstractVector{T}, C::Vector{Float64}, U::Matrix{Float64}, V::Matrix{Float64}; args... )
	const Nxi = length(xi)
	Y = Array(Complex{Float64}, size(U,1), Nxi)

	for nxi in 1:Nxi
		@inbounds Y[:,nxi] = FourDaubScaling( xi[nxi], C, U, V; args... )
	end

	return Y'
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

