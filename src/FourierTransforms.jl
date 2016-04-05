# ------------------------------------------------------------
# Exact Fourier transforms of Haar wavelet/scaling function

@doc """
	FourHaarScaling(xi[, J, k])

The Fourier transform of the Haar scaling function on scale `J` and translation `k` evaluated at `xi`.
If not supplied, `J = 0` and `k = 0`.
"""->
function FourHaarScaling{T<:Real}(xi::T)
	if xi == zero(xi)
		return one(Complex{Float64})
	else
		return (1.0 - exp(-2.0*pi*im*xi)) / (2.0*pi*im*xi)
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
		return (1.0 - exp(-pi*im*xi))^2 / (2.0*pi*im*xi)
	end
end


# ------------------------------------------------------------
# Fourier transform of Daubechies wavelets

@doc """
	feval(xi::Float64, C::Vector{Float64})

*F*ilter *eval*uation at `xi` of the filter `C`, i.e., compute
the linear combination 

	sum( C[n]*exp(2*pi*n*xi) )
"""->
function feval(xi::Float64, C::Vector{Float64})
	y = zero(Complex{Float64})
	for n = 1:length(C)
		y += C[n]*cis(2.0*pi*(n-1)*xi)
	end
	return y
end

@doc """
	FourDaubScaling(xi, C; ...)
	FourDaubScaling(xi, N; ...)

The Fourier transform of the Daubechies `N` scaling function evaluated
	at array `xi`. 
`N` is the number of zeros at -1.
For a single `xi`, the filter coefficients `C` must be supplied instead
of `N`.

The function is computed as an 'infinite' product;
to control this there are optional arguments:

- `prec`: Include factors that are numerically smaller than 1-prec.
- `maxCount`: The maximum number of factors.
"""->
function FourDaubScaling{T<:Real}( xi::T, C::Vector{Float64}; prec=eps(), maxCount=100)
	@assert isapprox(sum(C), 1.0)
	@assert prec >= eps()
	@assert maxCount >= 1

	const almost1 = 1.0 - prec
	Y = one(Complex{Float64})
	count = 1
	while count <= maxCount
		xi /= 2.0
		y = feval(xi, C)
		Y *= y

		# Convergence check: |y(xi) - 1| is small for small xi. But y is
		# exactly 1 in all even integers, so prevent premature exit
		abs(xi) >= one(xi) ? continue : count += 1

		if abs(y) >= almost1
			break
		end
	end

	return Y
end

function FourDaubScaling{T<:Real}( xi::DenseArray{T}, C::Vector{Float64}; args... )
	# TODO: eachindex
	Y = Array(Complex{Float64}, size(xi))
	for m = 1:length(xi)
		@inbounds Y[m] = FourDaubScaling( xi[m], C; args... )
	end

	return Y
end

function FourDaubScaling{T<:Real}( xi::DenseArray{T}, N::Int; args... )
	# Filter coefficients
	C = ifilter(N)
	scale!(C, 1/sum(C))

	return FourDaubScaling( xi, C; args... )
end


@doc """
	FourDaubWavelet(xi, N[, J, k]; ...)

The Fourier transform of the Daubechies `N` wavelet function evaluated at `xi`.
`N` is the number of zeros at -1.

The optional arguments are passed to `FourDaubScaling`.
"""->
function FourDaubWavelet{T<:Real}( xi::T, C::Vector{Float64}; args... )
	@assert isapprox(sum(C), 1.0)

	xi /= 2
	Y = FourDaubScaling(xi, C; args...)

	# High pass filter
	Y *= feval( -xi-0.5, C )
	Y *= cis( 2*pi*xi )

	return Y
end

function FourDaubWavelet{T<:Real}( xi::DenseArray{T}, N::Integer; args... )
	# Filter coefficients
	C = ifilter(N)
	scale!(C, 1/sum(C))

	Y = Array(Complex{Float64}, size(xi))
	for m = 1:length(xi)
		@inbounds Y[m] = FourDaubWavelet( xi[m], C; args... )
	end

	return Y
end

@doc """
	FourScalingFunc( xi, wavename::String, J=0, k=0; ... )

Compute the Fourier transform at `xi` of the scaling function `wavename` with
scale `J` and translation `k`.
"""->
function FourScalingFunc( xi, wavename::AbstractString, J::Integer=0, k::Integer=0; args... )
	@assert J >= 0 "Scale must be a non-negative integer"

	lowername = lowercase(wavename)
	if lowername == "haar" || lowername == "db1"
		return FourHaarScaling(xi, J, k)
	elseif isdaubechies(lowername)
		vm = van_moment(lowername)
		return FourDaubScaling(xi, vm, J, k)
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

	for i = 1:vm
		@inbounds UV[i, 1:vm+2*i-1] = bfilter(F, i-1) / sqrt2
	end

	return UV[:,1:vm], UV[:,vm+1:end]
end

@doc """
	FourDaubScaling( xi, F:::ScalingFilters; ... ) -> Matrix

The Fourier transform of the Daubechies `N` boundary wavelet transform
defined by the filters `F` evaluated at `xi`.
"""->
function FourDaubScaling( xi::Number, F::ScalingFilters; prec=sqrt(eps()), maxcount=50 )
	# TODO: Assertions in macro?
	@assert prec >= eps()
	@assert maxcount >= 1

	const U, V = UVmat(F.left)
	const vm = van_moment(F)

	# The notation is copied directly from the article of Poon & Gataric (see references in docs).
	# v1(xi) is the vector of desired Fourier transforms and v10 = v1(0)
	const Vcol = size(V,2)
	const v20 = ones(Float64, Vcol)
	const v10 = (eye(U) - U) \ V*v20

	const C = F.internal / sum(F.internal)
	const v2index = vm + [0:(Vcol-1);]
	v2(xi) = FourDaubScaling(xi, C)*cis(-2*pi*v2index*xi)

	# For large j,
	# v1(xi) \approx U^j*v10 + sum_{l=0}^{j-1} U^l*V*v2(xi/2^{l+1}) 
	v1 = U^maxcount*v10
	for l = 0:maxcount-1
		v1 += U^l*V*v2(xi/2^(l+1))
	end

	return v1
end

function FourDaubScaling( xi::AbstractVector, F::ScalingFilters; args... )
	Nxi = length(xi)
	Y = Array(Complex{Float64}, van_moment(F), Nxi)

	for nxi = 1:Nxi
		@inbounds Y[:,nxi] = FourDaubScaling(xi[nxi], F; args...)
	end

	return Y'
end

# ------------------------------------------------------------
# Dilation and translation. 

for name in [:FourHaarScaling, :FourHaarWavelet]
	@eval begin
		@vectorize_1arg Real $name

		function $name{T<:Real}(xi::DenseArray{T}, J::Int)
			M = length(xi)
			C2 = 2.0^(-J)
			C = 2.0^(-J/2)
			y = Array(Complex{Float64}, M)

			for m = 1:M
				@inbounds y[m] = C*$name( C2*xi[m] )
			end
			return y
		end

		function $name{T<:Real}(xi::DenseArray{T}, J::Int, k::Int)
			y = $name(xi, J)
			D = exp( -2.0*pi*im*2.0^(-J)*k*xi )
			had!(y, D)
			return y
		end
	end
end

for name in [:FourDaubScaling, :FourDaubWavelet]
	@eval begin
		function $name{T<:Real}(xi::DenseArray{T}, N::Integer, J::Integer; args...)
			scale_xi = scale( 2.0^(-J), xi )
			y = $name( scale_xi, N; args... )
			scale!( 2.0^(-J/2), y )
			return y
		end

		function $name{T<:Real}(xi::DenseArray{T}, N::Integer, J::Int, k::Integer; args...)
			y = $name(xi, N, J; args... )
			D = exp( -2.0*pi*im*2.0^(-J)*k*xi )
			had!(y, D)
			return y
		end
	end
end

