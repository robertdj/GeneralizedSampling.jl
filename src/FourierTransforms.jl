# Fourier transform of Daubechies scaling function

@doc """
### FourDaubScaling(xi, N)

The Fourier transform of the Daubechies `N` scaling function evaluated at `xi`.
`N` is the number of zeros at -1.

The function is computed as an 'infinite' product;
to control this there are optional arguments:

- `prec`: Include factors that are numerically smaller than 1-prec.
- `M`: The maximum number of factors.
"""->
function FourDaubScaling( xi::Real, N::Int; prec=eps(), M::Int=20 )
	xi /= 2.0
	Y = y = DaubLowPass(xi, Val{N})

	m = 1
	while abs(y) <= 1-prec && m <= M
		xi /= 2.0
		y = DaubLowPass(xi, Val{N})
		Y *= y

		m += 1
	end

	return Y
end

# Vectorized version
function FourDaubScaling{T<:Real}( xi::Array{T}, N::Int; arg... )
	map( xi -> FourDaubScaling(xi, N; arg...), xi )
end


function DaubLowPass(xi::Real, ::Type{Val{0}})
	# For Fourier transform with exp(-i*xi*x)
	return 0.5*(1.0 + exp(im*xi))
end

