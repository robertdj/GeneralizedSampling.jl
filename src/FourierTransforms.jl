# Fourier transform of Daubechies scaling function

@doc """
	FourDaubScaling(xi, N)

The Fourier transform of the Daubechies `N` scaling function evaluated at `xi`.
`N` is the number of zeros at -1.
"""->
function FourDaubScaling( xi::Real, N::Int; prec=sqrt(eps()) )
	xi /= 2
	y = DaubLowPass(xi, Val{N})
	Y = copy(y)

	n = 1
	while abs(y) <= 1-prec && n <= 30
		xi /= 2
		y = DaubLowPass(xi, Val{N})
		Y *= y

		n += 1
	end

	return Y
end

# Vectorized version
function FourDaubScaling{T<:Real}(xi::Array{T}, N::Int)
	map( xi -> FourDaubScaling(xi, N), xi )
end


function DaubLowPass(xi::Real, ::Type{Val{0}})
	return 0.5*(1 + exp(im*xi))
end

