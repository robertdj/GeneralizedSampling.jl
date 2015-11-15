# ------------------------------------------------------------
# Scaling function

function HaarScaling{T<:Real}(xi::T)
	0 <= xi < 1 ? 1 : 0
end

# ------------------------------------------------------------
# Dilation, translation and vectorization

function HaarScaling{T<:Real}(xi::T, J::Int)
	2.0^(J/2)*HaarScaling(2.0^J*xi)
end

function HaarScaling{T<:Real}(xi::T, J::Int, k::Int)
	2.0^(J/2)*HaarScaling(2.0^J*xi - k)
end

function HaarScaling{T<:Real}(xi::AbstractArray{T}, J::Int, k::Int)
	map( xi -> HaarScaling(xi,J,k), xi )
end


# ------------------------------------------------------------

@doc """
	weval(coef, J::Int, L::Int)

Evaluate `coef` vector in the Haar basis at scale `J` in the points
`0:2^-L:1-2^-L`.
"""->
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
			y[n] += coef[m]*HaarScaling( x[n], J, m-1 )
		end
	end

	return x, y
end

