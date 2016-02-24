# ------------------------------------------------------------
# Boundary low pass filter coefficients for Daubechies wavelets
# http://www.pacm.princeton.edu/~ingrid/publications/54.txt

typealias BoundaryFilter Array{Any,1}

@doc """
	lfilter(N::Int) -> BoundaryFilter

Return the left boundary filters for the scaling functions with `N` vanishing moments.
"""->
function lfilter(N::Int)
	@assert 2 <= N <= 3 "Boundary filters are not available with this number of vanishing moments"
	return LEFT_FILTERS[N]
end

@doc """
	lfilter(BoundaryFilter, k::Int) -> Vector

Return the left boundary filter for the `k`'th scaling function (0 <= `k` < the number of vanishing moments).
"""->
function lfilter(F::BoundaryFilter, k::Int)
	@assert 0 <= k < length(F)
	return F[k+1]
end

# TODO: Obsolete?
@doc """
	lfilter(N::Int, k::Int) -> Vector

Return the left boundary filter for the `k`'th scaling function with `N` vanishing moments.
By construction, `k = 0,...,N-1`.
"""->
function lfilter(N::Int, k::Int)
	@assert 0 <= k < N
	return lfilter(N)[k+1]
end

@doc """
	van_moment(F::BoundaryFilter) -> Integer

Return the number of vanishing moments of the boundary scaling functions defined by `F`.
"""->
function van_moment(F::BoundaryFilter)
	return length(F)
end

const LEFT_FILTERS = Dict{Int, Array{Any,1}}(
2 => Any[
[ 0.6033325119E+00 ; 0.6908955318E+00 ; -0.3983129977E+00 ]
,
[ 0.3751746045E-01 ; 0.4573276599E+00 ; 0.8500881025E+00 ; 0.2238203570E+00 ; -0.1292227434E+00 ]
]
,
3 => Any[
[ 0.3888997639E+00 ; -0.8820782813E-01 ; -0.8478413085E+00 ; 0.3494874367E+00 ]
,
[ -0.6211483178E+00 ; 0.5225273932E+00 ; -0.2000080031E+00 ; 0.3378673486E+00 ; -0.3997707705E+00 ; 0.1648201297E+00 ]
,
[ -0.9587863831E-02 ; 0.3712255319E-03 ; 0.3260097101E+00 ; 0.8016481645E+00 ; 0.4720552620E+00 ; -0.1400420809E+00 ; -0.8542510010E-01 ; 0.3521962365E-01 ]
]
)


