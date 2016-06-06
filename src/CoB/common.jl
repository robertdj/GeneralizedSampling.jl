@doc """
	size(CoB) -> (M,N)
	size(CoB, d) -> M or N

`M` is the number of samples and `N` is the number of reconstruction functions.
"""->
function Base.size(T::CoB)
	( size(T,Val{1}), size(T,Val{2}) )
end

function Base.size(T::CoB, d::Int)
	#d == 1 || d == 2 || throw(DomainError())
	size(T, Val{d})
end

