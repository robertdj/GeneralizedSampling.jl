#=
Compute 1D NFFT on each column/row of a matrix.
# TODO: Not very elegant, try to mimick the mapreduce approach of e.g. sum
=#

function NFFT.nfft!{T}(p::NFFTPlan{1}, f::AbstractMatrix{T}, fHat::StridedMatrix{T}, ::Type{Val{1}})
	( Nrows = size(f,2) ) == size(fHat,2) || throw(DimensionMismatch())

	for rowidx in 1:Nrows
		frow = slice(f, :, rowidx)
		fhatrow = slice(fHat, :, rowidx)

		NFFT.nfft!(p, frow, fhatrow)
	end

	return fHat
end

function NFFT.nfft!{T}(p::NFFTPlan{1}, f::AbstractMatrix{T}, fHat::StridedMatrix{T}, ::Type{Val{2}})
	( Nrows = size(f,1) ) == size(fHat,2) || throw(DimensionMismatch())

	for rowidx in 1:Nrows
		frow = slice(f, rowidx, :)
		fhatrow = slice(fHat, :, rowidx)

		NFFT.nfft!(p, frow, fhatrow)
	end

	return fHat
end


function NFFT.nfft_adjoint!{T}(p::NFFTPlan{1}, fHat::AbstractMatrix{T}, f::StridedMatrix{T}, ::Type{Val{1}})
	( Nrows = size(f,2) ) == size(fHat,2) || throw(DimensionMismatch())

	for rowidx in 1:Nrows
		frow = slice(f, :, rowidx)
		fhatrow = slice(fHat, :, rowidx)

		NFFT.nfft_adjoint!(p, fhatrow, frow)
	end

	return f
end

function NFFT.nfft_adjoint!{T}(p::NFFTPlan{1}, fHat::AbstractMatrix{T}, f::StridedMatrix{T}, ::Type{Val{2}})
	( Nrows = size(f,1) ) == size(fHat,2) || throw(DimensionMismatch())

	for rowidx in 1:Nrows
		frow = slice(f, rowidx, :)
		fhatrow = slice(fHat, :, rowidx)

		NFFT.nfft_adjoint!(p, fhatrow, frow)
	end

	return f
end

