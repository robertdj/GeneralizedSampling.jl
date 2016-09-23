function myfft!(A::StridedVector{Complex{Float64}}, P::FFTPlan, B::AbstractVector{Complex{Float64}})
	length(A) == P.M || throw(DimensionMismatch())
	length(B) == P.N || throw(DimensionMismatch())

	fill!(P.tmpVec, ComplexZero)
	for l in 1:P.N
		idx = 1 + P.q*(l - 1)
		@inbounds P.tmpVec[idx] = B[l] * P.pre_phaseshift[l]
	end

	P.forwardFFT * P.tmpVec

	copy!(A, 1, P.tmpVec, 1, P.M)
	had!(A, P.post_phaseshift)

	return A
end

@inline myfft!(A::StridedArray, P::NFFT.NFFTPlan, B::AbstractArray) = NFFT.nfft!(P, B, A)

