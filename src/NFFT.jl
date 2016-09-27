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

function myfft!(A::StridedVector{Complex{Float64}}, P::FFTPlan, B::AbstractMatrix{Complex{Float64}})
	fill!(P.tmpVec, ComplexZero)
    for l_2 in 1:P.N[2]
        idx_2 = 1 + P.q[2]*(l_2 - 1)
        for l_1 in 1:P.N[1]
            idx_1 = 1 + P.q[1]*(l_1 - 1)
            @inbounds P.tmpVec[idx_1, idx_2] = B[l_1, l_2] * P.pre_phaseshift[l_1, l_2]
        end
	end

	P.forwardFFT * P.tmpVec

    #= Aview = P.tmpVec[1:P.M[1], 1:P.M[2]] =#
	#= had!(A, Aview, P.post_phaseshift) =#

    idx = 0
    for l_2 in 1:P.M[2]
        for l_1 in 1:P.M[1]
            idx += 1
            A[idx] = P.tmpVec[l_1, l_2] * P.post_phaseshift[idx]
        end
    end

	return A
end

@inline myfft!(A::StridedArray, P::NFFT.NFFTPlan, B::AbstractArray) = NFFT.nfft!(P, B, A)

