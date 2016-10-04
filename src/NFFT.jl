@inline myfft!(A::StridedArray, P::NFFT.NFFTPlan, B::AbstractArray) = NFFT.nfft!(P, B, A)
@inline myfft_adjoint!(A::StridedArray, P::NFFT.NFFTPlan, B::AbstractArray) = NFFT.nfft_adjoint!(P, B, A)


function myfft!(A::StridedVector{Complex{Float64}}, P::FFTPlan, B::AbstractVector{Complex{Float64}})
    length(A) == prod(P.M) || throw(DimensionMismatch())
	size(B) == P.N || throw(DimensionMismatch())

    fill!(P.FFTvec, ComplexZero)
    for l in 1:P.N[1]
        idx = 1 + P.q[1]*(l - 1)
        P.FFTvec[idx] = B[l] * P.pre_phaseshift[l]
    end

    P.forwardFFT * P.FFTvec

    copy!(A, 1, P.FFTvec, 1, P.M[1])
    had!(A, P.post_phaseshift)

    return A
end

function myfft!(A::StridedVector{Complex{Float64}}, P::FFTPlan, B::AbstractMatrix{Complex{Float64}})
    length(A) == prod(P.M) || throw(DimensionMismatch())
	size(B) == P.N || throw(DimensionMismatch())

    fill!(P.FFTvec, ComplexZero)
    for l_2 in 1:P.N[2]
        idx_2 = 1 + P.q[2]*(l_2 - 1)
        for l_1 in 1:P.N[1]
            idx_1 = 1 + P.q[1]*(l_1 - 1)
            P.FFTvec[idx_1, idx_2] = B[l_1, l_2] * P.pre_phaseshift[l_1, l_2]
        end
    end

    P.forwardFFT * P.FFTvec

    idx = 0
    for l_2 in 1:P.M[2]
        for l_1 in 1:P.M[1]
            idx += 1
            A[idx] = P.FFTvec[l_1, l_2] * P.post_phaseshift[idx]
        end
    end

    return A
end


function myfft_adjoint!(A::StridedVector{Complex{Float64}}, P::FFTPlan, B::AbstractVector{Complex{Float64}})
	size(A) == P.N || throw(DimensionMismatch())
    length(B) == prod(P.M) || throw(DimensionMismatch())

    fill!(P.FFTvec, ComplexZero)
    for l in 1:P.M[1]
        P.FFTvec[l] = B[l] * conj(P.post_phaseshift[l])
    end

    P.backwardFFT * P.FFTvec

    idx = 1
    for l in 1:P.N[1]
        A[l] = P.FFTvec[idx] * conj(P.pre_phaseshift[l])
        idx += P.q[1]
    end

    return A
end

