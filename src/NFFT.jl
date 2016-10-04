@inline myfft!(A::StridedArray, P::NFFT.NFFTPlan, B::AbstractArray) = NFFT.nfft!(P, B, A)
@inline myfft_adjoint!(A::StridedArray, P::NFFT.NFFTPlan, B::AbstractArray) = NFFT.nfft_adjoint!(P, B, A)


@generated function myfft!{D}(A::StridedVector{Complex{Float64}}, P::FFTPlan, B::AbstractArray{Complex{Float64}, D})
    quote
        length(A) == prod(P.M) || throw(DimensionMismatch())
        size(B) == P.N || throw(DimensionMismatch())

        fill!(P.FFTvec, ComplexZero)
        FFTvec = P.FFTvec
        pre_phaseshift = P.pre_phaseshift
        @nloops $D l d->1:P.N[d] d->idx_d = 1 + P.q[d]*(l_d - 1) begin
            (@nref $D FFTvec idx) = (@nref $D B l) * (@nref $D pre_phaseshift l)
        end

        P.forwardFFT * FFTvec

        idx = 0
        @nloops $D l d->1:P.M[d] begin
            idx += 1
            A[idx] = (@nref $D FFTvec l) * P.post_phaseshift[idx]
        end

        return A
    end
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

