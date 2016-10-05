@inline myfft!(A::StridedArray, P::NFFT.NFFTPlan, B::AbstractArray) = NFFT.nfft!(P, B, A)
@inline myfft_adjoint!(A::StridedArray, P::NFFT.NFFTPlan, B::AbstractArray) = NFFT.nfft_adjoint!(P, B, A)


@generated function myfft!{D}(A::StridedVector{Complex{Float64}}, 
                              P::FFTPlan{D}, B::AbstractArray{Complex{Float64}, D})
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


@generated function myfft_adjoint!{D}(A::StridedArray{Complex{Float64}, D}, 
                                      P::FFTPlan{D}, B::AbstractVector{Complex{Float64}})
    quote
        size(A) == P.N || throw(DimensionMismatch())
        length(B) == prod(P.M) || throw(DimensionMismatch())

        fill!(P.FFTvec, ComplexZero)
        FFTvec = P.FFTvec
        idx_0 = 0
        @nloops $D l d->1:P.M[d] begin
            idx_0 += 1
            (@nref $D FFTvec l) = B[idx_0] * conj(P.post_phaseshift[idx_0])
        end

        P.backwardFFT * FFTvec

        pre_phaseshift = P.pre_phaseshift
        @nexprs $D d -> idx_d = 1
        @nloops $D l d->1:P.N[d] d->idx_{d-1}=1 d->idx_d+=P.q[d] begin
            (@nref $D A l) = (@nref $D FFTvec idx) * conj( (@nref $D pre_phaseshift l) )
        end

        return A
    end
end

