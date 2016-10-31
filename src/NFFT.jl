myfft!(fHat::StridedArray, P::NFFT.NFFTPlan, f::AbstractArray) = NFFT.nfft!(P, f, fHat)
myfft_adjoint!(f::StridedArray, P::NFFT.NFFTPlan, fHat::AbstractArray) = NFFT.nfft_adjoint!(P, fHat, f)


@generated function myfft!{D}(fHat::StridedVector{Complex{Float64}}, 
                              P::FFTPlan{D}, f::AbstractArray{Complex{Float64}, D})
    quote
        length(fHat) == prod(P.M) || throw(DimensionMismatch())
        size(f) == P.N || throw(DimensionMismatch())

        fill!(P.FFTvec, ComplexZero)
        FFTvec = P.FFTvec
        pre_phaseshift = P.pre_phaseshift
        @nloops $D l d->1:P.N[d] d->idx_d = 1 + P.q[d]*(l_d - 1) begin
            (@nref $D FFTvec idx) = (@nref $D f l) * (@nref $D pre_phaseshift l)
        end

        P.forwardFFT * FFTvec

        idx = 0
        @nloops $D l d->1:P.M[d] begin
            idx += 1
            fHat[idx] = (@nref $D FFTvec l) * P.post_phaseshift[idx]
        end

        return fHat
    end
end


@generated function myfft_adjoint!{D}(f::StridedArray{Complex{Float64}, D}, 
                                      P::FFTPlan{D}, fHat::AbstractVector{Complex{Float64}})
    quote
        size(f) == P.N || throw(DimensionMismatch())
        length(fHat) == prod(P.M) || throw(DimensionMismatch())

        fill!(P.FFTvec, ComplexZero)
        FFTvec = P.FFTvec
        idx_0 = 0
        @nloops $D l d->1:P.M[d] begin
            idx_0 += 1
            (@nref $D FFTvec l) = fHat[idx_0] * conj(P.post_phaseshift[idx_0])
        end

        P.backwardFFT * FFTvec

        pre_phaseshift = P.pre_phaseshift
        @nexprs $D d -> idx_d = 1
        @nloops $D l d->1:P.N[d] d->idx_{d-1}=1 d->idx_d+=P.q[d] begin
            (@nref $D f l) = (@nref $D FFTvec idx) * conj( (@nref $D pre_phaseshift l) )
        end

        return f
    end
end

