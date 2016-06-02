function NFFT.nfft!{T}(p::NFFTPlan{2}, f::AbstractVector{T}, fHat::Vector{T}, d::Integer)
	tmpVec = slice(p.tmpVec, :, 1)
	fill!(tmpVec, zero(T))

	@inbounds apodization!(p, f, tmpVec, d)
	fft!(tmpVec)
	fill!(fHat, zero(T))
	@inbounds convolve!(p, tmpVec, fHat, d)
	return fHat
end

function apodization!{T}(p::NFFTPlan{2}, f::AbstractVector{T}, g::AbstractVector{T}, d::Integer)
  n = p.n[d]
  N = p.N[d]
  const offset = round( Int, n - N / 2 ) - 1
  for l in 1:N
    g[((l+offset)% n) + 1] = f[l] * p.windowHatInvLUT[d][l]
  end
end

function convolve!{T}(p::NFFTPlan{2}, g::AbstractVector{T}, fHat::AbstractVector{T}, d::Integer)
  n = p.n[d]

  for k in 1:p.M # loop over nonequispaced nodes
    c = floor(Int,p.x[d,k]*n)
    for l in (c-p.m):(c+p.m) # loop over nonzero elements

      idx = ((l+n)% n) + 1
      idx2 = abs(((p.x[d,k]*n - l)/p.m )*(p.K-1)) + 1
      idx2L = floor(Int,idx2)

      fHat[k] += g[idx] * (p.windowLUT[d][idx2L] + ( idx2-idx2L ) * (p.windowLUT[d][idx2L+1] - p.windowLUT[d][idx2L] ) )
    end
  end
end

