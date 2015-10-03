immutable Freq2wave1D
	# Sampling
	samples::Vector{Float64}
	weights::Union{Bool, Vector{Float64}}

	# Reconstruction
	wave::String
	column1::Vector{Complex{Float64}}
	# TODO: J is redundant; remove?
	J::Int

	# Multiplication with FFT: T*x = diag*NFFT(x)
	diag::Vector{Complex{Float64}}
	FFT::NFFT.NFFTPlan{1,Float64}
end


immutable Freq2wave2D
	# Sampling
	samples::Matrix{Float64}
	weights::Union{Bool, Vector{Float64}}

	# Reconstruction
	wave::String
	column1::Vector{Complex{Float64}}
	# TODO: J is redundant; remove?
	J::Int

	# Multiplication with FFT: T*x = diag*NFFT(x)
	diag::Vector{Complex{Float64}}
	FFT::NFFT.NFFTPlan{2,Float64}
end

Freq2wave = Union(Freq2wave1D, Freq2wave2D)


#include("CoB/common.jl")
include("CoB/freq2wave.jl")

