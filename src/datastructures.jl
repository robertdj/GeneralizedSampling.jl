type freq2wave
	# Sampling
	samples::Vector{Float64}
	weights::Union{Bool, Vector{Float64}}
	M::Int

	# Reconstruction
	#= wave::Function =#
	wave::String
	basis::Vector{Complex{Float64}}
	J::Int
end

@doc """
	freq2wave(samples::Vector, wave::String, J::Int)

Make change of basis type for switching between frequency responses and wavelets.

- `samples` are the sampling locations.
- `wave` is the name of the wavelet; see documentation for possibilities.
- `J` is the scale of the wavelet transform to reconstruct.
"""->
function freq2wave(samples::Vector, wave::String, J::Int)
	if isuniform(samples)
		weights = false
	else
		weights = weights(samples, B)
	end

	M = length(samples)

	#= basis = wave(samples, J) =#
	basis = FourHaarWavelet(samples, J, 0)

	C = freq2wave(samples, weights, M, wave, basis, J)
end


function Base.show(io::IO, C::freq2wave)
	if C.weights == false
		U = " uniform frequency samples"
	else
		U = " non-uniform frequency samples"
	end

	println(io, "Change of basis matrix:")
	println(io, "From ", C.M, U, " to ", 2^C.J, " ", C.wave, " wavelets")
end

