type Freq2wave1D
	# Sampling
	samples::Vector{Float64}
	weights::Union{Bool, Vector{Float64}}

	# Reconstruction
	wave::String
	column::Vector{Complex{Float64}}
	J::Int

	# Multiplication with FFT: T*x = diag*NFFT(x)
	diag::Vector{Complex{Float64}}
	FFT::NFFT.NFFTPlan

	#Freq2wave1D() = new()
end


@doc """
	freq2wave(samples::Vector, wave::String, J::Int)

Make change of basis for switching between frequency responses and wavelets.

- `samples` are the sampling locations.
- `wave` is the name of the wavelet; see documentation for possibilities.
- `J` is the scale of the wavelet transform to reconstruct.
"""->
function freq2wave(samples::Vector, wave::String, J::Int)
	M = length(samples)
	# TODO: Warning if J is too big

	# NFFTPlans: Frequencies must be in the torus [-1/2, 1/2)
	N = 2^J
	xi = scale(samples, 1/N)
	xi_frac = frac(xi)
	p = NFFTPlan(xi_frac, N)

	# Evaluate the first column in change of basis matrix
	func = string("Four", wave, "Scaling")
	basis = eval(parse(func))( samples, J, 0 )
	diag = basis .* cis(-pi*N*xi_frac)

	if isuniform(samples)
		weights = false
	else
		weights = weights(samples, B)
	end

	T = Freq2wave1D(samples, weights, wave, basis, J, diag, p)

	return T
end


function Base.show(io::IO, C::Freq2wave1D)
	println(io, "1D change of basis matrix")

	typeof(C.weights) == Bool ?  U = " " : U = " non-"

	println(io, "From: ", length(C.samples), U, "uniform frequency samples")
	print(io, "To: ", 2^C.J, " ", C.wave, " wavelets")
end

