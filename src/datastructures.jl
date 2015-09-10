type Freq2wave1D
	# Sampling
	samples::Vector{Float64}
	weights::Union{Bool, Vector{Float64}}

	# Reconstruction
	wave::String
	column::Vector{Complex{Float64}}
	# TODO: J is redundant; remove?
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
	# TODO: Parse strings such as "db4"
	func = string("Four", wave, "Scaling")
	basis = eval(parse(func))( samples, J, 0 )
	diag = basis .* cis(-pi*N*xi_frac)

	if isuniform(samples)
		weights = false
	else
		# TODO: Bandwidth?
		weights = weights(samples, B)
	end

	return Freq2wave1D(samples, weights, wave, basis, J, diag, p)
end


function Base.show(io::IO, T::Freq2wave1D)
	println(io, "1D change of basis matrix")

	typeof(T.weights) == Bool ?  U = " " : U = " non-"

	println(io, "From: ", length(T.samples), U, "uniform frequency samples")
	print(io, "To: ", 2^T.J, " ", T.wave, " wavelets")
end


# ------------------------------------------------------------
# Basic operations for Freq2wave1D

@doc """
Compute 

	Freq2wave1D * vector
"""->
function *(T::Freq2wave1D, x::Vector{Complex{Float64}})
	@assert 2^T.J == length(x)

	Tx = nfft(T.FFT, x)
	had!(Tx, T.diag)

	return Tx
end

#= function *{T<:Number}(T::Freq2wave1D, x::Vector{T}) =#
#= 	z = complex(float(x)) =#
function *(T::Freq2wave1D, x::Vector{Float64})
	z = complex(x)
	return T*z
end

