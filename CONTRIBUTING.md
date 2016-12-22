# Contributing

Pull requests with bug fixes are very welcome.


## Extending functionality

The theory of generalized sampling is very generic, but the package currently focuses on Fourier bases/frames and wavelet bases.

If you would like to contribute with functions for computing new representations the following is needed:

- A `type` for your change of basis matrix (like `Freq2Wave`)
- The functions `Base.A_mul_B!` and `Base.Ac_mul_B!` for in-place multiplication should be overloaded for your new type.
- Miscellaneous functions for your new type would be appreciated. E.g. `size`, `eltype` and `collect`.
- Include the new functions in a subfolder of `src` with an indicative name.

Remember to include tests for the new types/functions and make sure the tests pass.

