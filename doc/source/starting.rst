Getting Started
===============

Installation
------------

The *GeneralizedSampling* package is available through the Julia package system by running ``Pkg.add("GeneralizedSampling")``.


Background
----------

When sampling a function/signal, one often has to use a basis that is inferior for reconstruction, i.e., a basis requires a large number of elements to give a decent approximation.
Think e.g. of the `Gibbs phenomenon <https://en.wikipedia.org/wiki/Gibbs_phenomenon>`_ that occurs when using the Fourier basis for reconstructing discontinuous functions.

Generalized sampling is a general technique for transforming samples of a function wrt. one basis into samples wrt. another basis, i.e., essentially performing a change of basis.

The *GeneralizedSampling* package currently implements transformations from the Fourier to the Wavelet basis on :math:`L^2([0,1])`.


Basic Usage
-----------

The primary element in *GeneralizedSampling* is a change of basis type (`CoB`) computed from the **sample locations**, the name of wavelet used for reconstruction and the scale `J` of the wavelet space:

.. code-block:: julia

    T = freq2wave(samples, wavename, J)

The :code:`T` object behaves in many ways like the ordinary matrix it resemples.
In particular, multiplication, multiplication with its adjoint and the backslash operator for least squares solutions work as expected.
So if the Fourier transform :code:`Ghat` of a function :code:`G` is sampled in the locations used for :code:`T`, :code:`G`'s representation :code:`w` in the wavelet domain of scale `J` is approximated by

.. code-block:: julia

    f = Ghat(samples)
    w = T \ f

To evaluate ``w`` in the wavelet basis, the `WaveletPlots package <https://github.com/robertdj/WaveletPlots.jl>`_ can be used:

.. code-block:: julia

    using WaveletPlots
    y = weval(w, wavename)


**Note**: 

- The theory of generalized sampling promises that the change of basis matrix :math:`T` is numerically stable, *if* the number of samples are sufficiently high compared to :math:`J`.
- The condition number of :math:`T` is not directly available. To compute the condition number the change of basis matrix has to be computed explicitly with ``collect(T)``.
- The change of basis matrix may very well be too large to compute explicitly; check ``size(T)`` before collecting.

