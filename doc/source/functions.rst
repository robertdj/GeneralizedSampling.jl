Functions
=========

The functions provide by *GeneralizedSampling* are divided into three categories: 
High level functions related to change of basis types, functions for computing Fourier transforms used in construction change of basis objects and miscellaneous.

Furthermore, the type hierarchy of change of basis objects is introduced.


.. _allowedwavelets:

Wavelets
~~~~~~~~

Currently *GeneralizedSampling* supports reconstruction in Daubechies wavelets/scaling functions.
As the reconstruction happens on :math:`[0,1]` the functions near the boundaries needs to be modified -- which can happen in multiple ways.
We have chosen the boundary wavelets from :cite:`Cohen:Daubechies:Vial:1993`, which has the same number of vanishing moments as the internal/non-boundary wavelets.

The allowed wavelets are named "haar", "db1", "db2", ..., "db8".


Change of basis
---------------

.. function:: freq2wave(samples, wavename, J, B, ...)

    Compute a frequency-to-wavelet change of basis object.

    - ``samples`` are the sampling **locations**. For 1D samples it a vector and for 2D samples it is matrix with 2 columns.
    - ``wavename`` is a string as described in in the :ref:`allowedwavelets` section.
    - ``J`` is the scale of the reconstruction wavelets (the :math:`V_J` space in multiresolution terminology). Note that :math:`2^J` has to be larger than the length of the wavelet's support.
    - If ``samples`` are *not* a uniform grid, a bandwidth ``B`` has to be supplied that is larger than ``maxabs(samples)``. Note that if ``B`` is too large the *density*  of the samples may also be too large, which degenerates the condition number of ``T``
    - The optional arguments are passed to the functions computing :ref:`fourier` of the wavelet (if needed).

    As mentioned in :ref:`starting`, the benefit of generalized sampling is that the computations are numerically stable.
    However, some assumptions must be fulfilled:
    
    - There should be more samples than reconstruction coefficients. The ratio between samples and the number of reconstruction coefficients that ensures a numerically stable matrix is called the *staple sampling rate*. For uniform samples the stable sampling rate is well described -- see :ref:`references`. For non-uniform samples the staple sampling rate also depends on the *density* of the samples, which is defined as the minimum radius that gives a covering of the bandwidth area with equal sized circles centered at the sampling points.
    - The samples should be distributed around the origin, i.e., only positive samples does not work.

.. function:: collect(T)

    Compute the explicit change of basis matrix from ``T``.
    Note that this is time consuming and possible impossible to hold in memory for large sampling/reconstrution sets.

.. function:: \

    Overload of the usual backslash operator:
    ``x = T \ b``
    computes the least squares solution of :math:`Tx \approx b`.

.. function:: \*, '\*

    Overload of multiplication with ``T`` and ``T'``, the adjoint of ``T``.

.. function:: size(T)

    Get the size of ``T`` as a tuple :math:`(M,N)`.

.. function:: size(T, d)

    Get the size along dimension ``d`` of ``T``.

.. function:: wsize(T)

    Get the size of the reconstructed wavelets as a tuple.
    In 1D the result is :math:`(N,)` and in 2D the result is :math:`(N,N)`.

.. function:: isuniform(T)

    Returns ``true`` if the samples used for ``T`` are uniform and ``false`` otherwise.

.. function:: hasboundary(T)

    Returns ``true`` if the wavelet used for reconstruction in ``T`` has special functions near boundaries and ``false`` otherwise.

.. function:: van_moment(T)

    Get the number of vanishing moments of the wavelet used for reconstruction in ``T``.


Types
~~~~~

.. code-block:: julia

    CoB
    Freq2Wave <: CoB
    Freq2NoBoundaryWave <: Freq2Wave
    Freq2BoundaryWave <: Freq2Wave


.. _fourier:

Fourier transforms
------------------

Fourier transforms of the scaling functions are available.
The high level interface is

.. function:: FourScalingFunc(xi, wavename, J, k; ...)

    Evaluate the Fourier transform of ``wavename`` at ``xi``.

    - ``xi`` is either a real number of an array of real numbers.
    - ``wavename`` is a string as described in in the :ref:`allowedwavelets` section.
    - Optional ``J`` is the scale of the scaling function, which by default is 0.
    - Optional ``k`` is the translation of the scaling function, which by default is 0.

    The remaining arguments relate to the iterative computations of the Fourier transforms and are usually not needed. 
    Check the inside documentation for more info.

The lower level functions are available for each type of scaling function.
The common interface is that the first arguments are the point(s) at which to evaluate the Fourier transform and specification(s) of scaling function. 
The ``J`` and ``k`` arguments are the scale and translation, respectively, as above.

.. function:: FourHaarScaling(xi, J, k)

    Fourier transform of the Haar scaling function (which is exact).

.. function:: FourDaubScaling(xi, C, J, k; ...)

    Fourier transform of the Daubechies scaling function defined by the filter vector ``C``..
    The filter ``C`` must sum to 1.

.. function:: FourDaubScaling(xi, N, J, k; ...)

    Fourier transform of the Daubechies scaling function with ``N`` vanishing moments.

.. function:: FourDaubScaling(xi, N, side, J; ...)

    Fourier transform of the Daubechies scaling function with ``N`` vanishing moments adapted to the left (``'L'``) or right (``'R'``) boundaries.


Miscellaneous
-------------


