Functions
=========

The functions provide by *GeneralizedSampling* are divided into three categories: 
High level functions related to change of basis types, functions for computing Fourier transforms used in construction change of basis objects and miscellaneous.

Furthermore, the type hierarchy of change of basis objects is introduced.


.. _allowedwavelets:

Wavelets
--------

Currently *GeneralizedSampling* supports reconstruction in Daubechies wavelets/scaling functions.
As the reconstruction happens on :math:`[-1/2,1/2]` the functions near the boundaries needs to be modified -- which can happen in multiple ways.
We have chosen the boundary wavelets from :cite:`Cohen:Daubechies:Vial:1993`, which has the same number of vanishing moments as the internal/non-boundary wavelets.

The allowed wavelets are named "haar", "db1", "db2", ..., "db8".


.. _CoB:

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
    
    - There should be more samples than reconstruction coefficients. The ratio between samples and the number of reconstruction coefficients that ensures a numerically stable matrix is called the *stable sampling rate*. For uniform samples the stable sampling rate is well described -- see :ref:`references`. For non-uniform samples the stable sampling rate also depends on the *density* of the samples, which is defined as the minimum radius that gives a covering of the bandwidth area with equal sized circles centered at the sampling points.
    - The samples should be distributed around the origin, i.e., only positive samples does not work.

    Example:

    .. code-block:: julia
    
        samples = grid(2^7, 0.5)
        T = freq2wave(samples, "db2", 6)

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
-----

The abstract change of basis supertype is denoted ``CoB``.

The specific change of basis types implemented are from Fourier to wavelet bases. 
They are collectively denoted ``Freq2Wave`` and are a subtype of ``CoB``:

.. code-block:: julia

    Freq2Wave <: CoB

The computations for wavelets with boundary correction are more involved than for those without and therefore two subtypes of ``Freq2Wave`` are introduced for both 1D and 2D:

.. code-block:: julia

    Freq2NoBoundaryWave1D <: Freq2Wave
    Freq2BoundaryWave2D <: Freq2Wave
    Freq2NoBoundaryWave1D <: Freq2Wave
    Freq2BoundaryWave2D <: Freq2Wave


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

    As an example, the following command computes the Fourier transform of the Daubechies 2 scaling function and plots the real and imaginary part using `Winston <https://github.com/nolta/Winston.jl>`_:

    .. code-block:: julia
    
        x = linspace(-5, 5, 1000)
        y = FourDaubScaling(x, "db2")
        using Winston
        plot(x, real(y), x, imag(y))

.. function:: FourScalingFunc(xi, wavename, side, J, k; ...)

    As above, but for the boundary scaling functions.
    ``side`` is either ``'L'`` or ``'R'``.

    Note that these Fourier transforms are for the scaling functions that in the time domain are translated to fir the reconstruction interval :math:`[-1/2,1/2]`, i.e., their Fourier transforms are phase shifted.

The lower level functions are available for each type of scaling function, but not documented here. 
Check the documentation in Julia with the usual ``?function`` where ``function`` is ``FourHaarScaling`` or ``FourDaubScaling``.


Miscellaneous
-------------

Functions that are used for internal documentation are not documented here; they all have documentation available from within Julia.

To generate sampling locations from a uniformly spaced grid there are functions in 1D and 2D.

.. function:: grid(M, D)

    Return a vector of ``M`` locations evenly distributed around the origin with distance `D`.
    By default, ``D = 1``.

.. function:: grid( (M,N), D )

    Return a matrix with 2 columns containing the ``x``- and ``y``-values of a uniformly distributed grid of locations around the origin with distance ``D``. 
    There are ``M`` different locations in the 1st dimension and ``N`` different locations in the 2nd dimension.

.. function:: isuniform(points)

    Returns ``true`` if ``points`` are located on a uniform grid such as the output from ``grid`` and ``false`` otherwise.


For a configuration of sampling locations ``xi`` the density correcting weights and its density are available as

.. function:: weights(xi, K)

.. function:: density(xi, K)

The bandwidth ``K`` is explained in :ref:`CoB` and must be at least ``maxabs(xi)``.

When dealing with wavelets with boundary corrections, computations differs for the internal and boundary parts.
To this end, the ``split`` function is available to help divide a vector or matrix of coefficients into the parts related to internal/boundary functions.

.. function:: split(x, B)

    Returns three vector *slices* of the ``B`` leftmost, the internal and the ``B`` rightmost entries of ``x``, respectively.

.. function:: split(A, B)

    Returns slices of the outer parts of ``A`` and its internal parts. 
    The outer parts are each of the four :math:`B\times B` corners and each of the four non-corner sides (with one dimension equal to ``B``).

