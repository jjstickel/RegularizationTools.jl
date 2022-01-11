# RegularizationTools.jl

*A Julia package to perform Tikhonov regularization for small to moderate size problems.*

RegularizationTools.jl bundles a set routines to compute the [regularized  Tikhonov inverse](https://en.wikipedia.org/wiki/Tikhonov_regularization) using standard linear algebra techniques.  

## Package Features
- Computes the Tikhonov inverse solution with optional boundary constraints.
- Computes optimal regularization parameter using generalized cross validation or the L-curve.
- Solves problems up to ~1000 equations.
- Supports higher order regularization schemes out of the box.
- Supports specifying an *a-priori* estimate of the solution.
- Supports user specified smoothing matrices.
- User friendly interface.
- Extensive documentation.

## About
Tikhonv regularization is also known as Phillips-Twomey-Tikhonov regularization or ridge regression (see Hansen, 2000 for a review). The Web-of-Sciences database lists more than 4500 peer-reviewed publications mentioning "Tikhonov regularization" in the title or abstract, with a current publication rate of ≈350 new papers/year. 

# Quick Start

The package computes the regularized Tikhonov inverse ``{\rm x_{\lambda}}`` by solving the minimization problem 

```math
{\rm {\rm x_{\lambda}}}=\arg\min\left\{ \left\lVert {\bf {\rm {\bf A}{\rm x}-{\rm b}}}\right\rVert _{2}^{2}+\lambda^{2}\left\lVert {\rm {\bf L}({\rm x}-{\rm x_{0}})}\right\rVert _{2}^{2}\right\} 
```

subject to the optional constraint ``{\rm x}_{l}<{\rm x}<{\rm x}_{u}``. Here ``{\rm x_{\lambda}}`` is the regularized estimate of ``{\rm x}``,
``\left\lVert \cdot\right\rVert _{2}`` is the Euclidean norm, ``{\rm {\bf L}}`` is the Tikhonov filter matrix, ``\lambda`` is the regularization parameter, and ``{\rm x_{0}}`` is a vector of an *a-priori* guess of the solution. The initial guess can be taken to be ``{\rm x_{0}}=0`` if no *a-priori* information is known. The solve function searches for the optimal ``\lambda`` and returns the inverse. Optionally, ``{\rm x}_{l}`` and ``{\rm x}_{u}`` can be used to impose boundary constraints on the solution ``{\rm x_{\lambda}}``.

The following script is a minimalist example how to use this package.

```@example
using RegularizationTools, MatrixDepot, Lazy, DelimitedFiles
random(n) = @> readdlm("random.txt") vec x -> x[1:n] 

# This is a test problem for regularization methods
r = mdopen("shaw", 100, false)       # Load the "shaw" problem from MatrixDepot
A, x  = r.A, r.x                     # A is size(100,100), x is length(100)

y = A * x                            # y is the true response 
b = y + 0.2y .* random(100)          # response with superimposed noise
x₀ = 0.4x                            # some a-priori estimate x₀

# Solve 2nd order Tikhonov inversion (L = uppertridiag(−1, 2, −1)) with intial guess x₀
xλ = invert(A, b, Lₖx₀(0, x₀))
include("theory/helpers.jl")   # hide
standard_plot(y, b, x, xλ, x₀) # hide
```

!!! info
	Per [documentation](https://docs.julialang.org/en/v1/stdlib/Random/), the way random numbers are generated in Julia is considered an implementation detail. Bug fixes and speed improvements may change the stream of numbers that are generated after a version change. To keep the examples in the documentation of this package consistent between Julia updates, the function ```random(n)``` is used to load n pre-generated normally distributed random numbers from a text file.

# Installation

The package can be installed from the Julia package prompt with

```julia
julia> ]add RegularizationTools
```

The closing square bracket switches to the package manager interface and the ```add``` command installs the package and any missing dependencies. To return to the Julia REPL hit the ```delete``` key.

To load the package run

```julia
julia> using RegularizationTools
```

For optimal performance, also install the [Intel MKL linear algebra](https://github.com/JuliaComputing/MKL.jl) library.

## Related work
* [MultivariateStats](https://multivariatestatsjl.readthedocs.io/en/stable/index.html): Implements ridge regression without *a priori* estimate and does not include tools to find the optimal regularization parameter.
* [RegularizedLeastSquares](https://tknopp.github.io/RegularizedLeastSquares.jl/latest/): Implements optimization techniques for large-scale scale linear systems.

## Author and Copyright
Markus Petters, Department of Marine, Earth, and Atmospheric Sciences, NC State University
.

# Contribtions
[Jonathan Stickel](https://github.com/jjstickel)
- Several bug fixes
- Improvement of derivative matrix support

# Citation and Acknowledgements
*If you use this package and publish your work, please cite:*

Petters, M. D.: Revisiting matrix-based inversion of scanning mobility particle sizer (SMPS) and humidified tandem differential mobility analyzer (HTDMA) data, Atmos. Meas. Tech., 14, 7909–7928, https://doi.org/10.5194/amt-14-7909-2021, 2021.

This work was supported by the United States Department of Energy, Office of Science, Biological and Environment Research, Grant number DE-SC0021074.
