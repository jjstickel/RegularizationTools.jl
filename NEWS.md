# RegularizationTools.jl NEWS

#### Notes on release changes and ongoing development
---

## v0.5.0
- Significant performance gain in "setupRegularizationProblem" by switching from GSVD to QR factorization approach
- Update documentation
- Update tests to be independent of random number generation
- Update Project and Manifest to be current. Added CSV and DataFrames as deps for unit tests
- Explicitly account for the unregularized component of x which is not affected by the regularization scheme (fix for issue #7 raised by Jonathan Stickel). Note that the bug/issue affected solutions for very larger regularization parameters, but had no/neglible effect on calculations near the optimal solution for regularization. 
- Update code for smoothing matrix generation (base on example in issue #7 contributed by Jonathan Stickel). The updated code also directly allows computation of arbitrary order of the smoothing matrix.
- Added smoothing example from issue #7 to examples folder

## v0.4.1
- Minor release associated with julia 1.6.0 
- Update Project and Manifest
- Make tests of invert function resilient against changes in random number generator
- Change order in scattering example to 2

## v0.4.0
- Add high-level API invert function to simplify notation

## v0.3.1
- Add memoization
- Merge CompatHelper: add new compat entry for "Memoize" at version "0.4"
- Fix various bugs related to non-square matrices

## v0.3.0
- Add abstract interface to generate design matrix from a forward model.
- Add documentation and examples how to use the interface
- Add fallback in solve function in case cholesky factorization fails.
- Fix bug when creating L matrix for non-square design matrices.

## v0.2.0
- Add constraint minimization solver to enfore upper and lower bound for some problems.
- Add documentation.

## v0.1.1
- Fix standard form conversion error. Solution is now based on generalized SVD. 
- Touch up documentation.
- RegularizationTools is now part of the general registry

## v0.1.0
- Initial Release
