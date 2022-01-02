# RegularizationTools

*A Julia package to perform Tikhonov regularization for small to moderate size problems.*

| **Documentation**                                           |  **Citations** |
|:-------------------------------------------------------------------------------:|:-------------------------------------------------------------------------------------------------------------------------------------------:|
| [![][docs-stable-img]][docs-stable-url] | [![DOI](https://img.shields.io/badge/DOI-10.5194%2Famt--14--7909--2021-blue)](https://doi.org/10.5194/amt-14-7909-2021)  |

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

## Documentation

- [**STABLE**][docs-stable-url] &mdash; **documentation of the most recently tagged version.**

## Project Status
The package is being developed for, Julia `1.5` and above.

## Citations

*If you use this package and publish your work, please cite:*

Petters, M. D.: Revisiting matrix-based inversion of scanning mobility particle sizer (SMPS) and humidified tandem differential mobility analyzer (HTDMA) data, Atmos. Meas. Tech., 14, 7909–7928, https://doi.org/10.5194/amt-14-7909-2021, 2021.

This work was supported by the United States Department of Energy, Office of Science, Biological and Environment Research, Grant number DE-SC0021074.

[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-dev-url]: https://mdpetters.github.io/RegularizationTools.jl/latest/

[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[docs-stable-url]: https://mdpetters.github.io/RegularizationTools.jl/stable/
