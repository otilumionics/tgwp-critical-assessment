# Overview

A collection of [Julia](https://julialang.org/ "Julia website")
scripts to generate numerical data and plots for the paper "Thawed
Gaussian wave packet dynamics: a critical assessment of three
propagation schemes"[^1].

[^1]: Ilya G. Ryabinkin, Rami Gherib, and Scott N. Genin,
<https://arxiv.org/abs/2405.01729>

# Prerequisites

To be able to run the code you need to install
[Julia](https://julialang.org/) and
[Gnuplot](http://www.gnuplot.info/) binaries to your system. In any
modern Linux distribution, such as Ubuntu and derivatives, Fedora,
Debian, Red Hat and derivatives or others, these two programs are
available in repositories. Specifically, the code was tested on
`Julia` versions 1.10.* and 1.5.3 and `Gnuplot` version 5.2.8.

In `Julia` the following packages must be installed:

```text
FFTW
FastGaussQuadrature
NLopt
OrdinaryDiffEq
Plots
Polynomials
PyPlot
SpecialPolynomials

```

This could be done interactively by hitting `]` in REPL and using `add
<PackageName>` command, or programmatically as

```julia
using Pkg

Pkg.add("<PackageName>")
```

The following package is optional, needed only if tests are executed:

```text
PrettyTables
```

# Source structure

```code
.
├── api
│   ├── methods.jl
│   ├── types.jl
│   └── uapprox.jl
├── dynamics
│   └── dynamics.jl
├── models
│   ├── ammonia.jl
│   ├── morse.jl
│   └── quartic.jl
├── output
│   ├── data
│   └── figs
├── tests
│   ├── ammonia_maser.jl
│   ├── classical_harmonic_oscillator.jl
│   ├── classical_morse_oscillator.jl
│   └── H2_Fourier_grid_Hamiltonian.jl
├── tools
│   ├── gaussians.jl
│   ├── spectra.jl
│   └── wp.jl
├── constants.jl
├── init.jl
├── LICENSE
├── README.md
├── sample_output.debug.log
├── sample_output.log
├── sec_4_1.jl
├── sec_4_2.jl
├── sec_4_3.jl
├── sec_4_4.jl
├── sec_4_5.jl
├── sec_A_1.jl
├── sec_S1.jl
└── simulations.jl
```

# Installation & running

There is no specific installation routine. Just download files in a
separate directory. The main script is `sumulations.jl`. It can be
interactively loaded into `Julia` REPL as

```julia
julia> include("simulations.jl")
```

or executed non-interactively as

```bash

julia simulations.jl 2>&1
```

After execution the sub-directories `./output/data` and
`./output/figs` will be populated with numerical outputs and plots,
respectively. There are 7 figures in the paper; all except the 5-th
are produced in `EPS` format, the figure 5 is produced as `PDF` to
support transparency. `EPS` figures can be converged to `PDF` by
running

```bash
for file in *.eps; do epstopdf ${file}; done
```

`epstopdf` utility is distributed with `TeXLive`; it is available for
Windows, Linux, and Mac computers.

# Customization

`sumulations.jl` can be customized to provide more or less output
including plots. All interactive plots are enabled only with debug
printout. To switch between debug mode, user can (un)comment out line
10 in `simulations.jl`:

```julia
Logging.global_logger(Logging.ConsoleLogger(stderr, Logging.Debug))
```

# Sample output

Sample outputs (excluding plots) can be found in `sample_output.log`
and `sample_output.debug.log`.
