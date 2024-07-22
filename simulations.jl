using Statistics, Printf, Logging
using OrdinaryDiffEq
using LinearAlgebra: norm
using Plots

pyplot(dpi=300)               #, size=0.5.*(1440, 900))#, thickness_scaling=2.2)
closeall()

## Comment out the following line to disable numerous debug messages and plots
# Logging.global_logger(Logging.ConsoleLogger(stderr, Logging.Debug))

## some useful definitions
include("init.jl")

## Physical constants & conversion factors
include("constants.jl")

## Abstractions
include("api/types.jl")
include("api/methods.jl")
include("api/uapprox.jl")

## Tools
include("tools/wp.jl")
include("tools/gaussians.jl")
include("tools/spectra.jl")

## Model potentials
include("models/morse.jl")      # Morse potential
include("models/ammonia.jl")    # Ammonia double-well potential
include("models/quartic.jl")    # Quartic double-well potential

## Dynamical models
include("dynamics/dynamics.jl")

# ## Tests
# include("tests/H2_Fourier_grid_Hamiltonian.jl")
# include("tests/ammonia_maser.jl")
# include("tests/classical_harmonic_oscillator.jl")
# include("tests/classical_morse_oscillator.jl")

## paths
datasave_path = "./output/data"
figures_path  = "./output/figs"

################################# Calculations #################################
figure = mkfigure()             # For online plotting

include("sec_4_1.jl")
include("sec_4_2.jl")
include("sec_4_3.jl")
include("sec_4_4.jl")
include("sec_4_5.jl")
include("sec_A_1.jl")
include("sec_S1.jl")            # Not included in the manuscript

println("Hit Enter to continue...")
readline()

nothing
# EOF
