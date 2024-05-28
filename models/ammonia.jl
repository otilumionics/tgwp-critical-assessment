module AmmoniaModel
using SpecialPolynomials
using ..Constants

using Main: Potential, Exact
import Main: U

# Module does not export, only adds methods

# Atomic units below
const ħ = 1
const m = 2.561*Constants.Da2me # Reduced mass, ≈ 4668 mₑ

struct AmmoniaPotential <:Potential
    k::Float64
    b::Float64
    c::Float64
#=
           Ammonia double-well potential: parabola + a Gaussian barrier

V(x) = 1/2k x^2 + b exp(-cx^2)

References:

 [1] J. Chem. Phys. 36, 1914–1918 (1962), https://doi.org/10.1063/1.1701290 (orig)
 [2] https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Quantum_Tutorials_(Rioux)/04%3A_Spectroscopy/4.04%3A_The_Ammonia_Inversion_and_the_Maser
=#
    AmmoniaPotential(;
                     k = 0.07598,
                     b = 0.05684,
                     c = 1.3696) = new(k, b, c)
end

# Derivatives of b*exp(-c*x^2)
y(x, b, c; der=0) = b*(-sqrt(c))^der*basis(Hermite, der)(sqrt(c)*x)*exp(-c*x^2)

function (p::AmmoniaPotential)(x::Number; der::Int=0)
    # For convenience...
    k = p.k
    b = p.b
    c = p.c

    if der == 0
        return k*x^2/2 + b*exp(-c*x^2)
    elseif der == 1
        return k*x + b*(-2c*x)*exp(-c*x^2)
    elseif der == 2
        return k + b*(-2c + (-2c*x)^2)*exp(-c*x^2)
    elseif der >= 3
        return y(x,b,c; der)
    else
        error("wrong derivative order, der=$(der)")
    end

    return nothing              # Should not be reached

end

"""
    U(approx::Exact{AmmoniaPotential}, q::Number, w::Number; der=(q=0,w=0))

Analytically known extended potential `U(q, w)` and its first-order partial
derivatives -- original derivation/implementation.
"""
function U(approx::Exact{AmmoniaPotential}, q::Number, w::Number; der=(q=0,w=0))
    # @debug approx

    # For convenience...
    k = (approx.p).k
    b = (approx.p).b
    c = (approx.p).c

    t = 1 + 2c*w^2

    if der == (q=0, w=0)
        return ħ^2/m/w^2/8 + k*q^2/2 + k*w^2/2 + b/√t*exp(-c*q^2/t)
    elseif der == (q=1, w=0)
        return k*q + b/√t*(-2c*q/t)*exp(-c*q^2/t)
    elseif der == (q=0, w=1)
        return -ħ^2/m/w^3/4 + k*w +
            b/(t*√t)*exp(-c*q^2/t)*(-2 + 4c*q^2/t)*c*w
    else
        error("not implemented")
    end

    return nothing              # Should not be called
end

end # module AmmoniaModel
