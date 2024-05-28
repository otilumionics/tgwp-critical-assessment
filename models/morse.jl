module MorseModel
using  Main: Potential, Exact
using ..Constants
import Main: U, findmin, energy_levels

# Module does not export, only adds methods

# Atomic units below
const ħ = 1
const m = 1836.15267344         # Reduced mass: 1 proton mass

struct MorsePotential <: Potential
    D::Float64                  # Dissociation limit (if V=0 at min), aka Dₑ
    a::Float64                  # Location of a minimum, aka xₑ
    b::Float64                  # aka β
#=
                  Table II of J. Chem. Phys 130, 134713(2009)
  Montmorillonite, dioctahedral phyllosilicate

  V(x) = D(1 - exp(-b(x-a)))² - D // Shifted by (-D) compared to the original
=#
    MorsePotential(;
                   D = 132.2491/Eh2kcal_per_mol,
                   a = 0.9450/bohr2angs,
                   b = 2.1815*bohr2angs) = new(D, a, b)
end

# n-th derivative of exp(-b(x-a))
y(x, b, a, n::Int) = (-b)^n*exp(-b*(x-a))

function (p::MorsePotential)(x::Number; der::Int=0)
    # @debug p

    if der == 0
        return p.D*((y(x, p.b, p.a, 0) - 1)^2 - 1)
    elseif der >= 1
        return p.D*(y(x, 2p.b, p.a, der) - 2y(x, p.b, p.a, der))
    else
        error("wrong derivative order, der=$(der)")
    end

    return nothing              # Should not be called
end

"""
    U(approx::Exact{MorsePotential}, q::Number, w::Number; der=(q=0,w=0))

Analytically known extended potential and its first-order partial derivatives
for the Morse model. The potential itself was originally derived in Ref.[^1],
derivatives are done here.

[^1]: F. Arickx, J. Broeckhove, E. Kesteloot, L. Lathouwers, and P. V. Leuven,
      “Dynamics of wave packets and the time-dependent variational principle,”
      Chem. Phys. Lett. 128, 310–314 (1986).
"""
function U(approx::Exact{MorsePotential}, q::Number, w::Number; der=(q=0,w=0))
    # @debug approx

    # For convenience...
    D = (approx.p).D
    a = (approx.p).a
    b = (approx.p).b

    if der == (q=0, w=0)
        return ħ^2/m/w^2/8 +
            D*(y(q, 2b, a + w^2*b, 0) - 2y(q, b, (a + w^2/2*b), 0))
    elseif der == (q=1, w=0)
        return D*(y(q, 2b, a + w^2*b, 1) - 2y(q, b, (a + w^2/2*b), 1))
    elseif der == (q=0, w=1)
        return -ħ^2/m/w^3/4 +
            D*(-2w*b*y(q, 2b, a + w^2*b, 1) + 2w*b*y(q, b, (a + w^2/2*b), 1))
    else
        error("not implemented")
    end

    return nothing              # Should not be called
end

"""
    function minimum(V::MorsePotential; guess=nothing, optimizer=nothing)

Morse potential has only one minimum at x = a. Named parameters `guess` and
`optimizer` are dummy.
"""
@inline findmin(V::MorsePotential; guess=nothing, optimizer=nothing) = (-V.D, V.a)

@inline harmonic_frequency(V::MorsePotential, m::Number) = V.b*√(2V.D/m)

"""
    energy_levels(V::MorsePotential, n::Integer; start=0, stop=1, length=1,
    ħ=MorseModel.ħ, m=MorseModel.m)

Energy levels of the Morse potential. `n` is the number of levels to compute.
Named parameters `start`, `stop`, and `length` are dummy.
"""
function energy_levels(V::MorsePotential, n::Integer;
                       start=0, stop=1, length=1, # Dummy args
                       ħ=ħ, m=m)
    ωₑ   = harmonic_frequency(V, m)
    ω(v) = ħ*ωₑ*(v + 1/2)       # Harmonic eigenlevels, v = 0, 1, ...

    vmax = trunc(Int, √(2V.D*m)/V.b/ħ - 1//2)

    if (n-1) <= vmax
        [ω(k) - ω(k)^2/4V.D for k=0:(n-1)]
    else
        error("Morse potential does not support $(n+1) bounded levels")
    end
end

"""
    analytic_trajectory(p₀, x₀, m, V::MorsePotential, t)

Analytic solution for the Morse potential. We implemented the solution given in
Ref.[^1], Eqs. 36 and 37. However, it seems that this solution is *not* exact;
to make it exact one should replace 1/(4s²-1) factors with just 1/4s². Ref.[^1]
is not the only paper that shows analytic solution for the Morse oscillator; the
earlier study on this subject is Ref.[^2]. We did not test solution from there.

[^1]: Eric M. Heatwole and Oleg V. Prezhdo. "Analytic dynamics of the Morse
      oscillator derived by semiclassical closures," J. Chem. Phys. 130, 244111
      (2009)

[^2]: W. C. DeMarcus. "Classical motion of a Morse oscillator," Am. J. Phys. 46,
      733-734 (1978)
"""
function analytic_trajectory(p₀::Number, x₀::Number,
                             m::Number,  V::MorsePotential, t)
    # Map Morse parameters to more conventional names
    D, xₑ, β = V.D, V.a, V.b

    # Initial energy
    E₀ = p₀^2/2m + V(x₀)

    # Various Morse mappings
    ν  = √(8m*D)/β
    s  = √(-2m*E₀)/β
    @inline y(x) = ν*exp(-β*(x - xₑ))
    y₀ = y(x₀)
    ω  = harmonic_frequency(V, m)
    # @assert ω == ħ*β^2*ν/2m
    γ  = 2s*ω/ν
    γt = γ*t
    βs = β*s

    u = zeros(2)
    den  =  p₀/βs*sin(γt) + cos(γt) +   y₀*ν/4s^2*(1-cos(γt))
    u[1] = (p₀*cos(γt) - βs*sin(γt) + β*y₀*ν/4s*sin(γt))/den # p(t)
    u[2] = log(den*ν/y₀)/β + xₑ                              # x(t)

    return u
end

end # module MorseModel
