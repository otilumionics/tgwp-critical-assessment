module QuarticModel
using Polynomials
using Main: Potential

# Module does not export, only adds methods

# Atomic units here and below
const m = 4668.0                # Reduced mass, similar to ammonia (not used)

struct QuarticPotential <: Potential
    params::Vector{Float64}

    function QuarticPotential(params::Vector{Float64})
        length(params) == 5 || error("Quartic potential is determined by order-4 polynomial")
        new(params)
    end
end

#=
    Symmetric quartic potential, derived to match ammonia potential --
    positions of minima & barrier height
=#
QuarticPotential() =            # Argument-less constructor with defaults
    QuarticPotential([
        0.05684,                # a0
        0.0,                    # a1*x
       -0.035132518834223594,   # a2*x²
        0.0,                    # a3*x³
        0.03353434564701764     # a4*x⁴
    ])

# Potential & its derivatives
function (p::QuarticPotential)(x::Number; der::Int=0)
    # @debug p

    poly = Polynomial(p.params)
    return derivative(poly, der)(x)
end

# NB: There is no "exact" definition for U: it coincides with a Taylor approx at 4th order

end # module QuarticModel
