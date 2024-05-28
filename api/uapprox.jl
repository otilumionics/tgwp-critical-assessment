"""
    Numerical approximations for the extended potential U(q, w)
"""
module ExtendedPotential
using FastGaussQuadrature: gausshermite
using Main: Potential, GaussHermite, Taylor
import Main: U

const ħ = 1

"""
    U(a::GaussHermite{T, N}, q::Number, w::Number; der=(q=0,w=0)) where {T<:Potential, N}

Approximation of U(q, w) using Gauss-Hermite quadratures. Implements Eq.~(29) of
the manuscript.
"""
function U(a::GaussHermite{T, N}, q::Number, w::Number;
           der=(q=0,w=0)) where {T<:Potential, N}
    # `a` stands for an "approximation"
    m = a.m
    V = a.p
    (xi, wi) = gausshermite(N)

    args = sqrt(2)*w*xi .+ q
    if der == (q=0, w=0)
        return ħ^2/m/w^2/8 + 1/√π*(wi'*V.(args))
    elseif der == (q=1, w=0)
        # ∂U/∂q
        return 1/√π*(wi'*V.(args, der=1))
    elseif der == (q=0, w=1)
        # ∂U/∂w; FIXME: ħ in normalization constant?
        return -ħ^2/m/w^3/4 + √(2/pi)*(wi'*(xi.*V.(args, der=1)))
    else
        error("Unknown/unimplemented partial derivative: der=$(der)")
    end

    return nothing              # Should not be reached!
end

"""
    U(a::Taylor{T, N}, q::Number, w::Number; der=(q=0,w=0)) where {T<:Potential, N}

Approximation of U using Taylor series, see Eq.~(17) of the manuscript.
"""
function U(a::Taylor{T, N}, q::Number, w::Number;
           der=(q=0,w=0)) where {T<:Potential, N}
    # `a` stands for an "approximation"
    m = a.m
    V = a.p

    if der == (q=0, w=0)
        U = V(q) + ħ^2/m/w^2/8
        for n=1:(N÷2)
           # U += V(q, der=2n)*(w^2/2)^n/factorial(n)
            U += V(q, der=2n)*prod(w^2/2/k for k=1:n)
        end
        return U
    elseif der == (q=1, w=0)
        # ∂U/∂q
        dU_q = V(q, der=1)
        for n=1:(N÷2)
            # dU_q += V(q, der=2n+1)*(w^2/2)^n/factorial(n)
            dU_q += V(q, der=2n+1)*prod(w^2/2/k for k=1:n)
        end
        return dU_q
    elseif der == (q=0, w=1)
        # ∂U/∂w
        dU_w = -ħ^2/m/w^3/4
        for n=1:(N÷2)
            dU_w += V(q, der=2n)*w*(w^2/2)^(n-1)/factorial(n-1)
        end
        return dU_w
    else
        error("Unknown/unimplemented partial derivative: der=$(der)")
    end

    return nothing              # Should not be reached
end

end # module ExtendedPotential
