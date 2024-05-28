"""
    Gaussian (coherent state) products, overlap and norm integrals
"""
module GaussianToolbox

export gwp, gaussian_overlap, gaussian_norm2

const ħ = 1                     # Plank's constant

function A_def(u, w)::Complex
    μ = w^2
    α = 2w*u/ħ

    return ħ/4μ*(1im + α)
end

r_def(A, p, q)::Complex = q - p/2A

"""
    gwp(x::Real, p::Real, q::Real, u::Real, w::Real, λ::Real)

Gaussian wave packet in a slightly generalized parametrization of Arickx et al.,
Chem. Phys. Lett. 128, 310 (1986).

    G = exp(-1/4w^2*(1-2iuw/ħ)*(x-q)^2 + ip/ħ*(x-q) + i/ħλ)

"""
gwp(x::Real,
    p::Real, q::Real, u::Real, w::Real, λ::Real)::Complex =
    gwp(x, A_def(u, w), p, q, λ)

function gwp(x::Real,
             A::Complex, p::Real, q::Real, λ::Real)::Complex
    r = r_def(A, p, q)

    return _gwp(x, A, r, λ)*exp(-1im/ħ*p^2/4A)
end

# Only for the sake of comparison
gwp_alt(x::Real,
        A::Complex, p::Real, q::Real, λ::Real)::Complex =
            exp(1im/ħ*(A*(x - q)^2 + p*(x - q) + λ))

_gwp(x::Real, A::Complex, r::Complex, λ::Real)::Complex =
    exp(1im/ħ*(A*(x - r)^2 + λ))

function gaussian_product(x::Real,
                          A0::Complex, p0::Real, q0::Real, λ0::Real,
                          A::Complex,  p::Real,  q::Real , λ::Real)::Complex
    r0 = r_def(A0, p0, q0)
    r  = r_def(A,  p,  q)

    return exp(1im/ħ*(p0^2/4conj(A0) - p^2/4A))*gaussian_product(x,
                                                                 A0, r0, λ0,
                                                                 A,  r,  λ)
end

function gaussian_product(x::Real,
                          A0::Complex, r0::Complex, λ0::Real,
                          A::Complex,  r::Complex,  λ::Real)::Complex
    ξ = A - conj(A0)
    b = (A*r - conj(A0*r0))/ξ
    c = conj(A0)*A*(r - conj(r0))^2/ξ

    return exp(1im/ħ*(ξ*(x - b)^2 - c + (λ - λ0)))
end

function gaussian_overlap(p0::Real, q0::Real, u0::Real, w0::Real, λ0::Real,
                          p::Real,  q::Real,  u::Real,  w::Real,  λ::Real)::Complex
    A0 = A_def(u0, w0)
    A =  A_def(u, w)

    return gaussian_overlap(A0, p0, q0, λ0, A, p, q, λ)
end

function gaussian_overlap(A0::Complex, p0::Real, q0::Real, λ0::Real,
                          A::Complex,   p::Real,  q::Real,  λ::Real)::Complex
    r0 = r_def(A0, p0, q0)
    r  = r_def(A,  p,  q)

    ξ = A - conj(A0)
    s = -1im/ħ*ξ
    c = conj(A0)*A*(r - conj(r0))^2/ξ

    return √(pi/s)*exp(1im/ħ*(p0^2/4conj(A0) - p^2/4A))* exp(1im/ħ*((λ - λ0) - c))
end

function gaussian_norm2(p::Real,  q::Real,  u::Real,  w::Real, λ::Real)::Real
    return √(2pi)*w
end

end # module GaussianToolbox
