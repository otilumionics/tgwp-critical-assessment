abstract type Potential end

# All subtypes of `Potential` (unless trivial with no fields) must provide
# parameter-less constructors
(::Type{<:Potential})() = error("Please, provide parameter-less constructor")

# Concrete subtypes of Potential are functors with the following signature
(::Potential)(x::Number; der::Int=0) = error("No concrete potential provided")

"""
A wrapper class to convert any univariative function `f(x)` and its first
derivative into an instance of `Potential`, see [`findmin`](@ref)
"""
struct GenericPotential <: Potential end

# Extended potential approximations
abstract type UApproximation end

struct Exact{T<:Potential} <:UApproximation
    p::T
    m::Float64
end

# To use with the argument-less constructor holding default values
Exact{T}(m) where T<:Potential = Exact{T}(T(), m)

struct GaussHermite{T<:Potential, N} <:UApproximation
    p::T
    m::Float64
end

# To use with the argument-less constructor holding default values
GaussHermite{T, N}(m) where {T<:Potential, N} = GaussHermite{T,N}(T(), m)

struct Taylor{T<:Potential, N} <:UApproximation
    p::T
    m::Float64
end

# To use with the argument-less constructor holding default values
Taylor{T, N}(m) where {T<:Potential, N} = Taylor{T,N}(T(), m)

# Signature for the extended potential U(q,w) function
U(a::UApproximation, q::Number, w::Number; der=(q=0,w=0)) = NaN
