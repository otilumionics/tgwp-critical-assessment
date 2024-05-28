module WavepacketTools
using Main: Potential, UApproximation, U, findmin

"""
    struct UWrapper{T<:UApproximation, Q<:Number} <: Potential
        a::T
        x::Q
    end

Wrapper for converting the extended potential `U(x, w)` into a univariative
function of `w` (instance of `Potential`) that parametrically depends on `x`
"""
struct UWrapper{T<:UApproximation, Q<:Number} <: Potential
    a::T
    x::Q
end
(obj::UWrapper)(w::Number; der::Int=0) = U(obj.a, obj.x, w; der=(q=0,w=der))

"""
    minimal_uncertainty(x::Number)

Find the minimum-uncertainty width of a Gaussian wave packet centered at `x`
using the extended potential approximation `a`
"""
function minimal_uncertainty(a::UApproximation, x::Number; guess)
    findmin(UWrapper(a, x), guess=guess)
end

end # module WavepacketTools
