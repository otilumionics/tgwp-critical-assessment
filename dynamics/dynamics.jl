module Dynamics

"""
    globals = Dict{Symbol, Any}(...)

Global parameters available in the module [`Dynamics`](@ref). They are (with
default values):

    :m => nothing, # Particle's mass

    :V => nothing, # an instance of `Potential`

    :U => nothing, # extended potential U(q,w)

    :U_1 => nothing, # First Taylor approximation to U(q,w)
"""
globals = Dict{Symbol, Any}(
    :m   => nothing,
    :V   => nothing,
    :U   => nothing,
    :U_1 => nothing,
)

function setglobals!(; kwarg...)
    for (key, value) in kwarg
        if haskey(globals, key)
            globals[key] = value
        else
            error("Attempt to set unknown variable: `$(repr(key))`")
        end
    end

    return nothing
end

function resetglobals!()
    for (key, value) in globals
        globals[key] = nothing
    end

    return nothing
end

## Energy functions
"""
    E_cl(p::Number, x::Number)

Classical energy E = p^2/2m + V(x) of a particle with a momentum `p` and
coordinate` x. Particle's mass `m` and explicit form of a potential `V` is taken
from a dictionary [`globals`](@ref).
"""
function E_cl(p::Number, x::Number)
    # Classical energy
    V, m = globals[:V], globals[:m]

    return p^2/2m + V(x)
end

function E_hq(p::Number, q::Number, u::Number, w::Number)
    # Hemiquantal energy
    U, m = globals[:U], globals[:m]

    return p^2/2m + u^2/2m + U(q, w)
end

function E_sc(p::Number, q::Number, u::Number, w::Number)
    # Semiclassical energy
    U_1, m = globals[:U_1], globals[:m]

    return p^2/2m + u^2/2m + U_1(q, w)
end

## EOM for dynamical models
function classical!(dy, y, params, t)
    p, q = y
    V, m = globals[:V], globals[:m]

    dy[1] = -V(q, der=1)        # ṗ = -∂V/∂q
    dy[2] =  p/m                # q̇

    return nothing
end

function tdvp!(dy, y, params, t)
    p, q, u, w, λ = y
    U, m = globals[:U], globals[:m]

    # Hamilton equations for a (p,q) pair
    # y[1] = p, y[2] = q
    dy[1] = -U(q, w, der=(q=1,w=0)) # ṗ = -∂U/∂q
    dy[2] =  p/m                    # q̇

    # Hamilton equations for a (u,w) pair
    # y[3] = u, y[4] = w
    dy[3] = -U(q, w, der=(q=0,w=1)) # u̇ = -∂U/∂w
    dy[4] =  u/m                    # ẇ

    # Phase dynamics
    # y[5] = λ
    dy[5] = (dy[4]*u - dy[3]*w)/2 + p*dy[2] - E_hq(p, q, u, w) # Action

    return nothing
end

function ext_semiclassical!(dy, y, params, t)
    p, q, u, w, λ = y
    U_1, m = globals[:U_1], globals[:m]

    # Hamilton equations for a (p,q) pair
    # y[1] = p, y[2] = q
    dy[1] = -U_1(q, w, der=(q=1,w=0)) # ṗ = -∂U¹/∂q
    dy[2] =  p/m                      # q̇

    # Hamilton equations for a (u,w) pair
    # y[3] = u, y[4] = w
    dy[3] = -U_1(q, w, der=(q=0,w=1)) # u̇ = -∂U¹/∂w
    dy[4] =  u/m                      # ẇ

    # Phase dynamics
    # y[5] = λ
    dy[5] = (dy[4]*u - dy[3]*w)/2 + p*dy[2] - E_sc(p, q, u, w) # Action

    return nothing
end

function heller!(dy, y, params, t)
    p, q, u, w, λ = y
    V, U_1, m = globals[:V], globals[:U_1], globals[:m]

    # Hamilton equations for a (p,q) pair use classical V
    # y[1] = p, y[2] = q
    dy[1] = -V(q, der=1)        # ṗ = -V'(q)
    dy[2] =  p/m                # q̇

    # Hamilton equations for a (u,w) pair use U_1
    # y[3] = u, y[4] = w
    dy[3] = -U_1(q, w, der=(q=0,w=1)) # u̇ = -∂U_1/∂w
    dy[4] =  u/m                      # ẇ

    # Phase dynamics, uses U_1(q, w)
    # y[5] = λ
    dy[5] = (dy[4]*u - dy[3]*w)/2 + p*dy[2] - E_sc(p, q, u, w) # Action

    return nothing
end

end # module Dynamics
