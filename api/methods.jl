using NLopt
using LinearAlgebra: norm, Symmetric, eigvals!, eigen
import Base: findmin

"""
    findmin(V::Potential; guess, optimizer=:LD_LBFGS)

Find a minimum of a (univariative) potential `V` numerically starting from a
`guess`. Return a pair `(Vmin, xmin)`

The default optimizer is NLopt.:LD_LBFGS (a form of the BFGS scheme), which
requires the first-order derivative of `V`. Any univariative function `f(x)`
with a known first derivative can be converted into an instance of `Potential`
using the parameter-less wrapper class `GenericPotential`. If, however, `f(x)`
has no known expression for the first derivative, one may employ automatic
differentiation via dual numbers as

```julia
using DualNumbers

function (::GenericPotential)(x::Number; der::Int=0)
    der == 0 && return f(x)
    der == 1 && return dualpart(f(dual(x, 1)))

    return nothing
end

julia> findmin(GenericPotential(), ...)
```

Some models may implement specialized versions. For example, if their minima are
known analytically.
"""
function findmin(V::Potential; guess, optimizer=:LD_LBFGS)
    function obj_f(x::Vector, grad::Vector)
        if length(grad) > 0
            grad[1] = V(x[1], der=1)
        end
        return V(x[1])
    end

    opt = Opt(optimizer, 1)
    min_objective!(opt, obj_f)

    minf, argmin, ret_status = optimize(opt, [guess])
    @debug "Optimizer returned status: $(repr(ret_status))"
    @debug "Minimum found:" guess argmin[1] minf

    return (minf, argmin[1])
end

"""
    findmin(a::UApproximation; guess::NamedTuple{(:q,:w)}, optimizer=:LD_LBFGS)

Find a minimum of the extended potential `U(q, w)` evaluated under approximation
`a` starting from a `guess`.
"""
function findmin(a::UApproximation;
                 guess::NamedTuple{(:q,:w)}, optimizer=:LD_LBFGS)
    function obj_f(x::Vector, grad::Vector)
        if length(grad) > 0
            grad[1] = U(a, x..., der=(q=1, w=0))
            grad[2] = U(a, x..., der=(q=0, w=1))
        end
        return U(a, x..., der=(q=0,w=0))
    end
    opt = Opt(:LD_LBFGS, 2)
    min_objective!(opt, obj_f)

    minf, argmin, status = optimize(opt, [guess...])
    @debug "Optimizer returned status: $(repr(status))"
    @debug "Minimum found:" guess argmin minf

    return (minf, (q=argmin[1], w=argmin[2]))
end

"""
    energy_levels(V, nstates::Integer; start, stop, length, ħ=1, m=1Da,
    extended=false)

Eigenvalues of the one-dimensional Schrödinger equation with the Hamiltonian Ĥ =
T̂ + V̂ via the Fourier grid Hamiltonian (FGH) method [^1], [^2]. We follow a more
recent work [^1], where matrix elements of the kinetic energy are computed in
O(N) process, but occasionally refer to Ref.[^2] as it provides better
explanation of the method.

`V` is function, `nstates` is the number of levels requested, `start`, `stop`,
and `length` are parameters of the internal grid, which are similar (but not
identical!) to those of the function [`range`](@ref). `ħ` is the numerical value
of the Plank constant, and `m` is a reduced mass, both determine the energy
scale for the returned energy levels. If `ħ=1` and `m` is given in electron mass
(mₑ) then levels are in Hartree (Eₕ) provided that `V` accepts an argument in
bohrs. If `extended=true` a vector of eigenvalues, a matrix eigenvectors and
a grid specification is returned.

Some models may implement specialized versions (with `V` as subtypes of
`Potential`) -- for example, if eigenvalues are known analytically.

[^1]: J. Stare and G. G. Balint-Kurti, “Fourier grid Hamiltonian method for
      solving the vibrational Schrödinger equation in internal coordinates:
      theory and test applications,” J. Phys. Chem A 107, 7204-7214 (2003).

[^2]: C. C. Marston and G. G. Balint‐Kurti. “The Fourier grid Hamiltonian method
      for bound state eigenvalues and eigenfunctions,” J. Chem. Phys. 91,
      3571-3576 (1989).
"""
function energy_levels(V, nstates::Integer;
                       start, stop, length, # Internal grid parameters
                       ħ=1, m=1Constants.Da2me,
                       extended=false)
    # Convenient remapping to align with notation in Ref.[1] & [2]
    L, N = (stop - start), length

    # As emphasized in Ref.[2], a grid must have the odd number of points
    isodd(N) || error("even number of points in a grid")
    n  = (N-1) ÷ 2              # Eq. 17 of Ref.[2]

    # NOTE: Explanation given in Ref.[1] after Eq.(9) or in Ref.[2] near
    # Eq.(36) about the number of grid points and its relation to the length and
    # step size is INCORRECT! The correct description follows.
    #
    # Spatial grid: from `x0 = start` to x_{n-1}. `x_n = x0 + L = stop` is an
    # implicit point because the Fourier basis is periodic, x0 ≡ x_n. Hence, the
    # total number of grid points (including the implicit one) is (N+1) and the
    # step size is:
    Δx = L/N
    x  = collect(start .+ (0:(N-1))*Δx) # Explicit points, Eq.(36) of Ref.[2]

    # Equidistant k (wave-number) grid with a step size
    Δk = 2pi/L                  # Eq.(37) of Ref.[2]
    ω = pi/N

    # T is a symmetric
    # [Toeplitz matrix](https://en.wikipedia.org/wiki/Toeplitz_matrix) in the
    # FGH method; Form its first column by Eq.(13) of Ref.[1]
    tfc = Vector{Float64}(undef, N)
    tfc[1] = n*(n+1)/3          # i=j
    for m=1:(N-1)               # m = (i-j), sub-diagonals
        ωm = ω*m
        tfc[m+1] =
            inv(2N*sin(ωm)^2) * (
                (-1)^m*(n+1)*cos(ωm)  +
                (n+1)*cos(2*(n+1)*ωm) -
                sin(2*(n+1)*ωm)*cot(ωm))
    end
    tfc .*= ħ^2/2m*Δk^2

    H = Matrix{Float64}(undef, N, N)
    for j=1:N, i=j:N
        H[i,j] = tfc[(i-j)+1]   # The lower triangle (uplo=:L) of T stored in H
    end

    # @debug begin
    #     T = collect(Symmetric(H, :L))

    #     # Test of the formulae from Ref.[2], equation numbers refer to that paper
    #     Talt = zeros(Float64, N, N)
    #     # k-space grid is (-n:n)Δk, but only a (+)half is used as T(-k) = T(k)
    #     let k = collect( (1:n)*Δk) # Eq. 37, a positive-k half
    #         # Kinetic energy is diagonal in k representation
    #         Tdiag =  ħ^2/2m*k.^2   # Eq. 22
    #         for j=1:N, i=1:N
    #             tij = 0.0
    #             for l=1:n
    #                 tij += cos(2ω*l*(i-j))*Tdiag[l] # First term in rhs of Eq.26
    #             end
    #             Talt[i,j] += 2/N*tij
    #         end
    #     end
    #     "Alternative formula for T from Ref.[2]"
    # end (T - Talt) norm(T - Talt)

    # Main diagonal of H is T[i,i] + V(x[i])
    for i=1:N
        H[i,i] += V(x[i])
    end

    if extended                 # Return eigenvalues, eigenvectors, and the grid
        energies, eigenvectors = eigen(Symmetric(H, :L), 1:nstates)
        return energies, eigenvectors, x
    end

    return eigvals!(Symmetric(H, :L), 1:nstates) # Don't care about preserving H
end
