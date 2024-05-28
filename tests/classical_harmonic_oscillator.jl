section("Classical mechanics of the harmonic oscillator")
# Hamilton function: H = p²/2m + kx²/2 = p²/2m + mω²x²/2
let m = 2.0, k = 3.0
    ω = √(k/m)                  # Harmonic frequency

    # Initial conditions
    p₀, x₀ = 0.0, 1.0
    tspan  = (0.0, 6.0)

    # Analytic solution
    function analytic_trajectory(u₀, params, t)
        p₀, x₀ = u₀[1], u₀[2]

        u = zeros(2)
        u[1] = -x₀*m*ω*sin(ω*t) + p₀ *cos(ω*t)      # p(t)
        u[2] =  x₀*cos(ω*t)     + p₀/(m*ω)*sin(ω*t) # x(t)

        return u
    end

    Dynamics.setglobals!(V = GenericPotential(), m = m)

    function (p::GenericPotential)(x::Number; der::Int=0)
        der == 0 && return k*x^2/2
        der == 1 && return k*x

        return nothing
    end

    @info "Initial energy: " Dynamics.E_cl(p₀, x₀)

    # Setting up the ODE
    f       = ODEFunction(Dynamics.classical!; analytic=analytic_trajectory)
    problem = ODEProblem(f, [p₀, x₀], tspan)

    println("Running dynamics...")
    @time solution_cl = solve(problem, Tsit5();
                              saveat = 0.1,
                              abstol = 1.0e-8,
                              reltol = 1.0e-8)
    println()

    println("Plotting...")
    plot(solution_cl.t,
         hcat(solution_cl[1, :],
              [u[1] for u in solution_cl.u_analytic],
              solution_cl[2, :],
              [u[2] for u in solution_cl.u_analytic]),
         title="Harmonic oscillator solution",
         lw = [1 2 1 2],
         linestyles = [:solid :dash :solid :dash],
         labels = ["p(t), numeric" "p(t), analytic" "q(t), numeric" "q(t), analytic"],
         reuse=false)
    gui()

    plot(solution_cl.t,
         solution_cl[2, :] .- [u[2] for u in solution_cl.u_analytic],
         title="Errors in numerical solution for x(t)",
         labels=nothing,
         reuse=false)
    gui()

    Dynamics.resetglobals!()
end
println()
