section("Classical mechanics of the Morse oscillator")
let V = MorseModel.MorsePotential(), m = MorseModel.m
    # Initial conditions
    p₀, x₀ = 0.0, 1.400         # au, bohr, similar to Fig. 7
    tspan  = (0.0, 4000.0)

    analytic_trajectory(u₀, params, t) =
        MorseModel.analytic_trajectory(u₀[1], u₀[2], params[1], params[2], t)

    Dynamics.setglobals!(V = V, m = m)
    @info("Components of the initial energy: ",
          p₀^2/2m,
          V(x₀),
          Dynamics.E_cl(p₀, x₀))

    # Setting up the ODE
    f       = ODEFunction(Dynamics.classical!; analytic=analytic_trajectory)
    problem = ODEProblem(f, [p₀, x₀], tspan, (m, V))

    println("Running dynamics...")
    @time solution_cl = solve(problem, Tsit5();
                                     saveat = 0.1,
                                     abstol = 1.0e-8,
                              reltol = 1.0e-8)
    println()

    println("Plotting...")
    plot(solution_cl.t,
         hcat(solution_cl[2, :],
              [u[2] for u in solution_cl.u_analytic]),
         title="Morse oscillator solution",
         lw = [1 2],
         linestyles = [:solid :dash],
         ylims = (1.3, 3.5),    # Similar to Fig. 7.
         labels = ["q(t), numeric" "q(t), analytic"],
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
