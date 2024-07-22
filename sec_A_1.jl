figure = mkfigure()             # Restart figure numbering
section("Width dynamics in the ammonia model, Sec. A1")
let Vamm = AmmoniaModel.AmmoniaPotential(), m=AmmoniaModel.m

    def_full = Exact{AmmoniaModel.AmmoniaPotential}(m)
    def_u1   = Taylor{AmmoniaModel.AmmoniaPotential, 2}(m)

    U(q, w; kwarg...)   = Main.U(def_full, q, w; kwarg...)
    U_1(q, w; kwarg...) = Main.U(def_u1, q, w; kwarg...)

    Dynamics.setglobals!(V = Vamm, U = U, U_1 = U_1, m = m)
    # Solver parameters
    solver_params = (integrator = Tsit5(),
                     abstol     = 1.0e-12,
                     reltol     = 1.0e-12)
    @info "Solver parameters: " solver_params
    println()

    # Initial values, in a.u.
    inicn = (p0 = 0.0, q0 = -0.900,
             u0 = 0.0, w0 =  0.150,
             λ0 = 0.0)
    @info "Initial values: " inicn
    println()
    @debug("Potentials values at initial the configuration: ",
           Vamm(inicn.q0),
           U(inicn.q0, inicn.w0),
           U_1(inicn.q0, inicn.w0))

    # Time grid
    nt = 2^14
    t0 = 0.0
    dt = 2*2.015
    t_max = t0 + (nt-1)*dt
    tspan = (0.0, t_max)
    @info begin
        """The number of t steps, nt = $(nt)
           Time grid parameters:
        """
        end tspan dt
    println()

    println("Running dynamics...")
    sol_hq, sol_hlr =
        map([Dynamics.tdvp!,    # TDVP-based aka hemiquantal dynamics
             Dynamics.heller!]) do dyn_kind
                 @time solve(ODEProblem(dyn_kind, collect(inicn), tspan),
                             solver_params.integrator;
                             saveat=dt,
                             abstol=solver_params.abstol,
                             reltol=solver_params.reltol,
                             maxiters=10_000_000)
             end
    println()
    @assert sol_hq.t == sol_hlr.t

    @time sol_sc = solve(
        ODEProblem(Dynamics.ext_semiclassical!, collect(inicn), tspan),
        solver_params.integrator;
        saveat=dt,
        abstol=solver_params.abstol,
        reltol=solver_params.reltol)

    @debug("length of time arrays: ",
           length(sol_hq.t),
           length(sol_sc.t),
           length(sol_hlr.t))

    # Saving trajectory data into file. Need to handle a problem that
    # semiclassical dynamics may end prematurely
    ammonia_trajectory_sam_data =
        joinpath(datasave_path, "ammonia_trajectory_sam.dat")
    open(ammonia_trajectory_sam_data, "w") do io
        println(io, """
#      t                         p(t)                                    q(t)                                    u(t)                                    w(t)                                    λ(t)
#------------ --------------------------------------  --------------------------------------  --------------------------------------  --------------------------------------  --------------------------------------
#                 TDVP         Ext.sc       Heller        TDVP         Ext.sc       Heller        TDVP         Ext.sc       Heller        TDVP         Ext.sc       Heller        TDVP         Ext.sc       Heller """
                )
        for (i, t) in enumerate(sol_hq.t)
            @printf(io, " %12.5f %12.7f %12.7f %12.7f  %12.7f %12.7f %12.7f  %12.7f %12.7f %12.7f  %12.7f %12.7f %12.7f  %12.7f %12.7f %12.7f\n",
                    t,
                    sol_hq[1, i], get(view(sol_sc, 1, :), i, NaN), sol_hlr[1, i], # p(t)
                    sol_hq[2, i], get(view(sol_sc, 2, :), i, NaN), sol_hlr[2, i], # q(t)
                    sol_hq[3, i], get(view(sol_sc, 3, :), i, NaN), sol_hlr[3, i], # u(t)
                    sol_hq[4, i], get(view(sol_sc, 4, :), i, NaN), sol_hlr[4, i], # w(t)
                    sol_hq[5, i], get(view(sol_sc, 5, :), i, NaN), sol_hlr[5, i]) # λ(t)
        end
    end

    @debug begin
        io = IOBuffer()
        println(io, "SAM in ammonia: all components of the solution")
        for (k, varname) in enumerate(["p" "q" "u" "w"])
            title = "$(varname)-dynamics"
            println(io, "Figure A$(figure()): ", title)
            # @assert (sol_hq.t == sol_sc.t) && (sol_hq.t == sol_hlr.t)
            plot(sol_hq.t,
                 hcat(sol_hq[k, :],
                      sol_hlr[k, :]),
                 title=title,
                 labels=["TDVP" "Heller"],
                 xlabel="t, a.u",
                 ylabel="$(varname)(t)",
                 reuse=false)

            plot!(sol_sc.t,
                  sol_sc[k, :],
                  title=title,
                  label="Ext. Semiclassical",
                  xlabel="t, a.u",
                  ylabel="$(varname)(t)",
                  reuse=false)
            gui()
        end
        String(take!(io))
    end

    # Fig A1: width w(t) dynamics
    fig_A1_file      = joinpath(figures_path, "fig_A1.eps")
    gnuplot_script = """
#!/usr/bin/gnuplot
set terminal postscript eps enhanced color \\
    font "Helvetica" 24 lw 2.0 dl 1.0 \\
    size 1.2*5, 1.2*3.5

set key top left Left vertical reverse enhanced autotitles nobox \\
spacing 1.00 width -1.5 height 0.5

unset label

set xrange [ 0.0 : 15000.0 ]
set xtics 2500
set mxtics 5

set xlabel "Time, a.u."

set yrange [ 0.0 : 1.0 ]
set ytics 0.2
set mytics 2
set ylabel "width, bohr"

p_pos=1
q_pos=2
u_pos=3
w_pos=4
l_pos=5

sp=w_pos                        # Showing w(t)
    plot '$(ammonia_trajectory_sam_data)' u 1:(column((sp-1)*3+2)) w l lw 3 dt 1 lc "blue"   t "TDVP     ", \\
     ''                                 u 1:(column((sp-1)*3+3)) w l lw 3 dt 6 lc "red"    t "Ext. semicl", \\
     ''                                 u 1:(column((sp-1)*3+4)) w l lw 2 dt 3 lc "black"  t "Heller     "
"""
    run(pipeline(IOBuffer(gnuplot_script), `/usr/bin/env gnuplot`, fig_A1_file))

end # let
println()
