section("Tunneling dynamics in ammonia, Sec. 4.4")
let Vamm = AmmoniaModel.AmmoniaPotential(), m=AmmoniaModel.m

    def_full = Exact{AmmoniaModel.AmmoniaPotential}(m)
    def_u1   = Taylor{AmmoniaModel.AmmoniaPotential, 2}(m)

    U(q, w; kwarg...)   = Main.U(def_full, q, w; kwarg...)
    U_1(q, w; kwarg...) = Main.U(def_u1, q, w; kwarg...)

    # Minimal uncertainty wavepacket
    xrange = LinRange(-1.25, 1.25, 101)
    U_min_uncert = similar(xrange)
    w_min_uncert = similar(xrange)
    with_logger(NullLogger()) do # run silently
        for (i, x) in enumerate(xrange)
            U_min_uncert[i], w_min_uncert[i] =
                WavepacketTools.minimal_uncertainty(def_full, x; guess=0.1)
        end
    end

    ammonia_pots_data = joinpath(datasave_path, "ammonia_pots.dat")
    open(ammonia_pots_data, "w") do io
        println(io, "#  q, bohr     w_min_uncert, bohr   V, Ha          U_min_incert, Ha")
        for (k, q) in enumerate(xrange)
            @printf(io, "%9.4f   %12.7f         %12.7f   %12.7f\n",
                    q,
                    w_min_uncert[k],
                    Vamm(q),
                    U_min_uncert[k])
        end
    end

    if length(xrange) == 101    # Magic constant
        Umu_local_max, _ = findmax(U_min_uncert[40:62])
    else
        error("finding local maxima of U_min_uncert assumed 101-point grid")
    end
    nlvl   = 6
    levels = energy_levels(Vamm, nlvl, m=m, start=-5.0, stop=5.0, length=99)
    @debug "energy levels" levels

    # Fig 5: V, minimal uncertainty potential, and tunneling region.
    fig5_file      = joinpath(figures_path, "fig5.pdf")
    # Gnuplot: terminal here is set to `pdfcairo` to support transparency
    gnuplot_script = """
#!/usr/bin/gnuplot
set terminal pdfcairo enhanced color \\
    font "Liberation Sans, 18" lw 2.0 dl 1.2 \\
    size 1.2*5, 1.2*3.5

set key top right Left reverse enhanced autotitles nobox \\
spacing 1.00 width -2 height 1.5

set xrange [ -1.250 : 1.250 ]
set xtics auto
set mxtics 2

set xlabel "Coordinate, bohr"

# Take advantage of Julia representation of arrays via string interpolation
array levels[$nlvl] = $(levels)
f(i,x) = x < -1.15 || x > 1.15 ? NaN : levels[i]

set yrange [ 0.045 : 0.070]
set ytics 0.005
set mytics 5
set ylabel "Energy, E_h"

plot     '$(ammonia_pots_data)' u 1:3 w l lw 1.5 dt 1 lc "blue"   t "V_A     ", \\
         ''                     u 1:4 w l lw 1.5 dt (3,5) lc "red"    t "min. uncertainty", \\
for [i=1:$(nlvl)] ''            u 1:(f(i,\$1)) w l lc -1 dt (15,15) not, \\
         $(Vamm(0.0))                         w filledcurves above y1=$(Umu_local_max) lc "black" fs transparent solid 0.3 noborder not
"""
    run(pipeline(IOBuffer(gnuplot_script), `/usr/bin/env gnuplot`, fig5_file))

    Dynamics.setglobals!(V = Vamm, U = U, U_1 = U_1, m = m)
    # Solver parameters
    solver_params = (integrator = Tsit5(),
                     abstol     = 1.0e-12,
                     reltol     = 1.0e-12)
    @info "Solver parameters: " solver_params
    println()

    # Initial values, in a.u.
    inicn = (p0 = 0.0, q0 = -1.000,
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
    ammonia_trajectory_data =
        joinpath(datasave_path, "ammonia_trajectory.dat")
    open(ammonia_trajectory_data, "w") do io
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
        println(io, "Tunneling in ammonia: all components of the solution")
        for (k, varname) in enumerate(["p" "q" "u" "w"])
            title = "$(varname)-dynamics"
            println(io, "Figure $(figure()): ", title)
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

    # Fig 6. Coordinate dynamics in tunneling regime
    fig6_file      = joinpath(figures_path, "fig6.eps")
    gnuplot_script = """
#!/usr/bin/gnuplot
set terminal postscript eps enhanced color \\
    font "Helvetica" 24 lw 2.0 dl 1.0 \\
    size 1.2*5, 1.2*3.5

set key top right Left horizontal reverse enhanced autotitles nobox \\
spacing 1.00 width -2.5 height 0.5

unset label

set xrange [ 0.0 : 30000.0]
set xtics auto
set mxtics 10

set xlabel "Time, a.u."

set yrange [ -1.2 : 1.2 ]
set ytics 0.5
set mytics 5
set ylabel "Coordinate, bohr"

set xzeroaxis

p_pos=1
q_pos=2
u_pos=3
w_pos=4
l_pos=5

sp=q_pos                        # Showing position
    plot '$(ammonia_trajectory_data)' u 1:(column((sp-1)*3+2)) w l lw 3 dt 1 lc "blue"   t "TDVP     ", \\
     ''                               u 1:(column((sp-1)*3+3)) w l lw 3 dt 6 lc "red"    t "Ext. semicl", \\
     ''                               u 1:(column((sp-1)*3+4)) w l lw 2 dt 3 lc "black"  t "Heller     "
"""
    run(pipeline(IOBuffer(gnuplot_script), `/usr/bin/env gnuplot`, fig6_file))

end # let
println()
