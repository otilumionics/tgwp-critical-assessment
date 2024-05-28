section("Small-amplitude motion + spectra for the Morse potential, Sec. 4.3")
let Vᴹ = MorseModel.MorsePotential(), m = MorseModel.m

    def_full = Exact{MorseModel.MorsePotential}(m)
    def_u1   = Taylor{MorseModel.MorsePotential, 2}(m)

    U(q, w; kwarg...)   = Main.U(def_full, q, w; kwarg...)
    U_1(q, w; kwarg...) = Main.U(def_u1, q, w; kwarg...)

    Dynamics.setglobals!(V = Vᴹ, U = U, U_1 = U_1, m = m)

    # Solver parameters
    solver_params = (integrator = Tsit5(),
                     abstol     = 1.0e-8,
                     reltol     = 1.0e-8)
    @info "Solver parameters: " solver_params
    println()

    # Initial values, in a.u.
    inicn = (p0 = 0.0, q0 =  1.600,
             u0 = 0.0, w0 =  0.125,
             λ0 = 0.0)
    @info "Initial values: " inicn
    println()

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

    @info("Initial energy in various definitions",
          E_cl = Dynamics.E_cl(inicn.p0, inicn.q0),
          E_hq = Dynamics.E_hq(collect(inicn)[1:4]...),
          E_sc = Dynamics.E_sc(collect(inicn)[1:4]...))
    println()

    @debug("Effective potential values",
           Vᴹ.D,
           Vᴹ(inicn.q0),
           U(inicn.q0, inicn.w0),
           U_1(inicn.q0, inicn.w0))

    println("Running dynamics...")
    sol_hq, sol_sc, sol_hlr =
        map([Dynamics.tdvp!,    # TDVP-based aka "hemiquantal" dynamics
             Dynamics.ext_semiclassical!,
             Dynamics.heller!]) do dyn_kind
                 @time solve(ODEProblem(dyn_kind, collect(inicn), tspan),
                       solver_params.integrator;
                       saveat=dt,
                       abstol=solver_params.abstol,
                       reltol=solver_params.reltol)
             end
    println()

    @debug("Components of the TDVP solution at the end", sol_hq.t[end],
           sol_hq[1,end], sol_hq[2,end], # (p,q)
           sol_hq[3,end], sol_hq[4,end]) # (u,w)

    @debug begin
        io = IOBuffer()
        println(io, "Online plots of all components of the solution")
        for (k, varname) in enumerate(["p(t)", "q(t)", "u(t)", "w(t)"])
            npt = ceil(Int, 5300/dt)
            title = "$(varname)"
            println(io, "Drawing $(figure()): ", title)
            @assert (sol_hq.t == sol_sc.t) && (sol_hq.t == sol_hlr.t)
            plot(sol_hq.t[1:npt],
                 hcat(sol_hq[k,  1:npt],
                      sol_sc[k,  1:npt],
                      sol_hlr[k, 1:npt]),
                 title=title,
                 labels=["TDVP" "Ext. Semiclassical" "Heller"],
                 xlabel="t, a.u",
                 reuse=false)
            gui()
        end
        println(io, "Drawing $(figure()): λ(t)/t")
        plot(sol_hq.t,
             hcat(sol_hq[5,  :]./sol_hq.t,
                  sol_sc[5,  :]./sol_hq.t,
                  sol_hlr[5, :]./sol_hq.t),
             xlims = (0, 5300),
             title="λ(t)/t",
             labels=["TDVP" "Ext. Semiclassical" "Heller"],
             xlabel="t, a.u",
             reuse=false)
        gui()
        String(take!(io))
    end

    # Saving trajectory data into file
    morse_trajectory_sam_data =
        joinpath(datasave_path, "morse_trajectory_sam.dat")
    open(morse_trajectory_sam_data, "w") do io
        println(io, """
#      t                         p(t)                                    q(t)                                    u(t)                                    w(t)                                    λ(t)
#------------ --------------------------------------  --------------------------------------  --------------------------------------  --------------------------------------  --------------------------------------
#                 TDVP         Ext.sc       Heller        TDVP         Ext.sc       Heller        TDVP         Ext.sc       Heller        TDVP         Ext.sc       Heller        TDVP         Ext.sc       Heller """

                )
        for (i, t) in enumerate(sol_hq.t)
            @printf(io, " %12.5f %12.7f %12.7f %12.7f  %12.7f %12.7f %12.7f  %12.7f %12.7f %12.7f  %12.7f %12.7f %12.7f  %12.7f %12.7f %12.7f\n",
                    t,
                    sol_hq[1, i], sol_sc[1, i], sol_hlr[1, i], # p(t)
                    sol_hq[2, i], sol_sc[2, i], sol_hlr[2, i], # q(t)
                    sol_hq[3, i], sol_sc[3, i], sol_hlr[3, i], # u(t)
                    sol_hq[4, i], sol_sc[4, i], sol_hlr[4, i], # w(t)
                    sol_hq[5, i], sol_sc[5, i], sol_hlr[5, i]) # λ(t)
        end
    end

    # Fig 2: Linear momentum p(t) dynamics
    fig2_file      = joinpath(figures_path, "fig2.eps")
    gnuplot_script = """
#!/usr/bin/gnuplot
set terminal postscript eps enhanced color \\
    font "Helvetica" 24 lw 2.0 dl 1.0 \\
    size 1.2*5, 1.2*3.5

set key top right Left horizontal reverse enhanced autotitles nobox \\
spacing 1.00 width -2.5 height 0.5

unset label

set xrange [ 0.0 : 5300.0 ]
set xtics auto
set mxtics 10

set xlabel "Time, a.u."

set yrange [ -10.0 : 10.0 ]
set ytics 5
set mytics 5
set ylabel "p(t), a.u."

set xzeroaxis

p_pos=1
q_pos=2
u_pos=3
w_pos=4
l_pos=5

sp=p_pos                        # Showing momentum
    plot '$(morse_trajectory_sam_data)' u 1:(column((sp-1)*3+2)) w l lw 3 dt 1 lc "blue"   t "TDVP     ", \\
     ''                                 u 1:(column((sp-1)*3+3)) w l lw 3 dt 6 lc "red"    t "Ext. semicl", \\
     ''                                 u 1:(column((sp-1)*3+4)) w l lw 2 dt 3 lc "black"  t "Heller     "
"""
    run(pipeline(IOBuffer(gnuplot_script), `/usr/bin/env gnuplot`, fig2_file))

    # Fig 3: width w(t) dynamics
    fig3_file      = joinpath(figures_path, "fig3.eps")
    gnuplot_script = """
#!/usr/bin/gnuplot
set terminal postscript eps enhanced color \\
    font "Helvetica" 24 lw 2.0 dl 1.0 \\
    size 1.2*5, 1.2*3.5

set key top left Left vertical reverse enhanced autotitles nobox \\
spacing 1.00 width -1.5 height 0.5

unset label

set xrange [ 0.0 : 5000.0 ]
set xtics 1000
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
    plot '$(morse_trajectory_sam_data)' u 1:(column((sp-1)*3+2)) w l lw 3 dt 1 lc "blue"   t "TDVP     ", \\
     ''                                 u 1:(column((sp-1)*3+3)) w l lw 3 dt 6 lc "red"    t "Ext. semicl", \\
     ''                                 u 1:(column((sp-1)*3+4)) w l lw 2 dt 3 lc "black"  t "Heller     "
"""
    run(pipeline(IOBuffer(gnuplot_script), `/usr/bin/env gnuplot`, fig3_file))

    @debug begin
        io = IOBuffer()
        title = "Energy components for the TDVP-based dynamics"
        println(io, "Drawing $(figure()): ", title)
        plot(sol_hq.t,
             hcat(Dynamics.E_cl.(sol_hq[1, :], sol_hq[2, :]),
                  Dynamics.E_sc.(sol_hq[1, :], sol_hq[2, :],
                                 sol_hq[3, :], sol_hq[4, :]),
                  Dynamics.E_hq.(sol_hq[1, :], sol_hq[2, :],
                                 sol_hq[3, :], sol_hq[4, :])),
             xlims=(0, 5000),
             ylims=(-Vᴹ.D, -0.1),
             title=title,
             labels=["E_cl"  "E_sc" "E_hq"],
             xlabel="t, a.u.",
             reuse=false)
        gui()

        title = "Energy components for Heller's dynamics"
        println(io, "Drawing $(figure()): ", title)
        plot(sol_hlr.t,
             hcat(Dynamics.E_cl.(sol_hlr[1, :], sol_hlr[2, :]),
                  Dynamics.E_sc.(sol_hlr[1, :], sol_hlr[2, :],
                                 sol_hlr[3, :], sol_hlr[4, :])),
             xlims=(0, 5000),
             ylims=(-Vᴹ.D, -0.1),
             title=title,
             labels=["E_cl"  "E_sc"],
             xlabel="t, a.u.",
             reuse=false)
        gui()
        String(take!(io))
    end

    # Spectra simulations
    freqs_eh =
        Spectra.freqs(sol_hq.t, dt; scale=2pi, shift=Vᴹ.D)*Constants.Eh2invcm

    C_hq = [Spectra.autocorrelation(t, inicn..., vars...)
            for (t, vars) in zip(sol_hq.t, sol_hq.u)]
    spectrum_hq = Spectra.spectrum(C_hq)

    C_sc = [Spectra.autocorrelation(t, inicn..., vars...)
            for (t, vars) in zip(sol_sc.t, sol_sc.u)]
    spectrum_sc = Spectra.spectrum(C_sc)

    C_hlr = [Spectra.autocorrelation(t, inicn..., vars...)
     for (t, vars) in zip(sol_hlr.t, sol_hlr.u)]
    spectrum_hlr = Spectra.spectrum(C_hlr)

    # Dumping spectral information to file
    morse_spectra_sam_data = joinpath(datasave_path, "morse_spectra_sam.dat")
    open(morse_spectra_sam_data, "w") do io
        println(io, """
# ω                      Intensity
#               TDVP          Ext.semicl    Heller """
                )
        for (i, ω) in enumerate(freqs_eh)
            @printf(io, " %12.5f %12.7e %12.7e %12.7e\n",
                    ω,
                    abs.(spectrum_hq[i]),
                    abs.(spectrum_sc[i]),
                    abs.(spectrum_hlr[i]))
        end
    end

    # Fig4: Spectra via gnuplot
    levels         = energy_levels(Vᴹ, 4, m=m)*Constants.Eh2invcm
    fig4_file      = joinpath(figures_path, "fig4.eps")
    gnuplot_script = """
#!/usr/bin/gnuplot
set terminal postscript eps enhanced color \\
    font "Helvetica" 24 lw 2.0 dl 1.0 \\
    size 1.2*5, 1.2*3.5

set multiplot
set size 1.00, 1.00
set origin 0, 0

set key top left Left reverse enhanced autotitles nobox \\
spacing 1.0 samplen 4.0 width 1.0 height 0.50 opaque

# Vertical lines with eigenenergies
unset arrow
set arrow from $(levels[1]), graph 0.00 to $(levels[1]), graph 0.98 lw 2 lc -1 dt 3 back nohead
set arrow from $(levels[2]), graph 0.00 to $(levels[2]), graph 0.98 lw 2 lc -1 dt 3 back nohead
set arrow from $(levels[3]), graph 0.00 to $(levels[3]), graph 0.98 lw 2 lc -1 dt 3 back nohead
set arrow from $(levels[4]), graph 0.00 to $(levels[4]), graph 0.98 lw 2 lc -1 dt 3 back nohead

set xrange [ 0.0 : 15000.0 ]
set xtics 2500
set mxtics 5

set xlabel "Frequency, cm^{-1}"

set yrange [ 0.0 : 0.60 ]
set ytics 0.15
set mytics 3
set ylabel "Intensity, arbitrary units"

set object 1 rect from graph 0.485,0.505 to graph 0.935,0.952 fc "light-cyan" front

plot '$(morse_spectra_sam_data)' u 1:2 w l lw 3 dt 1 lc "blue" t "TDVP", \\
     ''                          u 1:4 w l lw 3 dt 6 lc "red"  t "Heller"

# Semiclassical spectrum is plotted as inset on top of the object 1
set tics font ",20"

set size   0.50, 0.45
set origin 0.44, 0.50

unset key
set key top right Left horizontal reverse enhanced autotitles nobox \\
spacing 1.00 width -1.0 height 0.5 noopaque

unset object 1

unset arrow

set xtics 5000
set mxtics 5
unset xlabel

set yrange [ 0.0 : 0.60 ]
set ytics 0.2
set mytics 2
unset ylabel

plot '$(morse_spectra_sam_data)' u 1:3 w l lw 2 dt 1 lc "black" t "Ext. semicl"

unset multiplot
"""
    run(pipeline(IOBuffer(gnuplot_script), `/usr/bin/env gnuplot`, fig4_file))

    # Online spectrum plot similar to Fig. 4
    @debug begin
        io = IOBuffer()
        println(io, "Drawing $(figure()): Spectrum")
        let fq=Constants.Eh2invcm*(2pi*Spectra.freqs(sol_hq.t, dt) .+ Vᴹ.D)
            p = plot(fq, abs.(spectrum_hq),
                     xlabel = "Frequency, cm⁻¹", ylabel="Intensity",
                     xlim = (0.0, 10_000.0),
                     label = "TDVP",
                     reuse = false)
            plot!(fq, abs.(spectrum_sc),
                  label = "Ext. semiclassical")
            plot!(fq,abs.(spectrum_hlr),
                  label = "Heller")
            gui()
            String(take!(io))
        end
    end

    println("Expansion of the initial GWP in the egenfunction basis")
    G0(x)    = GaussianToolbox.gwp(x, collect(inicn)...)
    G0_norm2 = GaussianToolbox.gaussian_norm2(collect(inicn)...) # Analytic

    # Extended trapezoidal rule. It appears that something like Romberg
    # integration would be better, but it is fooled by long zero tails of the
    # integrand on large grids. On the other hand, trapezoidal rule converges
    # exponentially for such functions.
    trapezoidal_rule(v, dx=1) = (sum(v) - (v[1]+v[end])/2)*dx

    # Eigenfunctions of the Morse potential are easier numerically
    nlvl = 24
    Eᴹ, Ψᴹ, grid = energy_levels(x -> Vᴹ(x), nlvl;
                                 start=0.0, stop=17.0, length=257,
                                 m=m, extended=true)
    dx = grid[2]-grid[1]
    @debug "Grid spacing" dx
    @debug begin
        Eᴹ_ana = energy_levels(Vᴹ, nlvl, m=m) .- Vᴹ.D
        Δ = Eᴹ - Eᴹ_ana
        """ Eigenvalues: 
        *   Numerical     Analytic      Δ"""
    end hcat(Eᴹ, Eᴹ_ana, Δ)

    G0_grid = G0.(grid)./sqrt(G0_norm2) # Normalized GWP on a grid
    @debug("Normalization check (numeric) for the normalized initial GWP: ",
           trapezoidal_rule(abs2.(G0_grid), dx))

    # Re-normalize eigenfunctions
    Ψᴹ ./= sqrt(dx)
    @debug("Norm squared of the eigenvectors: ",
           Ψᴹ_norm2 = [trapezoidal_rule(abs2.(Ψᴹ[:,k]), dx) for k=1:nlvl])

    expansion_coeffs =
        [trapezoidal_rule(conj(Ψᴹ[:,k]).*G0_grid, dx) for k=1:nlvl]
    expansion_weights = abs2.(expansion_coeffs)
    @debug "Expansion weights: " expansion_weights
    @debug "Sum of weights (completeness test)" sum(expansion_weights)

    scaled_expansion_weights = expansion_weights/expansion_weights[1]
    @info "Scaled expansion weights (to the first one) " scaled_expansion_weights
    println()

    # @debug begin
    #     plot(grid, real(G0_grid), reuse=false)
    #     gui()
    #     plot(grid, Ψᴹ[:, 1:nlvl], reuse=false)
    #     gui()
    # end

    Dynamics.resetglobals!()
end # let
println()
