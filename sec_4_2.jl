section("Static limit for the Morse potential, Sec. 4.2")
let Vᴹ = MorseModel.MorsePotential(), m = MorseModel.m, ħ = MorseModel.ħ

    #= The TDVP-based model in the static limit minimizes the full U, the
    extended semiclassical model minimizes U1, while the Heller's model
    minimizes bare V for q and U1 for w. =#
    def_full = Exact{MorseModel.MorsePotential}(m)
    def_u1   = Taylor{MorseModel.MorsePotential, 2}(m)
    guess = (q=Vᴹ.a, w=0.1)

    Umin_full, _pos = findmin(def_full, guess=guess)
    qmin_full, wmin_full = _pos.q, _pos.w
    @info("Stationary WP for the TDVP-based model",
          qmin=qmin_full, wmin=wmin_full, Umin=Umin_full)
    println()

    Umin_extsc, _pos = findmin(def_u1, guess=guess)
    qmin_extsc, wmin_extsc = _pos.q, _pos.w
    @info("Stationary WP for the extended semicalssical model",
          qmin=qmin_extsc, wmin=wmin_extsc, Umin=Umin_extsc)
    println()

    @debug begin
        qmin_extsc_alt =
            -log1p(-3wmin_extsc^2*Vᴹ.b^2/(2 + 4wmin_extsc^2*Vᴹ.b^2))/Vᴹ.b + Vᴹ.a
        @assert qmin_extsc ≈ qmin_extsc_alt
        "Verify analytic expression for qmin for U1, Eq. 36: "
    end qmin_num=qmin_extsc qmin_ana=qmin_extsc_alt
    println()

    qmin_hl = Vᴹ.a
    Umin_hl, wmin_hl =
        with_logger(NullLogger()) do # run silently
            WavepacketTools.minimal_uncertainty(def_u1, qmin_hl;
                                                guess=guess.w)
        end
    @info("Stationary WP for Heller's model",
          qmin = qmin_hl, wmin = wmin_hl, Umin = Umin_hl)
    println()

    @debug begin
        ωₑ = MorseModel.harmonic_frequency(Vᴹ, m)
        "Numerical/analytic wmin: "
    end wmin_num = wmin_hl wmin_ana = √(ħ/(2m*ωₑ))

    println('-'^115)
    println("Global minimum of the full model lies                   ",
            (Umin_full + Vᴹ.D)*Constants.Eh2invcm, " cm⁻¹ ",
            "above the bottom of the well")

    println("Global minimum of the extended semiclassical model lies ",
            (Umin_extsc + Vᴹ.D)*Constants.Eh2invcm, " cm⁻¹ ")

    println("Global minimum of Heller's model lies                   ",
            (Umin_hl + Vᴹ.D)*Constants.Eh2invcm, " cm⁻¹ ",
            "// definition 1: from U1",)
    println("                                                        ",
            (U(def_full, qmin_hl, wmin_hl)+ Vᴹ.D)*Constants.Eh2invcm, " cm⁻¹ ",
            "// definition 2: from the full U")

    levels = energy_levels(Vᴹ, 1, m=m)
    println("The ground vibrational state is at:                     ",
            levels[1]*Constants.Eh2invcm, " cm⁻¹ ")
    println('-'^115)
    println()

    xrange = LinRange(1.0, 4.0, 201)

    # Minimum uncertainty potential (can be computed only on a grid)
    U_min_uncert =
        with_logger(NullLogger()) do # run silently
            map(xrange) do x
                Umin, _ = WavepacketTools.minimal_uncertainty(def_full, x;
                                                              guess=guess.w)
                Umin
            end
        end

    # Write data files
    fig1_data = joinpath(datasave_path, "morse_pots.dat")
    open(fig1_data, "w") do io
        println(io, "# x, bohr    V, Eh        Umin_un, Eh        GS lvl, Eh")
        for (i, x) in enumerate(xrange)
            @printf(io, " %8.5f %12.7f %12.7f %12.7f\n",
                    x,
                    Vᴹ(x),
                    U_min_uncert[i],
                    1.4 < x < 2.2 ? levels[1] - Vᴹ.D : NaN)
        end
    end

    # Plotting via gnuplot
    fig1_file      = joinpath(figures_path, "fig1.eps")
    # Forward reference to the end of Sec. 5D : line and label for
    # ext.semiclassical trapping. Must be defined in Julia code to benefit from
    # string interpolation
    inipt = 1.4
    gnuplot_script = """
#!/usr/bin/gnuplot
set terminal postscript eps enhanced color \\
    font "Helvetica" 24 lw 2.0 dl 1.0 \\
    size 1.2*5, 1.2*3.5

set key top right Left reverse enhanced autotitles nobox \\
spacing 1.00 width -2 height 1.5

f(x) = x < $(inipt) || x > 2.5 ? NaN : $(Vᴹ(inipt))

unset label
set label 1 at $(qmin_hl),   $(-Vᴹ.D)     "" point pt 7 ps 2.5
set label 2 at $(qmin_full), $(Umin_full) "" point pt 9 ps 3
set label 3 at $inipt, $(Vᴹ(inipt)) center "{/ZapfDingbats \\073}" font ",32"

set xrange [ $(xrange.start) : $(xrange.stop) ]
set xtics 0.5
set mxtics 5
set xlabel "Coordinate, bohr"

set yrange [ -0.225 : 0.0 ]
set ytics 0.05
set mytics 5
set ylabel "Energy, E_h"


plot '$(fig1_data)' u 1:2 w l lw 3 lc "blue"            t "V_M", \\
     ''             u 1:3 w l lw 3 lc "red"    dt 3     t "min. uncertainty", \\
     ''             u 1:4 w l lw 2 lc "black"  dt "-" not
     # f(x)                   lw 4 lc "gray"   lt 0   not
"""
    run(pipeline(IOBuffer(gnuplot_script), `/usr/bin/env gnuplot`, fig1_file))
end # let
println()
