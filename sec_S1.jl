section("Matching quartic potential to the ammonia model")
let Vamm = AmmoniaModel.AmmoniaPotential()
    # ... and let x₀ be the position of the rightmost (positive) minimum of Vamm
    x₀ = 1.0                    # guess
    Vamm_min, x₀ = findmin(Vamm, guess=x₀)

    # Vq = ax⁴ + bx³ + cx² + dx + e
    #
    # (1)  Vq(0) = e = Vamm(0)
    #
    # (2)  Vq' = 4ax³ + 3bx² + 2cx + d = 0  at x=(0, ±x₀)
    #      Vq' = 4ax(x² - x₀²) = 0
    #      b = 0
    #      c = -2ax₀²
    #      d = 0
    # (3)  Vq(±x₀) = Vamm_min
    #      -axₒ⁴ + e = -axₒ⁴ + Vamm(0) = Vamm_min
    #       a = (Vamm(0) - Vamm_min)/x₀⁴
    a = (Vamm(0) - Vamm_min)/x₀^4
    b = 0
    c = -2a*x₀^2
    d = 0
    e = Vamm(0)
    Vq = QuarticModel.QuarticPotential(reverse([a, b, c, d, e]))
    @info "Matched parameters for Vq = ax⁴ + bx³ + cx² + dx + e:" a b c d e

    x = LinRange(-1.25, 1.25, 101)
    @debug begin
        io = IOBuffer()
        title = "Model ammonia potential and its quartic approximation"
        println(io, "Drawing A$(figure()): ", title)
        plot(x, [Vamm.(x) Vq.(x)],
             title=title,
             labels=["Vamm" "Vq"],
             xlabel = "x, bohr",
             ylabel = "V(x), Eₕ",
             reuse=false)
        gui()
        String(take!(io))
    end
end
