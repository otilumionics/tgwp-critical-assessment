section("Accuracy of numerical approximations for the ammonia model, Sec. 4.1")
let Vamm = AmmoniaModel.AmmoniaPotential, m = AmmoniaModel.m

    println("Taylor series convergence radius, 1/√2c = ", 1/√(2Vamm().c))

    # Data for Table II
    ngh, ntl = 10, 10           # Number of Gauss-Hermite and Taylor terms
    Approximations = [(GaussHermite, ngh), (Taylor, ntl)]
    xrange = LinRange(-2.0, 2.0, 101)
    wlist  = [0.1, 0.6, 1.0]
    @info "Variables domain:" xrange wlist
    println()

    println("Table 2 numerical data in LaTeX format")
    println('-'^60)
    for w in wlist
        println("\\multicolumn{3}{c}{\$w=\\SI{$w}{\\bohr}\$} \\\\")
        for (Approximation, n) in Approximations
            deviation = map(xrange) do q
                U(Exact{Vamm}(m), q, w) - U(Approximation{Vamm, n}(m), q, w)
            end
            # @show deviation
            # @info("Stats for $(Approximation)",
            #       meanabs = round(   mean(abs, deviation), sigdigits=2),
            #       absmax  = round(maximum(abs, deviation), sigdigits=2))
            println("$Approximation & \\num{",
                    round(   mean(abs, deviation), sigdigits=2),
                    "} & \\num{",
                    round(maximum(abs, deviation), sigdigits=2),
                    "} \\\\")
        end
    end
    println('-'^60)
end # let
println()
