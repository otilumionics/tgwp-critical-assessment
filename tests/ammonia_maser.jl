using PrettyTables

section("Ammonia (NH₃) maser")
println("""
Reproduce the transition frequency table from
https://chem.libretexts.org/Bookshelves/Physical_and_Theoretical_Chemistry_Textbook_Maps/Quantum_Tutorials_(Rioux)/04%3A_Spectroscopy/4.04%3A_The_Ammonia_Inversion_and_the_Maser

NOTE 1. It seems there is a typo in `c` value. Recalculated from the original
paper [^1] `c = 1.3696`, not `1.36696`.

NOTE 2. The method for finding eigenvalues is not clear from the Web page above.
It cites a paper [^2] which is a short description of a software product.
However, we kept grid specification as x ∈ [-2.0, 2.0], the number of points N = 99.

[^1]: J. D. Swalen and J. A. Ibers, “Potential function for the inversion of
      ammonia,” J. Chem. Phys. 36, 1914–1918 (1962).

[^2]: J. C. Hansen J. Chem. Educ. Software, 8C Number 2 (1996).
""")
let
    V = AmmoniaModel.AmmoniaPotential()
    m = AmmoniaModel.m          # Effective mass for the inversion
    n = 4                       # Number of states
    levels_num = energy_levels(V, n, m=m,
                               # Grid parameters
                               start=-2.0, stop=2.0, length=99)
    transitions =
        ["E₀ → E₁" (levels_num[2] - levels_num[1])*Constants.Eh2invcm   0.79
         "E₂ → E₃" (levels_num[4] - levels_num[3])*Constants.Eh2invcm  36.0
         "E₁ → E₂" (levels_num[3] - levels_num[2])*Constants.Eh2invcm 932.5
         "E₀ → E₃" (levels_num[4] - levels_num[1])*Constants.Eh2invcm 968.3]
    header = (["Transition", "Theory", "Experiment"],
              ["",           "cm⁻¹",   "cm⁻¹"])
    pretty_table(transitions;
                 backend=Val(:text),
                 header=header,
                 formatters=ft_printf("%0.2f", 2:2))
end
println()
