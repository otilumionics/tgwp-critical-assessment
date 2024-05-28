using PrettyTables

section("Fourier grid Hamiltonian for H₂")
println("""
Reproduce Table II from Ref. [^1]. We report 12 decimals instead of 8.

Note that even the exact values do not match to the data in Table II of
Ref.[^1]. Moreover, values obtained here with the FGH method for N=129 are
*significantly* more accurate compared with the analytical ones than reported in
Ref.[^1]. Possible sources of discrepancy could be:

 1. Difference in parsing floating-point values by old and modern computers

 2. Limited precision in diagonalization or in matrix elements, e.g. the use of
 single-precision (32 bit) floats

 3. Different value of the proton mass.

[^1]: C. C. Marston and G. G. Balint-Kurti. “The Fourier grid Hamiltonian method
      for bound state eigenvalues and eigenfunctions,” J. Chem. Phys. 91,
      3571-3576 (1989).
""")
let
    V_H₂ = MorseModel.MorsePotential(D=0.1744, a=1.40201, b=1.02764)
    μ = 1836.15267344/2         # H₂ reduced mass, modern value
    n = 17                      # Number of vibrational  levels
    Ns = [65, 129]              # Fourier grid size(s)
    tp = 7.8126663425715        # Classical turning point for E₁₆ level, bohr
    k = length(Ns)
    levels_num = zeros(Float64, n, k)

    # Generic wrapper to call numerical minimization
    (::GenericPotential)(x::Number; der::Int=0) = V_H₂(x, der=der)

    for i=1:k
        levels_num[:, i] =
            energy_levels(GenericPotential(), n, m=μ,
                          start=0.0, stop=1.5*tp, length=Ns[i])
    end
    levels_num .+= V_H₂.D       # To make all E[i] > 0
    levels_exact = energy_levels(V_H₂, n, m=μ)

    header = (vcat("v",          # Vibrational quantum number
                  ["FGH method for N=$N" for N in Ns],
                   "Exact eigenvalues"),
              vcat("", fill("Eₕ", k+1)))
    pretty_table(hcat(collect(Integer, 0:(n-1)), levels_num, levels_exact);
                 backend = Val(:text),
                 alignment = :c,
                 columns_width = 0,
                 autowrap = true,
                 linebreaks = true,
                 header=header,
                 formatters=ft_printf("%.12f", 2:(k+2)))
    println()
end
println()
