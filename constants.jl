module Constants

export Eh2eV, Eh2kJ_per_mol, Eh2kcal_per_mol, Eh2K, Eh2invcm,
    au2fs, bohr2angs, Da2me

## Conversion factors
# Hartree equivalents are from https://en.wikipedia.org/wiki/Hartree#Other_relationships
const Eh2eV           = 27.211_386_245_988 # 1Eh in eV
const Eh2kJ_per_mol   = 2625.499_639_479_9 # 1Eh in kJ/mol
const Eh2kcal_per_mol = 627.509_474_063_1  # 1Eh in kcal/mol
const Eh2K            = 315_775.024_804_07 # 1Eh in K
const Eh2invcm        = 219_474.631_363_20 # 1Eh in cm^{-1}

# From https://en.wikipedia.org/wiki/Hartree_atomic_units
const au2fs           = 2.418_884_326_585_7e-2 # atomic unit of time in fs
const bohr2angs       = 0.529_177_210_903      # atomic unit of length in Å

# https://en.wikipedia.org/wiki/Dalton_(unit), 1/12 C-12
const Da2me           = 1_822.888_486_209 # 1Da in mass of electron

end # module Constants
