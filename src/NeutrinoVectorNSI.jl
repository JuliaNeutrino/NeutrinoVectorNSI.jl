"""
Vector (neutral-current) non-standard neutrino interactions with matter for
GLoBES.jl, as a propagation block that combines with the other blocks in
`NewPhysics`.
"""
module NeutrinoVectorNSI

using GLoBES
using GLoBES:
    PropagationBlock,
    flavour_pairs,
    pair_label,
    coupling_parameters,
    coupling_aliases,
    add_couplings!
import GLoBES: block_parameters, parameter_aliases, matter_correction!

export VectorNSI

"""
    VectorNSI(couplings...)

Neutral-current (vector) non-standard interactions with matter: the matter
term of the Hamiltonian becomes V [diag(1, 0, 0) + ε], with ε Hermitian in
units of the standard potential V = √2 G_F N_e (so ε_αβ is the effective
coupling to electrons, ε = Σ_f (N_f/N_e) ε^f).

Couplings are flavour pairs (`:ee`, `:eμ`, `:eτ`, `:μμ`, `:μτ`, `:ττ`, or
ASCII `:emu`, ...; all six by default). Parameters: `ε_αβ` and, off the
diagonal, the phase `δε_αβ` (ASCII `eps_etau`, `deltaeps_etau`, ...).

```julia
VectorNSI(:ee, :eτ)   # couplings: parameters ε_ee, ε_eτ and δε_eτ
VectorNSI()   # all six couplings
```
"""
struct VectorNSI <: PropagationBlock
    pairs::Vector{Tuple{Int, Int}}
    VectorNSI(couplings...) = new(flavour_pairs(couplings))
end

block_parameters(b::VectorNSI) = coupling_parameters(b.pairs, :ε, :δε)
parameter_aliases(b::VectorNSI) = coupling_aliases(b.pairs, :ε, :δε, :eps, :deltaeps)
matter_correction!(ε, b::VectorNSI, p, x) = add_couplings!(ε, b.pairs, x, :ε, :δε)

Base.show(io::IO, b::VectorNSI) =
    print(io, "VectorNSI(", join((pair_label(p...) for p in b.pairs), ", "), ")")

end
