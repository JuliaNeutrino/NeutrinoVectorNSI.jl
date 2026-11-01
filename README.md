# NeutrinoVectorNSI.jl

Vector (neutral-current) non-standard neutrino interactions with matter for [GLoBES.jl](https://github.com/JuliaNeutrino/GLoBES.jl). The matter potential becomes V [diag(1, 0, 0) + ε], with ε Hermitian and in units of V = √2 G_F N_e.

```julia
using Pkg
Pkg.develop(path = "/path/to/NeutrinoVectorNSI.jl")   # not registered yet
using GLoBES, NeutrinoVectorNSI
```

`VectorNSI` is a GLoBES.jl propagation block, with parameters `ε_αβ` and the phases `δε_αβ`:

```julia
model = NewPhysics(VectorNSI(:ee, :eτ))   # VectorNSI() for all six couplings
use_model!(model)                         # the probability engine of all experiments
p = OscParams(model, nufit_params("6.1"); ε_eτ = 0.2, δε_eτ = π / 2)
probability(p, :μ, :e; E = 2.5, L = 1300, ρ = 2.848)
```

It combines with the other blocks, e.g. [NeutrinoScalarNSI.jl](https://github.com/JuliaNeutrino/NeutrinoScalarNSI.jl) and [NeutrinoDecoherence.jl](https://github.com/JuliaNeutrino/NeutrinoDecoherence.jl), and rates, χ² and fits in GLoBES.jl use its parameters like the standard ones.

Tests: `julia --project=. test/runtests.jl`.
