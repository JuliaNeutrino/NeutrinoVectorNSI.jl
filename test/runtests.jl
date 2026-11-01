using NeutrinoVectorNSI, GLoBES, Test

const DUNE = (; L = 1300, ρ = 2.848)
const E = [0.5, 1.0, 2.5, 5.0]

@testset "parameters" begin
    m = NewPhysics(VectorNSI(:ee, :eτ))
    @test m.names == [:ε_ee, :ε_eτ, :δε_eτ]
    @test length(NewPhysics(VectorNSI()).names) == 6 + 3
    q = OscParams(m, nufit_params("6.1"); eps_ee = 0.2, deltaeps_etau = 1.0)
    @test (q.ε_ee, q.ε_eτ, q.δε_eτ) == (0.2, 0.0, 1.0)
    @test_throws ArgumentError VectorNSI(:ex)
    @test sprint(show, VectorNSI(:ee, :μτ)) == "VectorNSI(ee, μτ)"
end

@testset "probabilities" begin
    p = nufit_params("6.1")
    standard = [probability(p, :μ, :e; E, DUNE..., antineutrino = a) for a in (false, true)]
    m = NewPhysics(VectorNSI(:ee, :eτ, :μτ))
    use_model!(m)

    # ε = 0: the standard probabilities
    @test probability(OscParams(m, p), :μ, :e; E, DUNE...) ≈ standard[1] rtol = 1e-12
    @test probability(OscParams(m, p), :μ, :e; E, DUNE..., antineutrino = true) ≈ standard[2] rtol =
        1e-12

    # ε_ee = x is the standard potential times 1 + x, for ν and ν̄
    nsi = [
        probability(OscParams(m, p; ε_ee = 0.5), :μ, :e; E, DUNE..., antineutrino = a) for
        a in (false, true)
    ]
    standard_engine!()
    for (i, a) in enumerate((false, true))
        @test nsi[i] ≈ probability(p, :μ, :e; E, L = 1300, ρ = 1.5 * 2.848, antineutrino = a) rtol =
            1e-10
    end

    # No effect in vacuum
    use_model!(m)
    q = OscParams(m, p; ε_eτ = 0.3, δε_eτ = 1.2, ε_μτ = 0.05)
    standard_engine!()
    vacuum = probability(p, :μ, :e; E, L = 1300, ρ = 0.0)
    use_model!(m)
    @test probability(q, :μ, :e; E, L = 1300, ρ = 0.0) ≈ vacuum rtol = 1e-10

    # Off-diagonal NSI change P(μ → e) in matter, and probabilities stay unitary
    @test !(probability(q, :μ, :e; E, DUNE...) ≈ standard[1])
    total = sum(probability(q, :μ, β; E, DUNE..., antineutrino = true) for β in (:e, :μ, :τ))
    @test total ≈ ones(length(E)) rtol = 1e-12

    # A phase π flips the sign of the coupling
    @test probability(OscParams(q; ε_eτ = 0.3, δε_eτ = π), :μ, :e; E, DUNE...) ≈
          probability(OscParams(q; ε_eτ = -0.3, δε_eτ = 0.0), :μ, :e; E, DUNE...) rtol = 1e-10
    standard_engine!()
end
