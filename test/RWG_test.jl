using WaveguideModes
using WaveguideModes: WaveguideModes
using Unitful: ustrip
using Test

@testset "RWG" begin
    rwg1 = RWG("WR62")
    @test rwg1.a ≈ 0.015798799999999998
    @test rwg1.b ≈ 0.007899399999999999
    @test rwg1.ϵᵣ == 1.0
    @test rwg1.tanδ == 0.0 
    @test rwg1.σ == Inf
    @test rwg1.Rq == 0.0

    rwg2 = RWG("WR10")
    WaveguideModes.setup_modes!(rwg2, 75e9, 1)
    rwg3 = RWG("WR10", σ=1e100)
    WaveguideModes.setup_modes!(rwg3, 75e9, 1)
    @test rwg2.modes[1].Z ≈ rwg3.modes[1].Z
    @test rwg2.modes[1].γ ≈ rwg3.modes[1].γ

    rwg = RWG("WR10", σ=58e6, Rq=2000e-9)
    WaveguideModes.setup_modes!(rwg, 59.4e9, 1)
    @test rwg.modes[1].p == TE
    @test rwg.modes[1].m == 1
    @test rwg.modes[1].n == 0
    @test rwg.modes[1].γ ≈ 11.588267065714206 + 221.24049205550418im
    @test rwg.modes[1].Z ≈ 5.66501215556774 + 0.2871826803079236im

    rwg2 = RWG(a=0.00254/2, b=0.00254, σ=58e6, Rq=2000e-9)
    WaveguideModes.setup_modes!(rwg2, 59.4e9, 1)
    @test rwg2.modes[1].p == TE
    @test rwg2.modes[1].m == 0
    @test rwg2.modes[1].n == 1
    @test rwg2.modes[1].γ ≈ 11.588267065714206 + 221.24049205550418im
    @test rwg2.modes[1].Z ≈ 5.66501215556774 + 0.2871826803079236im

    a, b = lookup_rwg("WR10")
    f = 59u"GHz"
    γ1, Z1 = rwgte10gz(a, b, f)
    γ2, Z2 = rwgte10gz("WR10", f)
    γ3, Z3 = rwgte10gz(ustrip(u"m", a), ustrip(u"m", b), ustrip(u"Hz", f))
    @test γ1 ≈ γ2 ≈ γ3
    @test Z1 ≈ Z2 ≈ Z3
end
