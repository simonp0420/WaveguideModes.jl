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
    rwg2 = RWG("WR10", sigma=58u"MS/m", Rq=2000u"nm")
    rwg3 = RWG(a=0.00254, b=0.00127, sigma=58u"MS/m", Rq=2000u"nm")
    WaveguideModes.setup_modes!(rwg, 59.4e9, 1)
    WaveguideModes.setup_modes!(rwg2, 59.4e9, 1)
    WaveguideModes.setup_modes!(rwg3, 59.4e9, 1)
    @test rwg.modes[1].p == rwg2.modes[1].p == rwg3.modes[1].p == TE
    @test rwg.modes[1].m == rwg2.modes[1].m == rwg3.modes[1].m == 1
    @test rwg.modes[1].n == rwg2.modes[1].n == rwg3.modes[1].n == 0
    @test rwg.modes[1].γ == rwg2.modes[1].γ == rwg3.modes[1].γ ≈ 11.588267065714206 + 221.24049205550418im
    @test rwg.modes[1].Z == rwg2.modes[1].Z == rwg3.modes[1].Z ≈ 5.66501215556774 + 0.2871826803079236im

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

@testset "alphaploss" begin
    rwg = RWG("WR34", σ = 0.3 * 58u"MS/m")
    setup_modes!(rwg, 25e9, 1, method=:ploss)
    @test real(rwg.modes[1].γ) ≈ 0.09533851108191686
end

@testset "RWG Mode Structure and Validation" begin
    # Test RWGMode creation with valid indices
    mode_te10 = RWGMode(TE, 1, 0, 1000.0, 10e9, 0.1+0.5im, 50.0)
    @test mode_te10.p == TE
    @test mode_te10.m == 1
    @test mode_te10.n == 0
    
    mode_tm11 = RWGMode(TM, 1, 1, 1500.0, 10e9, 0.2+0.3im, 75.0)
    @test mode_tm11.p == TM
    @test mode_tm11.m == 1
    @test mode_tm11.n == 1
    
    # Test invalid TE modes
    @test_throws ArgumentError RWGMode(TE, 0, 0, 1000.0, 10e9, 0, 0)  # (0,0) not allowed for TE
    
    # Test invalid TM modes
    @test_throws ArgumentError RWGMode(TM, 0, 1, 1500.0, 10e9, 0, 0)   # m must be > 0 for TM
    @test_throws ArgumentError RWGMode(TM, 1, 0, 1500.0, 10e9, 0, 0)   # n must be > 0 for TM
end

@testset "RWG with Different Dimensions" begin
    # Test WR-series waveguides with various sizes
    wg_names = ["WR10", "WR34", "WR62"]
    
    for wg_name in wg_names
        rwg = RWG(wg_name)
        @test rwg.a > 0
        @test rwg.b > 0
        @test rwg.a ≠ rwg.b  # All WR guides have a ≠ b
        @test rwg.a > rwg.b  # Convention: a > b
    end
end

@testset "RWG setup_modes! - Mode Sorting" begin
    rwg = RWG("WR10")
    setup_modes!(rwg, 50e9, 8)
    
    # Modes should be sorted by cutoff frequency
    kcos = [m.kco for m in rwg.modes]
    @test issorted(kcos)
    
    # First mode should be TE10
    @test rwg.modes[1].p == TE
    @test rwg.modes[1].m == 1
    @test rwg.modes[1].n == 0
end

@testset "RWG Frequency Dependence" begin
    rwg = RWG("WR34")
    f1, f2 = 10e9, 30e9
    
    setup_modes!(rwg, f1, 3)
    modes_f1 = copy(rwg.modes)
    
    setup_modes!(rwg, f2, 3)
    modes_f2 = rwg.modes
    
    # Mode indices should be same
    for (m1, m2) in zip(modes_f1, modes_f2)
        @test (m1.m, m1.n, m1.p) == (m2.m, m2.n, m2.p)
    end
    
    # Propagation constants should differ
    @test any(m1.γ ≠ m2.γ for (m1, m2) in zip(modes_f1, modes_f2))
end

@testset "RWG Lossy Conductor - Conductivity Effect" begin
    a, b = 0.01, 0.005
    f = 32e9
    nmodes = 2
    
    rwg_lossless = RWG(a=a, b=b, σ=Inf)
    rwg_lossy = RWG(a=a, b=b, σ=58e6)
    
    setup_modes!(rwg_lossless, f, nmodes)
    setup_modes!(rwg_lossy, f, nmodes)
    
    # For lossy conductor, real part of γ should be larger
    for (m_ll, m_l) in zip(rwg_lossless.modes, rwg_lossy.modes)
        @test real(m_l.γ) ≥ real(m_ll.γ)
    end
end

@testset "RWG with Dielectric Loss" begin
    rwg_lossless = RWG("WR34", ϵᵣ=2.0, tanδ=0.0)
    rwg_lossy = RWG("WR34", ϵᵣ=2.0, tanδ=0.05)
    
    f = 20e9
    setup_modes!(rwg_lossless, f, 2)
    setup_modes!(rwg_lossy, f, 2)
    
    # Lossy dielectric should increase attenuation
    for (m_ll, m_l) in zip(rwg_lossless.modes, rwg_lossy.modes)
        @test real(m_l.γ) ≥ real(m_ll.γ)
    end
end

@testset "RWG setup_modes! - Method Comparison" begin
    rwg_auto = RWG("WR10", σ=58e6, Rq=500e-9)
    rwg_ploss = RWG("WR10", σ=58e6, Rq=500e-9)
    
    f = 50e9
    setup_modes!(rwg_auto, f, 2, method=:auto)
    setup_modes!(rwg_ploss, f, 2, method=:ploss)
    
    # Both should produce valid results
    @test all(!iszero(m.γ) for m in rwg_auto.modes)
    @test all(!iszero(m.γ) for m in rwg_ploss.modes)
    
    # Results may differ between methods
    # (just verify both produced results)
    @test length(rwg_auto.modes) == length(rwg_ploss.modes)
end

@testset "RWG Perfect Conductor (Lossless)" begin
    rwg = RWG("WR34")
    setup_modes!(rwg, 35e9, 3)
    
    # For perfect conductor and lossless dielectric, 
    # real part of γ should be zero
    for mode in rwg.modes
        @test iszero(real(mode.γ))
    end
end

@testset "RWG rwgkco Calculation" begin
    a, b = 0.01, 0.005
    
    # TE10 mode should have lower kco than TE20 or TE01
    kco_te10 = WaveguideModes.rwgkco(a, b, 1, 0)
    kco_te20 = WaveguideModes.rwgkco(a, b, 2, 0)
    kco_te01 = WaveguideModes.rwgkco(a, b, 0, 1)
    
    @test kco_te10 < kco_te20
    @test kco_te10 < kco_te01
    
    # kco should be symmetric in m and n when a and b are swapped
    kco_te10_swapped = WaveguideModes.rwgkco(b, a, 0, 1)
    @test kco_te10 ≈ kco_te10_swapped
end

@testset "RWG Update Function" begin
    mode = RWGMode(TE, 1, 0, 1000.0, 10e9, 0.1+0.5im, 50.0)
    
    f_new = 15e9
    γ_new = 0.15 + 0.6im
    Z_new = 55.0 + 0.2im
    
    mode_updated = WaveguideModes.update(mode; f=f_new, γ=γ_new, Z=Z_new)
    
    @test mode_updated.f == f_new
    @test mode_updated.γ == γ_new
    @test mode_updated.Z == Z_new
    
    # Mode indices should be unchanged
    @test mode_updated.m == mode.m
    @test mode_updated.n == mode.n
    @test mode_updated.p == mode.p
    @test mode_updated.kco == mode.kco
end

@testset "RWG Unitful Quantities" begin
    # Test with Unitful quantities for conductivity and roughness
    rwg = RWG(a=10u"mm", b=5u"mm", σ=58u"MS/m", Rq=500u"nm")
    @test rwg.a ≈ 0.01
    @test rwg.b ≈ 0.005
    @test rwg.σ ≈ 58e6
    @test rwg.Rq ≈ 500e-9
end

@testset "RWG High Frequency Behavior" begin
    # At very high frequencies, multiple modes should propagate
    rwg = RWG("WR10")
    f_low = 10e9
    f_high = 100e9
    
    setup_modes!(rwg, f_low, 5)
    modes_low = copy(rwg.modes)
    
    setup_modes!(rwg, f_high, 5)
    modes_high = rwg.modes
    
    # Number of propagating modes typically increases with frequency
    @test length(modes_high) ≥ length(modes_low)
    
    # At high frequency, TE10 mode should be well above cutoff
    te10_high = findfirst(m -> m.p == TE && m.m == 1 && m.n == 0, modes_high)
    @test !isnothing(te10_high)
    # Imaginary part (phase constant) should be large
    @test abs(imag(modes_high[te10_high].γ)) > abs(imag(modes_low[1].γ))
end

@testset "RWG Mode Update Consistency" begin
    rwg = RWG("WR34")
    
    # Create and update multiple times
    setup_modes!(rwg, 10e9, 3)
    modes_10 = [deepcopy(m) for m in rwg.modes]
    
    setup_modes!(rwg, 20e9, 3)
    modes_20 = [deepcopy(m) for m in rwg.modes]
    
    setup_modes!(rwg, 10e9, 3)
    modes_10_again = rwg.modes
    
    # Going back to 10 GHz should reproduce same results
    for (m_orig, m_again) in zip(modes_10, modes_10_again)
        @test m_orig.γ ≈ m_again.γ
        @test m_orig.Z ≈ m_again.Z
    end
end

@testset "RWG lookup_rwg function" begin
    # Test standard waveguide lookups
    a_wr10, b_wr10 = WaveguideModes.lookup_rwg("WR10")
    a_wr34, b_wr34 = WaveguideModes.lookup_rwg("WR34")
    
    @test a_wr10 < a_wr34  # WR10 is smaller than WR34
    
    # All dimensions should be positive
    @test all(>(0u"inch"), (a_wr10, b_wr10, a_wr34,b_wr34))
end

@testset "RWG Error Handling" begin
    rwg = RWG("WR34")
    setup_modes!(rwg, 20e9, 3)
    
    # Can't update with different number of modes
    @test_throws ArgumentError setup_modes!(rwg, 20e9, 4)
    
    # Empty waveguide requires nmodes > 0
    rwg_empty = RWG("WR34")
    @test_throws ArgumentError setup_modes!(rwg_empty, 20e9, 0)
    @test_throws ArgumentError setup_modes!(rwg_empty, 20e9, -1)
end

