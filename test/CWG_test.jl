using WaveguideModes
using WaveguideModes: WaveguideModes
using Unitful: ustrip
using Test

@testset "CWG Basic Structure" begin
    # Test CWG creation with default parameters
    cwg1 = CWG(a=0.01)
    @test cwg1.a == 0.01
    @test cwg1.l == 0.0
    @test cwg1.ϵᵣ == 1.0
    @test cwg1.tanδ == 0.0
    @test cwg1.σ == Inf
    @test cwg1.Rq == 0.0
    @test isempty(cwg1.modes)

    # Test CWG with custom parameters
    cwg2 = CWG(a=0.005, l=0.1, ϵᵣ=2.0, tanδ=0.01, σ=58e6, Rq=1000e-9)
    @test cwg2.a == 0.005
    @test cwg2.l == 0.1
    @test cwg2.ϵᵣ == 2.0
    @test cwg2.tanδ == 0.01
    @test cwg2.σ == 58e6
    @test cwg2.Rq == 1000e-9
end

@testset "CWGMode Basic Structure" begin
    # Test TE mode creation
    mode_te = CWGMode(TE, 1, 1, 5000.0, 10e9, 0.5+0.1im, 50.0)
    @test mode_te.p == TE
    @test mode_te.m == 1
    @test mode_te.n == 1
    @test mode_te.f == 10e9
    @test mode_te.γ == 0.5 + 0.1im
    @test mode_te.Z == 50.0

    # Test TM mode creation
    mode_tm = CWGMode(TM, 0, 1, 3000.0, 10e9, 0.3+0.05im, 75.0)
    @test mode_tm.p == TM
    @test mode_tm.m == 0
    @test mode_tm.n == 1

    # Test invalid mode indices
    @test_throws ArgumentError CWGMode(TE, 0, 0, 5000.0, 10e9, 0, 0)  # m=0, n=0 not allowed
    @test_throws ArgumentError CWGMode(TE, -1, 1, 5000.0, 10e9, 0, 0)  # m<0 not allowed
    @test_throws ArgumentError CWGMode(TM, 0, 0, 5000.0, 10e9, 0, 0)   # n=0 not allowed
end

@testset "CWG setup_modes! - Basic" begin
    cwg = CWG(a=0.01)
    setup_modes!(cwg, 10e9, 4)
    @test length(cwg.modes) == 4
    @test all(m.f == 10e9 for m in cwg.modes)
    
    # Modes should be sorted by cutoff frequency
    kcos = [m.kco for m in cwg.modes]
    @test issorted(kcos)
end

@testset "CWG setup_modes! - Mode Selection" begin
    cwg = CWG(a=0.01)
    setup_modes!(cwg, 10e9, 4, ms=[0, 1])
    
    # Check that modes are created with specified m values
    m_values = unique([m.m for m in cwg.modes])
    @test all(m in [0, 1] for m in m_values)
    
    # Verify mode indices are valid
    for mode in cwg.modes
        if mode.p == TE
            @test mode.m ≥ 0 && mode.n ≥ 1
        else  # TM
            @test mode.m ≥ 0 && mode.n ≥ 1
        end
    end
end

@testset "CWG setup_modes! - Update Existing Modes" begin
    cwg = CWG(a=0.01)
    f1 = 10e9
    f2 = 15e9
    
    setup_modes!(cwg, f1, 3)
    modes1 = copy(cwg.modes)
    
    # Update modes to different frequency
    setup_modes!(cwg, f2, 3)
    modes2 = cwg.modes
    
    # Same number of modes
    @test length(modes1) == length(modes2)
    
    # Same mode indices
    for (m1, m2) in zip(modes1, modes2)
        @test m1.m == m2.m
        @test m1.n == m2.n
        @test m1.p == m2.p
    end
    
    # Different frequency and possibly different γ and Z
    @test all(m.f == f2 for m in modes2)
    @test any(m1.γ ≠ m2.γ for (m1, m2) in zip(modes1, modes2))
end

@testset "CWG setup_modes! - Different Methods" begin
    cwg_yeap = CWG(a=0.01, σ=58e6)
    cwg_abe = CWG(a=0.01, σ=58e6)
    cwg_ploss = CWG(a=0.01, σ=58e6)
    
    f = 10e9
    setup_modes!(cwg_yeap, f, 2, method=:yeap)
    setup_modes!(cwg_abe, f, 2, method=:abe)
    setup_modes!(cwg_ploss, f, 2, method=:ploss)
    
    @test length(cwg_yeap.modes) == length(cwg_abe.modes) == length(cwg_ploss.modes)
    
    # Results should differ between methods (but not be zero)
    for (m_yeap, m_abe, m_ploss) in zip(cwg_yeap.modes, cwg_abe.modes, cwg_ploss.modes)
        @test m_yeap.f == m_abe.f == m_ploss.f
        # Results can differ between methods
        @test !iszero(m_yeap.γ) || !iszero(m_abe.γ)
    end
end

@testset "CWG with Lossy Dielectric" begin
    cwg = CWG(a=0.01, ϵᵣ=2.0, tanδ=0.01)
    setup_modes!(cwg, 10e9, 2)
    
    @test all(m.f == 10e9 for m in cwg.modes)
    # Real parts of γ should be non-zero due to dielectric loss
    @test any(real(m.γ) > 0 for m in cwg.modes)
end

@testset "CWG with Surface Roughness" begin
    cwg_smooth = CWG(a=0.01, σ=58e6, Rq=0.0)
    cwg_rough = CWG(a=0.01, σ=58e6, Rq=100e-9)
    
    setup_modes!(cwg_smooth, 15e9, 2)
    setup_modes!(cwg_rough, 15e9, 2)
    
    # Rougher surface should have higher attenuation for lossy conductors
    for (m_smooth, m_rough) in zip(cwg_smooth.modes, cwg_rough.modes)
        @test m_smooth.f == m_rough.f
        # Real part of γ (attenuation) should be higher with roughness
        @test real(m_rough.γ) ≥ real(m_smooth.γ)
    end
end

@testset "CWG cwgkcoa function" begin
    # Test some known Bessel function zeros
    # TE₁₁: first zero of J'₁ ≈ 1.841
    # TM₀₁: first zero of J₀ ≈ 2.405
    # TE₂₁: first zero of J'₂ ≈ 3.054
    
    kco_te11 = WaveguideModes.cwgkcoa(TE, 1, 1)
    kco_tm01 = WaveguideModes.cwgkcoa(TM, 0, 1)
    kco_te21 = WaveguideModes.cwgkcoa(TE, 2, 1)
    
    @test kco_te11 ≈ 1.841 rtol=0.01
    @test kco_tm01 ≈ 2.405 rtol=0.01
    @test kco_te21 ≈ 3.054 rtol=0.01
    
    # Higher order n should have higher kco for same p and m
    kco_tm01 = WaveguideModes.cwgkcoa(TM, 0, 1)
    kco_tm02 = WaveguideModes.cwgkcoa(TM, 0, 2)
    @test kco_tm02 > kco_tm01
end

@testset "CWG Lossless Case (Perfect Conductor)" begin
    cwg = CWG(a=0.01, σ=Inf, ϵᵣ=1.0, tanδ=0.0)
    setup_modes!(cwg, 20e9, 3)
    
    # For perfect conductor and lossless dielectric, attenuation should be zero
    for mode in cwg.modes
        if real(mode.γ) ≠ 0  # Allow for numerical errors
            @test abs(real(mode.γ)) < 1e-10
        end
    end
end

@testset "CWG Mode Impedance Properties" begin
    cwg = CWG(a=0.01)
    setup_modes!(cwg, 10e9, 3)
    
    for mode in cwg.modes
        # Wave impedance should be non-zero complex
        @test !iszero(mode.Z)
        # For passive waveguide, Z should have positive real part
        if real(mode.Z) < 0
            # This shouldn't happen in normal cases
            @test false
        end
    end
end

@testset "CWG Error Handling" begin
    cwg = CWG(a=0.01)
    setup_modes!(cwg, 10e9, 3)
    
    # Cannot call setup_modes with different nmodes when modes exist
    @test_throws ArgumentError setup_modes!(cwg, 15e9, 4)
    
    # nmodes must be positive
    cwg_empty = CWG(a=0.01)
    @test_throws ArgumentError setup_modes!(cwg_empty, 10e9, 0)
    @test_throws ArgumentError setup_modes!(cwg_empty, 10e9, -1)
end
