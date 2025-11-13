using WaveguideModes
using WaveguideModes: WaveguideModes
using Unitful: ustrip
using Test

@testset "RWG and CWG Comparison - Basic Properties" begin
    # Create equivalent guides with similar fundamental modes
    rwg = RWG("WR62")
    cwg = CWG(a=0.01)
    
    f = 10e9
    setup_modes!(rwg, f, 2)
    setup_modes!(cwg, f, 2)
    
    # Both should have modes
    @test length(rwg.modes) ≥ 1
    @test length(cwg.modes) ≥ 1
    
    # TE10 (RWG) and TE11 (CWG) are fundamentals of rectangular and circular guides
    te10_rwg = findfirst(m -> m.p == TE && m.m == 1 && m.n == 0, rwg.modes)
    te11_cwg = findfirst(m -> m.p == TE && m.m == 1 && m.n == 1, cwg.modes)
    
    @test !isnothing(te10_rwg)
    @test !isnothing(te11_cwg)
end

@testset "Mode Sorting Consistency" begin
    rwg = RWG("WR34")
    cwg = CWG(a=0.01)
    
    f = 15e9
    setup_modes!(rwg, f, 6)
    setup_modes!(cwg, f, 6)
    
    # RWG modes should be sorted by cutoff frequency
    rwg_kcos = [m.kco for m in rwg.modes]
    @test issorted(rwg_kcos)
    
    # CWG modes should be sorted by cutoff frequency
    cwg_kcos = [m.kco for m in cwg.modes]
    @test issorted(cwg_kcos)
end

@testset "Frequency Sweep - Both Guide Types" begin
    rwg = RWG("WR90")  
    cwg = CWG(a=0.01)
    
    frequencies = [10e9, 15e9, 20e9]
    
    for f in frequencies
        setup_modes!(rwg, f, 3)
        setup_modes!(cwg, f, 3)
        
        @test all(m.f == f for m in rwg.modes)
        @test all(m.f == f for m in cwg.modes)
        
        # At each frequency, fundamental mode should be above cutoff
        @test abs(imag(rwg.modes[1].γ)) > 0
        @test abs(imag(cwg.modes[1].γ)) > 0
    end
end

@testset "Loss Scaling - Both Guide Types" begin
    # Lossy vs lossless comparison
    rwg_lossless = RWG("WR34", σ=Inf)
    rwg_lossy = RWG("WR34", σ=58e6)
    cwg_lossless = CWG(a=0.01, σ=Inf)
    cwg_lossy = CWG(a=0.01, σ=58e6)
    
    f = 35e9
    setup_modes!(rwg_lossless, f, 2)
    setup_modes!(rwg_lossy, f, 2)
    setup_modes!(cwg_lossless, f, 2)
    setup_modes!(cwg_lossy, f, 2)
    
    # Loss should increase attenuation in both guide types
    for (m_ll, m_l) in zip(rwg_lossless.modes, rwg_lossy.modes)
        @test real(m_l.γ) ≥ real(m_ll.γ)
    end
    
    for (m_ll, m_l) in zip(cwg_lossless.modes, cwg_lossy.modes)
        @test real(m_l.γ) ≥ real(m_ll.γ)
    end
end

@testset "Mode Impedance - Comparison" begin
    rwg = RWG("WR34")
    cwg = CWG(a=0.01)
    
    f = 15e9
    setup_modes!(rwg, f, 3)
    setup_modes!(cwg, f, 3)
    
    # Wave impedances should be complex and non-zero
    for mode in rwg.modes
        @test !iszero(mode.Z)
        @test isa(mode.Z, Complex)
    end
    
    for mode in cwg.modes
        @test !iszero(mode.Z)
        @test isa(mode.Z, Complex)
    end
end

@testset "Different Conductivity Values" begin
    # Test with various conductivity levels
    σ_values = [Inf, 1e8, 58e6, 1e6]
    
    rwg_list = [RWG("WR34", σ=σ) for σ in σ_values]
    cwg_list = [CWG(a=0.01, σ=σ) for σ in σ_values]
    
    f = 20e9
    for (rwg, cwg) in zip(rwg_list, cwg_list)
        setup_modes!(rwg, f, 2)
        setup_modes!(cwg, f, 2)
        
        @test length(rwg.modes) ≥ 1
        @test length(cwg.modes) ≥ 1
    end
end

@testset "Dielectric Loss - Both Guides" begin
    rwg_no_loss = RWG("WR34", ϵᵣ=1.0, tanδ=0.0)
    rwg_with_loss = RWG("WR34", ϵᵣ=2.0, tanδ=0.05)
    cwg_no_loss = CWG(a=0.01, ϵᵣ=1.0, tanδ=0.0)
    cwg_with_loss = CWG(a=0.01, ϵᵣ=2.0, tanδ=0.05)
    
    f = 35e9
    setup_modes!(rwg_no_loss, f, 2)
    setup_modes!(rwg_with_loss, f, 2)
    setup_modes!(cwg_no_loss, f, 2)
    setup_modes!(cwg_with_loss, f, 2)
    
    # Dielectric loss should increase attenuation for RWG
    for (m1, m2) in zip(rwg_no_loss.modes, rwg_with_loss.modes)
        @test real(m2.γ) ≥ real(m1.γ)
    end
    
    # Dielectric loss should increase attenuation for CWG
    for (m1, m2) in zip(cwg_no_loss.modes, cwg_with_loss.modes)
        @test real(m2.γ) ≥ real(m1.γ)
    end
end

@testset "Type Predicates - TE/TM Detection" begin
    rwg = RWG("WR34")
    cwg = CWG(a=0.01)
    
    setup_modes!(rwg, 15e9, 3)
    setup_modes!(cwg, 15e9, 3)
    
    # Check TE/TM predicates for RWG
    for mode in rwg.modes
        if mode.p == TE
            @test WaveguideModes.isTE(mode)
            @test !WaveguideModes.isTM(mode)
        else
            @test WaveguideModes.isTM(mode)
            @test !WaveguideModes.isTE(mode)
        end
    end
    
    # Check TE/TM predicates for CWG
    for mode in cwg.modes
        if mode.p == TE
            @test WaveguideModes.isTE(mode)
            @test !WaveguideModes.isTM(mode)
        else
            @test WaveguideModes.isTM(mode)
            @test !WaveguideModes.isTE(mode)
        end
    end
end

@testset "Surface Roughness Effect - Both Guides" begin
    # RWG with and without roughness
    rwg_smooth = RWG("WR34", σ=58e6, Rq=0.0)
    rwg_rough = RWG("WR34", σ=58e6, Rq=100e-9)
    
    # CWG with and without roughness
    cwg_smooth = CWG(a=0.01, σ=58e6, Rq=0.0)
    cwg_rough = CWG(a=0.01, σ=58e6, Rq=100e-9)
    
    f = 35e9
    setup_modes!(rwg_smooth, f, 2)
    setup_modes!(rwg_rough, f, 2)
    setup_modes!(cwg_smooth, f, 2)
    setup_modes!(cwg_rough, f, 2)
    
    # Roughness should increase attenuation for RWG
    for (m_smooth, m_rough) in zip(rwg_smooth.modes, rwg_rough.modes)
        @test real(m_rough.γ) ≥ real(m_smooth.γ)
    end
    
    # Roughness should increase attenuation for CWG
    for (m_smooth, m_rough) in zip(cwg_smooth.modes, cwg_rough.modes)
        @test real(m_rough.γ) ≥ real(m_smooth.γ)
    end
end

@testset "RWG with All Constructor Variants" begin
    # Standard constructor
    rwg1 = RWG(a=0.01, b=0.005)
    @test rwg1.a == 0.01
    @test rwg1.b == 0.005
    
    # String spec constructor
    rwg2 = RWG("WR34")
    @test rwg2.a > 0
    @test rwg2.b > 0
    
    # Both should work for mode setup
    setup_modes!(rwg1, 15e9, 2)
    setup_modes!(rwg2, 15e9, 2)
    
    @test length(rwg1.modes) ≥ 1
    @test length(rwg2.modes) ≥ 1
end

@testset "Propagation Constant Consistency" begin
    # For perfect conductor and lossless dielectric,
    # imaginary part of γ should equal β
    rwg = RWG("WR34", σ=Inf, ϵᵣ=1.0, tanδ=0.0)
    cwg = CWG(a=0.01, σ=Inf, ϵᵣ=1.0, tanδ=0.0)
    
    f = 35e9
    k₀ = 2π * f / WaveguideModes.c₀
    
    setup_modes!(rwg, f, 2)
    setup_modes!(cwg, f, 2)
    
    for mode in rwg.modes
        # For lossless case, β = -imag(γ)
        β = -imag(mode.γ)
        @test !iszero(β)
    end
    
    for mode in cwg.modes
        β = -imag(mode.γ)
        @test !iszero(β)
    end
end

@testset "CWG Methods Consistency" begin
    # Test that different methods produce valid (if different) results
    cwg_yeap = CWG(a=0.01, σ=58e6)
    cwg_abe = CWG(a=0.01, σ=58e6)
    cwg_ploss = CWG(a=0.01, σ=58e6)
    
    f = 10e9
    setup_modes!(cwg_yeap, f, 2, method=:yeap)
    setup_modes!(cwg_abe, f, 2, method=:abe)
    setup_modes!(cwg_ploss, f, 2, method=:ploss)
    
    # All methods should produce results
    for (m_y, m_a, m_p) in zip(cwg_yeap.modes, cwg_abe.modes, cwg_ploss.modes)
        # All should have same frequency
        @test m_y.f == m_a.f == m_p.f
        
        # Mode indices should match
        @test (m_y.m, m_y.n, m_y.p) == (m_a.m, m_a.n, m_a.p) == (m_p.m, m_p.n, m_p.p)
        
        # γ values should be complex and non-zero (or very small for lossless)
        @test isa(m_y.γ, Complex)
        @test isa(m_a.γ, Complex)
        @test isa(m_p.γ, Complex)
    end
end

@testset "RWG kwarg Equivalence" begin
    # Test that different keyword forms produce same results
    f = 15e9
    
    # Using σ vs sigma
    rwg1 = RWG(a=0.01, b=0.005, σ=58e6)
    rwg2 = RWG(a=0.01, b=0.005, sigma=58e6)
    
    setup_modes!(rwg1, f, 2)
    setup_modes!(rwg2, f, 2)
    
    for (m1, m2) in zip(rwg1.modes, rwg2.modes)
        @test m1.γ == m2.γ
        @test m1.Z == m2.Z
    end
    
    # Using ϵᵣ vs epsr
    rwg3 = RWG(a=0.01, b=0.005, ϵᵣ=2.0)
    rwg4 = RWG(a=0.01, b=0.005, epsr=2.0)
    
    setup_modes!(rwg3, f, 2)
    setup_modes!(rwg4, f, 2)
    
    for (m3, m4) in zip(rwg3.modes, rwg4.modes)
        @test m3.γ == m4.γ
        @test m3.Z == m4.Z
    end
end

@testset "Mode Cutoff Behavior" begin
    rwg = RWG("WR34")
    cwg = CWG(a=0.01)
    
    # At very low frequency (near/below cutoff)
    f_low = 1e9
    setup_modes!(rwg, f_low, 3)
    setup_modes!(cwg, f_low, 3)
    
    # Some modes may be evanescent (cutoff)
    # Real part of γ should be significant for evanescent modes
    for mode in rwg.modes
        if real(mode.γ) > 0.1
            # This is an evanescent mode
            @test abs(imag(mode.γ)) < abs(real(mode.γ))
        end
    end
    
    for mode in cwg.modes
        if real(mode.γ) > 0.1
            # This is an evanescent mode
            @test abs(imag(mode.γ)) < abs(real(mode.γ))
        end
    end
end

@testset "Mode Count Variation with Frequency" begin
    rwg = RWG("WR34")
    
    # At low frequency, fewer modes propagate
    setup_modes!(rwg, 5e9, 10)
    nmodes_low = length(rwg.modes)
    
    # At high frequency, more modes can propagate
    setup_modes!(rwg, 30e9, 10)
    nmodes_high = length(rwg.modes)
    
    # More modes should propagate at higher frequency
    @test nmodes_high ≥ nmodes_low
end
