module WaveguideModes

using LinearAlgebra: det, norm, svd
using Unitful: Unitful, @u_str
using MetalSurfaceImpedance: Zsurface, effective_conductivity
using Printf: @printf
using Accessors: @set, @reset
using PrettyTables: pretty_table, ft_printf, tf_unicode_rounded, tf_html_default,
                    HtmlTableFormat
using AbstractPlutoDingetjes: is_inside_pluto
using FunctionZeros: besselj_zero, besselj_deriv_zero
using StaticArrays: @SVector, @SMatrix
using SimpleNonlinearSolve: SimpleNonlinearSolve as snls
using PRIMA: bobyqa

export @u_str
export TE, TM, setup_modes!,
       RWG, RWGMode, lookup_rwg, rwgte10gz, rwg_modes, rwg_modetable,
       CWG, CWGMode

using PhysicalConstants: PhysicalConstants
const c₀ = Unitful.ustrip(Float64, u"m/s", PhysicalConstants.CODATA2022.c_0)
const ϵ₀ = Unitful.ustrip(Float64, u"F/m", PhysicalConstants.CODATA2022.ε_0)
const μ₀ = Unitful.ustrip(Float64, u"H/m", PhysicalConstants.CODATA2022.μ_0)
const η₀ = Unitful.ustrip(Float64, u"Ω", PhysicalConstants.CODATA2022.Z_0)

abstract type Waveguide end
abstract type MetallicWaveguide <: Waveguide end
abstract type HomogeneousMetallicWaveguide <: MetallicWaveguide end

"""
    abstract type Mode end

Arbitrary electromagnetic mode.  Should contain the following fields:
- `γ::ComplexF64`: Attenuation constant in z direction [np/m]
"""
abstract type Mode end

"""
    abstract type MetalPipeMode <: Mode end

Electromagnetic mode supported by metal pipes, possibly filled with inhomogeneous dielectric.
"""
abstract type MetalPipeMode <: Mode end

"""
    abstract type HomogeneousMetalPipeMode <: MetalPipeMode end

Electromagnetic mode supported by homogeneously filled metal pipes.  Concrete subtypes should contain
at least the following fields:
- `p::TETM`: Mode classification (TE or TM).
- `f::Float64`: Frequency [Hz].
- `γ::ComplexF64`: Frequency-dependent attenuation constant [np/m].
- `Z::ComplexF64`: Frequency-dependent wave impedance [Ω].
"""
abstract type HomogeneousMetalPipeMode <: MetalPipeMode end

@enum TETM::Bool TE=true TM=false

isTE(mode::HomogeneousMetalPipeMode) = mode.p == TE
isTE(p::TETM) = p == TE
isTM(mode::HomogeneousMetalPipeMode) = mode.p == TM
isTM(p::TETM) = p == TM

"""
    update(mode::HomogeneousMetalPipeMode; f::Real, γ::Complex, Z::Complex) -> mode

Create a copy of an existing mode with updated frequency-dependent fields and return the updated mode.
"""
function update(mode::HomogeneousMetalPipeMode; f::Real, γ::Complex, Z::Complex)
    @reset mode.f = f
    @reset mode.γ = γ
    @reset mode.Z = Z
    return mode
end

"""
    mysqrt(x)

Same as `sqrt` unless `sqrt(x)` is pure negative imaginary in which
case it returns `-sqrt(x)` (i.e., positive pure imaginary).
"""
mysqrt(x) = sqrt(x)
function mysqrt(z::Complex)
    ans = sqrt(z)
    return iszero(real(ans)) && imag(ans) < 0 ? -ans : ans
end

"""
    is_html_environment()

Returns `true` if code is running under `IJulia` or `Pluto`, and `false` otherwise.
"""
function is_html_environment()
    return (isdefined(Main, :IJulia) && Main.IJulia.inited) || is_inside_pluto()
end

include("RWGModes.jl")

include("CWGModes.jl")

end # module
