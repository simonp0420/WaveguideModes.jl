module WaveguideModes

using LinearAlgebra: det, norm, svd, eigen, Hermitian
using Unitful: Unitful, @u_str
using MetalSurfaceImpedance: Zsurface, effective_conductivity
using Printf: @printf
using Accessors: @set, @reset
using PrettyTables: fmt__printf, HtmlTableFormat, MarkdownStyle, MarkdownTableFormat, 
                    HtmlTableStyle, pretty_table, TextHighlighter, TextTableStyle
using AbstractPlutoDingetjes: is_inside_pluto
using StaticArrays: @SVector, @SMatrix, SVector
using SimpleNonlinearSolve: SimpleNonlinearSolve as snls
using NonlinearSolve: NonlinearSolve as nls
using PRIMA: bobyqa
using FunctionZeros: besselj_zero, besselj_deriv_zero
using Crayons: @crayon_str

export @u_str
export TE, TM, setup_modes!,
       RWG, RWGMode, lookup_rwg, rwgte10gz, rwg_modes, rwg_modetable,
       CWG, CWGMode, cwg_modes, cwg_modetable

using PhysicalConstants: PhysicalConstants
const c₀ = Unitful.ustrip(Float64, u"m/s", PhysicalConstants.CODATA2022.c_0)
const ϵ₀ = Unitful.ustrip(Float64, u"F/m", PhysicalConstants.CODATA2022.ε_0)
const μ₀ = Unitful.ustrip(Float64, u"H/m", PhysicalConstants.CODATA2022.μ_0)
const η₀ = Unitful.ustrip(Float64, u"Ω", PhysicalConstants.CODATA2022.Z_0)

const np2dB = 20 * log10(exp(1)) # Convert neper to dB


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
    ModeDataTable <: Any

A struct containing data to be displayed for a few modes of a waveguide

## Fields
- `title`: A `String` containing the table title.
- `modedata`: An array of size `(nrows,ncols)` containing the modal data to be output. 
- `coltitles`: A vector of `ncols` strings containing the titles for each column.
- `colunits`:  A vector of `ncols` strings containing the units used for each column.
- `col_fmts`: A `Vector{String}` of length `ncols` containing the `@sprinf` columnar formats for the
  numeric entries in the table.
"""
@kwdef struct ModeDataTable
    title::String
    modedata::Array{Any}
    coltitles::Vector{String}
    colunits::Vector{String}
    colfmts::Vector{String}

    function ModeDataTable(title, modedata, coltitles, colunits, colfmts)
        ncols = size(modedata, 2)
        length(coltitles) == ncols || error("length(coltitles) ≠ size(modedata, 2)")
        length(colunits) == ncols || error("length(colunits) ≠ size(modedata, 2)")
        length(colfmts) == ncols || error("length(colfmts) ≠ size(modedata, 2)")
        return new(title, modedata, coltitles, colunits, colfmts)
    end
end

function Base.show(io::IO, ::MIME"text/plain", mdt::ModeDataTable)
    (;  title, modedata, coltitles, colunits, colfmts) = mdt
    maxtitlelen = 60
    title2 = String[]
    nexttline = ""
    for phrase in eachsplit(title, ',')
        if length(nexttline) + length(phrase) > maxtitlelen
            push!(title2, nexttline)
            nexttline = ""
        end
        nexttline *= string(phrase, ",")
    end
    push!(title2, nexttline)
    title = first(title2)
    subtitle = length(title2) > 1 ? join(title2[begin+1:end])[begin+1:end-1] : ""
    column_labels = [coltitles, colunits]
    formatters = [fmt__printf(colfmts[i], [i]) for i in 1:6]
    println(io)
    pretty_table(
        io,
        modedata;
        style = TextTableStyle(; title = crayon"bold", subtitle = crayon"bold"),
        backend = :text,
        title,
        subtitle,
        column_labels,
        formatters,
    )
end

function Base.show(io::IO, ::MIME"text/html", mdt::ModeDataTable)
    (;  title, modedata, coltitles, colunits, colfmts) = mdt
    maxtitlelen = 60
    title2 = String[]
    nexttline = ""
    for phrase in eachsplit(title, ',')
        if length(nexttline) + length(phrase) > maxtitlelen
            push!(title2, nexttline)
            nexttline = ""
        end
        nexttline *= string(phrase, ",")
    end
    push!(title2, nexttline)
    title = first(title2)
    subtitle = length(title2) > 1 ? join(title2[begin+1:end])[begin+1:end-1] : ""
    column_labels = [coltitles, colunits]
    formatters = [fmt__printf(colfmts[i], [i]) for i in 1:6]
    pretty_table(
        io,
        modedata;
        backend = :html,
        style = HtmlTableStyle(; 
                    title = ["font-weight" => "bold", "font-size" => "large"],
                    subtitle = ["font-weight" => "bold", "font-size" => "large"],
                    first_line_column_label = ["font-weight" => "bold"],
                    column_label = ["font-weight" => "regular", "font-style" => "italic"],
                ),        
        title,
        subtitle,
        column_labels,
        formatters,
    )
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

include("RWGModes.jl")

include("CWGModes.jl")

end # module
