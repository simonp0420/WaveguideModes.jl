module WaveguideModes

using LinearAlgebra: det, norm, svd
using Unitful: Unitful, @u_str
using MetalSurfaceImpedance: Zsurface, effective_conductivity
using Printf: @printf
using Accessors: @set, @reset
using PrettyTables: fmt__printf, HtmlTableFormat, MarkdownStyle, MarkdownTableFormat, 
                    HtmlTableStyle, pretty_table, TextHighlighter, TextTableStyle
using AbstractPlutoDingetjes: is_inside_pluto
using StaticArrays: @SVector, @SMatrix
using SimpleNonlinearSolve: SimpleNonlinearSolve as snls
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

Returns `true` if code is running under `IJulia`, `Pluto`, or `BonitoBook` and `false` otherwise.
"""
function is_html_environment()
    return (isdefined(Main, :IJulia) && Main.IJulia.inited) || is_inside_pluto() || 
    any(t -> occursin("IJulia", t) || occursin("Bonito", t) || occursin("Pluto", t), 
        string.(typeof.(Base.Multimedia.displays)))
end

"""
    display_mode_table(modedata::Array, titlestr::String, freq_unit, length_unit, col_fmts)

Display a mode table to the user's console, or to an HTML environment.

## Positional Arguments
- `modedata`: An array of size `(nrows,6)` containing the modal data to be output.  The columns
  contain ["Type (TE or TM)", "m", "n", "Cutoff Freq.", "Guide Wavelength", "Attenuation"]
- `thtlestr`: A `String` containing the table title.
- `freq_unit`: Either a `String` or an object (such as a Unitful unit) that can be converted
  to `String` to use as a column label for the frequency column.
- `length_unit`: Either a `String` or an object (such as a Unitful unit) that can be converted
  to `String` to use as a column label for the Guide Wavelength column.
- `col_fmts`: A `Vector{String}` of length 6 containing the `@sprinf` columnar formats for the
  numeric entries in the table.
"""
function display_mode_table(modedata::Array, titlestr::String, freq_unit, length_unit, col_fmts)
    size(modedata, 2) == 6 || throw(ArgumentError("Wrong number of columns in modedata"))
    maxtitlelen = 60
    title2 = String[]
    nexttline = ""
    for phrase in eachsplit(titlestr, ',')
        if length(nexttline) + length(phrase) > maxtitlelen
            push!(title2, nexttline)
            nexttline = ""
        end
        nexttline *= string(phrase, ",")
    end
    push!(title2, nexttline)
    title = first(title2)
    subtitle = length(title2) > 1 ? join(title2[begin+1:end])[begin+1:end-1] : ""
    column_labels = [["Type", "m", "n", "Cutoff Freq.", "Guide Wavelength", "Attenuation"]]
    column_units = ["", "", "", "[$freq_unit]", "[$length_unit]", "[dB/$length_unit]"]
    column_label_alignment = :c
    push!(column_labels, column_units)
    formatters = [fmt__printf(col_fmts[i], [i]) for i in 1:6]

    if is_html_environment()

        p = pretty_table(
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

    else

        println()
        p = pretty_table(
            modedata;
            style = TextTableStyle(; title = crayon"bold", subtitle = crayon"bold"),
            backend = :text,
            title,
            subtitle,
            column_labels,
            formatters,
            #highlighters = [TextHighlighter((d, i, j) -> iseven(i), crayon"negative")],
        )
    end
    return nothing
end

include("RWGModes.jl")

include("CWGModes.jl")

end # module
