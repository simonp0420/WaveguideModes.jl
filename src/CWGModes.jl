"""
    CWGMode <: HomogeneousMetalPipeMode

Struct representing an electromagnetic mode of a circular waveguide, homogeneously filled with dielectric.

## Fields:
- `p::TETM`: Mode type.
- `m::Int`: Index for modal variation along ϕ direction.
- `n::Int`: Index for modal variation along ρ direction.
- `kco::Float64`: Cutoff wavenumber [radian/meter].
- `f::Float64`: Frequency [Hz].
- `γ::ComplexF64`: Attenuation constant [neper/meter].
- `Z::ComplexF64`: Modal wave impedance normalized to η₀, the impedance of free space [unitless].
"""
struct CWGMode <: HomogeneousMetalPipeMode
    p::TETM   # Mode type
    m::Int    # ϕ index
    n::Int    # ρ index
    kco::Float64 # cutoff wavenumber [rad/meter]
    f::Float64 # Frequency [Hz]
    γ::ComplexF64 # Attenuation constant [np/meter]
    Z::ComplexF64 # Wave impedance normalized to η₀ [unitless]

    function CWGMode(
            p::TETM, m::Integer, n::Integer, kco::Real, f::Real, γ::Number, Z::Number)
        (m ≥ 0 && n ≥ 1) || throw(ArgumentError("Illegal (m,n) = $((m,n))"))
        return new(p, m, n, kco, f, γ, Z)
    end
end

"""
    CWGMode(; p, m, n, a, f=0, γ=0, Z=0)

Convenient keyword-argument constructor that defaults the frequency-dependent fields to zero and uses waveguide radius in meters.
"""
CWGMode(; p, m, n, a, f = 0, γ = 0, Z = 0) = CWGMode(p, m, n, cwgkcoa(p, m, n) / a, f, γ, Z)

"""
    cwgkcoa(p::TETM, m::Integer, n::Integer)

Compute the the product of cutoff wavenumber times radius for a circular waveguide.
## Required Positional Arguments
- `p`: Mode type (`TE` or `TM`).
- `m`: Mode number for azimuthal variation. Must be nonnegative.
- `n`: Mode number for radial variation. Must be positive.
## Return Value
- `kcoa`: The cutoff wavenumber of the mode times the waveguide radius.
"""
function cwgkcoa(p::TETM, m::Integer, n::Integer)
        m ≥ 0 || throw(ArgumentError("m must be nonnegative"))
        n > 0 || throw(ArgumentError("n must be positive"))
        return isTM(p) ? besselj_zero(m, n) : besselj_deriv_zero(m, n)
end


"""
    CWG <: HomogeneousMetallicWaveguide

A struct representing a section of metallic, circular waveguide, uniformly filled with dielectric.
It is assumed that the (possibly non-smooth) metal walls are either perfectly conducting or form a 
"good conductor", in the sense that the standard perturbational formulas for losses in waveguide walls
hold for all modes.


## Fields
- `a::Float64`: The (inner) radius of the waveguide [m].
- `l::Float64`: The length of the waveguide section (along z) [m].
- `ϵᵣ::Float64`: The dielectric constant for the material filling the guide [unitless].
- `tanδ::Float64`: The loss tangent for the material filling the guide [unitless].
- `σ::Float54`: The bulk conductivity of the waveguide metal walls [S/m].
- `Rq::Float64`: The metal wall RMS surface roughness [m].
- `modes::Vector{CWGMode}`: The modes treated in this waveguide.  If nonempty, the list 
  should be sorted in order of increasing cutoff frequency.

Both positional and keyword argument constructors for `CWG` are available.
"""
@kwdef struct CWG <: HomogeneousMetallicWaveguide
    a::Float64  # The (inner) radius of the waveguide [m]
    l::Float64 = 0.0  # Length [meter]
    ϵᵣ::Float64 = 1.0    # Relative permittivity of dielectric filling waveguide
    tanδ::Float64 = 0.0 # Loss tangent of dielectric filling waveguide
    σ::Float64 = Inf # Wall bulk conductivity [S/m]
    Rq::Float64 = 0.0 # Wall RMS surface roughness [m]
    modes::Vector{CWGMode} = CWGMode[]
end


"""
    setup_modes!(wg::CWG, f::Real, nmodes::Integer=length(c.modes); ms=1:20, method=:yeap)

Set up the modes for a uniform circular waveguide.

## Positional Arguments
- `wg`: A `CWG` instance.  If `wg.modes` is empty, then it will be appended to `length(n)`.  After this, or
  if it is already allocated, each mode in the vector will be replaced by one with updated values of `γ` and 
  `Z` corresponding to the frequency `f`.
- `f`: The frequency [Hz].
- `nmodes`: (optional) The number of modes to append to `wg.modes` if it is empty.  If `wg.modes` is nonempty, then
  `nmodes` must be equal to `length(wg.modes)`.

## Keyword Arguments
- `ms`: (optional) An `Int` or `AbstractVector{Int}` denoting the modal `m` values (the azimuthal mode numbers)
  to use when adding modes to the waveguide.  Ignored if `wg.modes` is nonempty. Values must be nonnegative.
- `method::Symbol`:  `:yeap` (default), `:abe`, or `:ploss` denoting whether to use method of Reference [1], [2],
  or [3] respectively.

## References
[1] K. H. Yeap et al: "Attenuation in Circular and Rectangular Waveguides", Electromagnetics, Vol. 37, No. 3,
2017, pp. 171-184.

[2] T. Abe and Y. Yamaguchi, "Propagation Constant Below Cutoff Frequency in a Circular Waveguide with
Conducting Medium", IEEE Trans. MTT, Vol. MTT-29, no. 7, July 1981, pp. 707-712.

[3] R. Collin, **Field Theory of Guided Waves, 2cd Edition**, Table 5.5.
"""
function setup_modes!(wg::CWG, f::Real, nmodes::Integer = length(wg.modes); 
                      ms::Union{Int,AbstractVector{Int}}=1:20,
                      method::Symbol=:yeap)
    (; a, σ, Rq, ϵᵣ, tanδ) = wg
    if isempty(wg.modes)
        nmodes ≤ 0 && throw(ArgumentError("nmodes must be > 0"))
        all(≥(0), ms) || throw(ArgumentError("ms values must be nonnegative"))
        modes = [CWGMode(; p, m, n, a) for p in (TE, TM) for m in ms for n in 1:nmodes]
        sort!(modes; by = x -> x.kco)
        deleteat!(modes, (nmodes + 1):length(modes))
        append!(wg.modes, modes)
    else
        nmodes == length(wg.modes) ||
            throw(ArgumentError("nmodes not equal to number of existing modes in wg"))
    end

    k₀ = 2π * f / c₀ # free-space wavenumber [rad/m]
    rootϵ = mysqrt(ϵᵣ * complex(1.0, -tanδ))
    ka = k₀ * a * rootϵ # wavenumber in dielectric times radius
    ηnorm = inv(rootϵ) # η/η₀
    η = η₀ * ηnorm # intrinsic impedance in dielectric
    Zs = Zsurface(f, σ, Rq, :normal) # surface impedance ([Ω/□])
    Zsoη = Zs / η
    # Update γ and Z
    for (q, mode) in pairs(wg.modes)
        (; m, p, kco) = mode
        umn = kco * a 
        if method == :yeap
            βa = cwgkza_yeap(p, m, umn, ka, Zsoη)
        elseif method == :abe
            βa = cwgkza_abe(p, m, umn, ka, Zsoη)
        elseif method == :ploss
            βa = cwgkza_ploss(p, m, umn, ka, Zsoη)
        else
            Throw(ArgumentError("Unknown method: $method"))
        end
        β = βa / a
        γ = im * β
        if p == TE
            Z = ka / βa * ηnorm
        else
            Z = βa / ka * ηnorm
        end
        wg.modes[q] = update(mode; f, γ, Z)
    end
    return wg
end


"""
    cwgkza_yeap(p::TETM, m, umn, ka, Zsn)

Compute product of z-directed complex wavenumber kz and waveguide radius a for a circular waveguide.

## Arguments
- `p::TETM`: Mode type (`TE` or `TM`).
- `m::Integer`: Mode index for variation in ϕ direction.
- `umn::Float64`: Equal to kco * a, where kco is the cutoff wavenumber for the mode of interest. (For
  circular guides this is a zero of the Bessel function J_m or a zero of its derivative).
- `ka::ComplexF64`: Product of wavenumber (for the dielectric filling the guide) times the waveguide radius.
- `Zsn::ComplexF64`: The ratio of the metal wall surface impedance to the intrinsic impedance of the dielectric
  material filling the waveguide.

## Return Value
- 'kza::Complex64': Product of the waveguide radius with the z-directed wavenumber [1/m] including both
  propagation (imaginary part) and attenuation (real part).

## Reference
K. H. Yeap et al: "Attenuation in Circular and Rectangular Waveguides", Electromagnetics, Vol. 37, No. 3,
2017, pp. 171-184
"""
function cwgkza_yeap(p::TETM, m::Integer, umn::Float64, ka::ComplexF64, Zsn::ComplexF64)
    kra = umn
    kra² = kra^2
    ka² = ka^2
    if iszero(Zsn)
        γa =  mysqrt(kra² - ka²)
    else
        kza = sqrt(ka² - kra²)
        f1 = kra / (2 * ka)
        zsum = inv(Zsn) + Zsn
        root = sqrt(4 + (2m * kza / kra^2)^2 - zsum^2)
        if p == TE
            δ = f1 * (-im * zsum + root) / ((m / umn)^2 - 1)
        else
            δ = inv(-f1 * (im * zsum + root))
        end
        γa = mysqrt((umn + δ)^2 - ka^2)
    end
    kza = -im * γa
    return kza
end


"""
    cwgkza_abe(p::TETM, m, umn, ka, Zsn)

Compute product of z-directed complex wavenumber kz and waveguide radius a for a circular waveguide.

## Arguments
- `p::TETM`: Mode type (`TE` or `TM`).
- `m::Integer`: Mode index for variation in ϕ direction.
- `umn::Float64`: Equal to kco * a, where kco is the cutoff wavenumber for the mode of interest. (For
  circular guides this is a zero of the Bessel function J_m or a zero of its derivative).
- `ka::ComplexF64`: Product of wavenumber (for the dielectric filling the guide) times the waveguide radius.
- `Zsn::ComplexF64`: The ratio of the metal wall surface impedance to the intrinsic impedance of the dielectric
  material filling the waveguide.

## Return Value
- 'kza::Complex64': Product of the waveguide radius with the z-directed wavenumber [1/m] including both
  propagation (imaginary part) and attenuation (real part).

## Reference
T. Abe and Y. Yamaguchi, "Propagation Constant Below Cutoff Frequency in a Circular Waveguide with
Conducting Medium", IEEE Trans. MTT, Vol. MTT-29, no. 7, July 1981, pp. 707-712.
"""
function cwgkza_abe(p::TETM, m::Integer, umn::Float64, ka::ComplexF64, Zsn::ComplexF64)
    umn² = umn^2
    ka² = ka^2
    if iszero(Zsn)
        γa = mysqrt(umn² - ka²)
    else
        h0a² = ka² - umn²
        h0a = sqrt(h0a²)
        kea = ka / Zsn
        v = sqrt(kea^2 - h0a²)
        m² = m^2
        if p == TE
            #Δu = (m² * (ka² - umn²) + umn² * umn²) / (-im * umn * v * (m² - umn²)) # Eq. (7)
            kea² = kea^2
            Δu = -im * ((m * h0a * v)^2 + (umn²)^2 * kea²) / (v * umn * (ka² + kea^2) * (m² - umn²))
        else
            Δu = - ka² * umn * v / (2v * ka² + im * umn² * (ka² + kea^2))
        end
        γa = mysqrt((umn + Δu)^2 - ka^2)
    end
    kza = -im * γza
    return kza
end


"""
    cwgkza_ploss(p::TETM, m, umn, ka, Zsn)

Compute product of z-directed complex wavenumber kz and waveguide radius a for a circular waveguide. 
Estimate the attenuation for modes above cutoff using the power loss method.

## Arguments
- `p::TETM`: Mode type (`TE` or `TM`).
- `m::Integer`: Mode index for variation in ϕ direction.
- `umn::Float64`: Equal to kco * a, where kco is the cutoff wavenumber for the mode of interest. (For
  circular guides this is a zero of the Bessel function J_m or a zero of its derivative).
- `ka::ComplexF64`: Product of wavenumber (for the dielectric filling the guide) times the waveguide radius.
- `Zsn::ComplexF64`: The ratio of the metal wall surface impedance to the intrinsic impedance of the dielectric
  material filling the waveguide.

## Return Value
- 'kza::Complex64': Product of the waveguide radius with the z-directed wavenumber [1/m] including both
  propagation (imaginary part) and attenuation (real part).

## Reference
R. Collin, **Field Theory of Guided Waves, 2cd Edition**, Table 5.5.
"""
function cwgkza_ploss(p::TETM, m::Integer, umn::Float64, ka::ComplexF64, Zsn::ComplexF64)
    Rsn = real(Zsn)
    umn² = umn^2
    ka² = ka^2
    γa = mysqrt(umn² - ka²)
    if iszero(Rsn) || real(ka) ≤ umn
        αca = 0.0
    else
        kcook² = (umn / real(ka))^2
        αca = Rsn / sqrt(1 - kcook²)
        if p == TE
            m² = m^2
            αca *= kcook² + m² / (umn² - m²)
        end
    end
    γa += αca
    kza = -im * γa
    return kza
end


"""
    cwg_modes(a, f; ms = 1:20, nmodes=10, ϵᵣ=1.0, tanδ=0.0, σ=Inf*u"S/m", Rq=0.0u"m") -> modedata

Compute cutoff frequency, guide wavelength, and attenuation constant for first few modes of a circular waveguide.

Note: This function is intended for programmatic use.  For interactive use, see `cwg_modetable`.  


## Required Positional Arguments
- `a`: The waveguide inner radius as any `Unitful` length quantity, e.g. `a=0.8128u"cm"`, `a=320u"mil"`, 
  or `a=0.320u"inch"`.
- `f`: The frequency as any `Unitful` frequency quantity, e.g. `f=3u"GHz"`,  `f=3e9u"Hz"` or `f=3000u"MHz"`.

## Optional Keyword Arguments
- `nmodes`: The number of modes to be treated.  Defaults to `10`.
- `ms`: An `Int` or `AbstractVector{Int}` denoting the modal `m` values (the azimuthal mode numbers)
  to consider. Values must be nonnegative.  Defaults to `1:20`.
- `ϵᵣ` or `epsr`: The dielectric constant for the material filling the waveguide.
- `tanδ` or `tandel`: The loss tangent for the dielectric material filling the waveguide.
- `σ` or `sigma`: Bulk conductivity of the metal walls in any `Unitful` quantity with dimensions equal to
  those of `u"S/m"`. Defaults to `Inf*u"S/m"` (perfect electric conductor).
- `Rq`: Surface RMS roughness of the metal walls in any `Unitful` length dimension, e.g. `6u"μm"`.
  Defaults to `0u"m"` (perfectly smooth).  The surface roughness is used with `σ` to compute a (lower) effective 
  conductivity via the [`MetalSurfaceImpedance`](https://github.com/simonp0420/MetalSurfaceImpedance.jl) package.

## Return Value
- `modedata::Matrix{Any}`: A matrix of size `(nmodes, 6)`, containing in each row `[p, m, n, fco, λg, α]`
  for a single mode. Here `p = "TE"` for a TE mode or `"TM"` for a TM mode, `m` and `n` are integer mode
  indices for the azimuthal and radial mode variation, resp., `fco` is the cutoff frequency expressed in 
  the same units as input argument `f` (but without attached units), `λg` is the guide wavelength expressed
  in the same units as `a` (but without attached units), and `α` is the mode attenuation constant in units 
  of dB/unitlength (without attached units), where unitlength is one unit of length in the same units 
  as `a` and `b`. The modes are listed in order of increasing cutoff frequency. 
"""
function cwg_modes(
        a::Unitful.Quantity{<:Real, Unitful.dimension(u"m")},
        f::Unitful.Quantity{<:Real, Unitful.dimension(u"Hz")};
        nmodes::Int = 10,
        ms::Union{Int,AbstractVector{Int}} = 1:20,
        ϵᵣ::Real = 1.0,
        epsr::Real = 1.0,
        tanδ::Real = 0.0,
        tandel::Real = 0.0,
        σ::Unitful.Quantity{<:Real, Unitful.dimension(u"S/m")} = Inf * u"S/m",
        sigma::Unitful.Quantity{<:Real, Unitful.dimension(u"S/m")} = Inf * u"S/m",
        Rq::Unitful.Quantity{<:Real, Unitful.dimension(u"m")} = 0.0u"m",
)

    ϵᵣ = max(ϵᵣ, epsr)
    tanδ = max(tanδ, tandel)
    σ = min(σ, sigma)

    # Verify inputs
    a > 0u"m" || throw(ArgumentError("a must be positive"))
    length_unit = Unitful.unit(a)
    f > 0u"Hz" || throw(ArgumentError("f must be positive"))
    freq_unit = Unitful.unit(f)
    nmodes > 0 || throw(ArgumentError("nmodes must be positive"))
    all(≥(0), ms) || throw(ArgumentError("ms must all be nonnegative"))
    ϵᵣ ≥ 1 || throw(ArgumentError("ϵᵣ must be greater than or equal to 1"))
    tanδ ≥ 0 || throw(ArgumentError("tanδ must be nonnegative"))
    σ > 0u"S/m" || throw(ArgumentError("σ must be posiive"))
    Rq ≥ 0u"m" || throw(ArgumentError("Rq must be nonnegative"))

    # Convert to Float64 in standard MKS units and remove units. 
    am = Unitful.ustrip(Float64, u"m", a)
    fHz = Unitful.ustrip(Float64, u"Hz", f)
    σSpm = Unitful.ustrip(Float64, u"S/m", σ)
    Rqm = Unitful.ustrip(Float64, u"m", Rq)

    # Create list of modes to consider:
    nmax = nmodes
    modes = [CWGMode(; p, m, n, a = am) for p in (TE, TM) for m in ms for n in 1:nmax]
    sort!(modes; by = x -> x.kco)
    deleteat!(modes, (nmodes + 1):length(modes))
    wg = CWG(; a = am, ϵᵣ, tanδ, σ = σSpm, Rq = Rqm, modes)
    setup_modes!(wg, fHz)

    modedata = Matrix{Any}(undef, nmodes, 6)

    cover2π = c₀ / (sqrt(ϵᵣ) * 2π)
    for (i, mode) in pairs(wg.modes)
        λgm = 2π / imag(mode.γ)
        λg = Unitful.ustrip(length_unit, λgm * u"m")
        dBpm = np2dB * real(mode.γ)  # [dB/m]
        dBplu = Unitful.ustrip(length_unit^-1, dBpm * u"m^-1") # [dB/length_unit]
        fcohz = mode.kco * cover2π
        fco = Unitful.ustrip(freq_unit, fcohz * u"Hz")
        modedata[i, :] .= (string(mode.p), mode.m, mode.n, fco, λg, dBplu)
    end

    return modedata
end


"""
    cwg_modetable(; a, f; kwargs...)

Pretty-print a table of circular waveguide mode properties: cutoff frequency, guide wavelength, and attenuation constant.

Guide wavelength and attenuation are accurately computed at frequencies below, at, or 
above cutoff, for both smooth and rough imperfectly conducting surfaces, by combining the Gradient method
of Reference [1] with the method of Yeap, et al. [2].
 For other modes the standard 
perturbational formulas are used in conjunction with the Gradient method to determine effective surface impedance.  


Note: This function is intended for interactive use.  For programmatic use, see `cwg_modes`.  
The table is printed to the user's console using the `PrettyTables` package.  

## Required Keyword Arguments
- `a`: The waveguide inner radius as any `Unitful` length quantity, e.g. `a=0.8128u"cm"`, `a=320u"mil"`, 
  or `a=0.320u"inch"`.
- `f`: The frequency expressed as a `Unitful` quantity with dimensions of inverse time.  Examples: `2.35u"GHz"`,
  `2350u"MHz"`, `2.35e9u"Hz"`, etc.

## Optional Keyword Arguments
- `ms`: An `Int` or `AbstractVector{Int}` denoting the modal `m` values (the azimuthal mode numbers)
  to use when adding modes to the waveguide. Values must be nonnegative. Defaults to `1:99`.
- `length_unit`:: A `Unitful:Units` object with dimension of length. Defaults to the units used for `a`. 
  The table produced as output will print values in the "guide wavelength" column using this unit; similarly
  the attenuation column will use units of dB/`length_unit`.
- `ϵᵣ` or `epsr`: Dielectric constant of material filling the waveguide. Default value is `1.0`.
- `tanδ` or `tandel`: Loss tangent of material filling the waveguide. Default value is `0.0`. 
- `σ` or `sigma`: The bulk conductivity of the waveguide metal walls expressed as a `Unitful` quantity
  with dimensions the same as those of `u"S/m"`.  Default value is `Inf*u"S/m"`.
- `Rq`: The RMS surface roughness of the waveguide walls expressed as a `Unitful` quantity having dimension
  of length. Default value is `0.0u"m"`.
- `col_fmts::String = ["%s", "%i", "%i", "%#.8g", "%8.4f", "%8.4f"]`: A vector of C-style format strings
  used to format the six columns of the table.

## References
- [1] D. N. Grujić, "Simple and Accurate Approximation of Rough Conductor Surface Impedance," 
  **IEEE Trans. Microwave Theory Tech.**, vol. 70, no. 4, pp. 2053-2059, April 2022.
- [2] K. H. Yeap et al: "Attenuation in Circular and Rectangular Waveguides", **Electromagnetics**, Vol. 37, No. 3,
  2017, pp. 171-184.
"""
function cwg_modetable(;
        a::Unitful.Quantity{<:Real, Unitful.dimension(u"m")},
        f::Unitful.Quantity{<:Real, Unitful.dimension(u"Hz")},
        ms::Union{Int, AbstractVector{Int}} = 0:99,
        nmodes::Int = 10,
        ϵᵣ::Real = 1.0,
        epsr::Real = 1.0,
        tanδ::Real = 0.0,
        tandel::Real = 0.0,
        σ::Unitful.Quantity{<:Real, Unitful.dimension(u"S/m")} = Inf * u"S/m",
        sigma::Unitful.Quantity{<:Real, Unitful.dimension(u"S/m")} = Inf * u"S/m",
        Rq::Unitful.Quantity{<:Real, Unitful.dimension(u"m")} = 0.0u"m",
        col_fmts = ["%s", "%i", "%i", "%#.8g", "%8.4f", "%8.4f"],
        length_unit = Unitful.unit(a))

    freq_unit = Unitful.unit(f)
    ϵᵣ = max(ϵᵣ, epsr)
    tanδ = max(tanδ, tandel)
    σ = min(σ, sigma)

    modedata = cwg_modes(a, f; ms, nmodes, ϵᵣ, tanδ, σ, Rq)

    mdis = ms isa Integer ? [ms] : ms
    tit = "CWG: m ∈ $mdis, a = $a, f = $f"
    !isinf(σ) && (tit *= ", σ = $σ")
    !iszero(Rq) && (tit *= ", Rq = $Rq")
    !isone(ϵᵣ) && (tit *= ", ϵᵣ = $ϵᵣ")
    !iszero(tanδ) && (tit *= ", tanδ = $tanδ")
    display_mode_table(modedata, tit, freq_unit, length_unit, col_fmts)
    return nothing
end
