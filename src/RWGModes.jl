"""
    RWGMode <: HomogeneousMetalPipeMode

Struct representing an electromagnetic mode of a rectangular waveguide, homogeneously filled with dielectric.

## Fields:
- `p::TETM`: Mode type.
- `m::Int`: Index for modal variation along x direction.
- `n::Int`: Index for modal variation along y direction.
- `kco::Float64`: Cutoff wavenumber [radian/meter].
- `f::Float64`: Frequency [Hz].
- `γ::ComplexF64`: Attenuation constant [neper/meter].
- `Z::ComplexF64`: Modal wave impedance normalized to η₀, the impedance of free space [unitless].
"""
struct RWGMode <: HomogeneousMetalPipeMode
    p::TETM   # Mode type
    m::Int    # x index
    n::Int    # y index
    kco::Float64 # cutoff wavenumber [rad/meter]
    f::Float64 # Frequency [Hz]
    γ::ComplexF64 # Attenuation constant [np/meter]
    Z::ComplexF64 # Wave impedance normalized to η₀ [unitless]

    function RWGMode(
            p::TETM, m::Integer, n::Integer, kco::Real, f::Real, γ::Number, Z::Number)
        if p == TE
            (m ≥ 0 && n ≥ 0 && (m, n) ≠ (0, 0)) ||
                throw(ArgumentError("Illegal (m,n) = $((m,n)) for TE mode"))
        else
            (m > 0 && n > 0) || throw(ArgumentError("Illegal (m,n) = $((m,n)) for TM mode"))
        end
        return new(p, m, n, kco, f, γ, Z)
    end
end

"""
    RWGMode(; p, m, n, a, b, f=0, γ=0, Z=0)

Convenient keyword-argument constructor that defaults the frequency-dependent fields to zero and uses waveguide dimensions in meters.
"""
RWGMode(; p, m, n, a, b, f = 0, γ = 0, Z = 0) = RWGMode(
    p, m, n, rwgkco(a, b, m, n), f, γ, Z)

"""
    RWG <: HomogeneousMetallicWaveguide

A struct representing a section of metallic, rectangular waveguide, uniformly filled with dielectric.
It is assumed that the (possibly non-smooth) metal walls are either perfectly conducting or form a 
"good conductor", in the sense that the standard perturbational formulas for losses in waveguide walls
hold for higher-order modes.  For the dominant TE10 and the TE01 modes, accurate formulas from the recent
literature are used. See the docstring for function `rwgte10gz` for more information.

## Fields
- `a::Float64`: The x-dimension of the waveguide [m].
- `b::Float64`: The y-dimension of the waveguide [m].
- `l::Float64`: The length of the waveguide section (along z) [m].
- `ϵᵣ::Float64`: The dielectric constant for the material filling the guide [unitless].
- `tanδ::Float64`: The loss tangent for the material filling the guide [unitless].
- `σ::Float54`: The bulk conductivity of the waveguide metal walls [S/m].
- `Rq::Float64`: The metal wall RMS surface roughness [m].
- `modes::Vector{RWGMode}`: The modes treated in this guide.  If provided, the list 
  should be sorted in order of increasing cutoff frequency.

Both positional and keyword argument constructors for `RWG` are available.
"""
@kwdef struct RWG <: HomogeneousMetallicWaveguide
    a::Float64  # x dimension [meter]
    b::Float64  # y dimension [meter]
    l::Float64 = 0.0  # Length [meter]
    ϵᵣ::Float64 = 1.0    # Relative permittivity of dielectric filling waveguide
    tanδ::Float64 = 0.0 # Loss tangent of dielectric filling waveguide
    σ::Float64 = Inf # Wall bulk conductivity [S/m]
    Rq::Float64 = 0.0 # Wall RMS surface roughness [m]
    modes::Vector{RWGMode} = RWGMode[]
end

"""
RWG(wgspec::AbstractString; kwargs...)

Alternate constructor allowing use of a specification string for the waveguide size, followed by the 
other optional keyword arguments of the standard constructor.

## Required Positional Input Argument
- `wgspec::AbstractString`: A standard (EIA, RCSC, or IEC) abbreviation for a rectangular waveguide size.
  Examples include `"WR650"` (EIA), `"WG7"` (RCSC), and `"R22"` (IEC).

## Optional Keyword Arguments and their Default Values:
The remaining fields of the standard keyword constructor. They are:
- `l::Float64`=0.0: The length of the waveguide section (along z) [m].
- `ϵᵣ::Float64`=1.0: The dielectric constant for the material filling the guide [unitless].
- `tanδ::Float64=0.0`: The loss tangent for the material filling the guide [unitless].
- `σ::Float54=Inf`: The bulk conductivity of the waveguide metal walls [S/m].
- `Rq::Float64=0.0`: The metal wall RMS surface roughness [m].
- `modes::Vector{RWGMode}=RWGMode[]`: The modes treated in this guide.  If provided, the list 
  should be sorted in order of increasing cutoff frequency.
"""
function RWG(wgspec::AbstractString; kwargs...)
    (a, b) = lookup_rwg(wgspec)
    am = Unitful.ustrip(u"m", a)
    bm = Unitful.ustrip(u"m", b)
    return RWG(; a = am, b = bm, kwargs...)
end

"""
    setup_modes!(wg::RWG, f::Real, nmodes::Integer=length(c.modes); method=:auto)

Set up the modes for a uniform rectangular waveguide.

## Positional Arguments
- `wg`: A `RWG` instance.  If `wg.modes` is empty, then it will be appended to `length(n)`.  After this, or
  if it is already allocated, each mode in the vector will be replaced by one with updated values of `γ` and 
  `Z` corresponding to the frequency `f`.
- `f`: The frequency [Hz].
- `nmodes`: (optional) The number of modes to append to `wg.modes` if it is empty.  If `wg.modes` is nonempty, then
  `nmodes` must be equal to `length(wg.modes)`.
## Keyword Arguments
- `method::Symbol`: (optional)  Either `:auto` (the default) or `:ploss`.  If `:auto`, then 
  for the TE10 and TE01 modes, the updated propagation constant is computed using the formulas from Reference 
  [1], while for all other modes, the variational approach of Collin [2] pp. 350-354 is used.  If `:ploss`,
  then the attenuation for modes above cutoff is calculated using the power loss method (Table 5.2 of [2]).
  Note that the power loss method is not valid for degenerate modes (i.e. modes where both `m` and `n` are
  nonzero).  It is included here for reference only. `:auto` should be used when an accurate answer is desired.

## Return Value
- `wg`: The updated `RWG` instance is returned.  

## References
[1] Lomakin, Konstantin, Gerald Gold, and Klaus Helmreich. "Analytical waveguide model precisely
predicting loss and delay including surface roughness." IEEE Trans. Microwave Theory Tech., Vol 66,
no. 6 (2018): 2649-2662.

[2] R. E. Collin, **Field Theory of Guided Waves** 2cd Edition, IEEE Press, 1991, Table 5.2, p. 351.
"""
function setup_modes!(wg::RWG, f::Real, nmodes::Integer = length(wg.modes); method::Symbol = :auto)
    (; a, b, σ, Rq, ϵᵣ, tanδ) = wg
    if isempty(wg.modes)
        nmodes ≤ 0 && throw(ArgumentError("nmodes must be > 0"))
        rlattice = 4 * sqrt(nmodes / (a * b * π))
        mmax = ceil(Int, rlattice * a)
        nmax = ceil(Int, rlattice * b)
        modes = [RWGMode(; p, m, n, a, b)
                 for p in (TE, TM), m in 0:mmax, n in 0:nmax
                 if
                 (p == TE && (m, n) ≠ (0, 0)) || (p == TM && m > 0 && n > 0)]
        sort!(modes; by = x -> x.kco)
        deleteat!(modes, (nmodes + 1):length(modes))
        append!(wg.modes, modes)
    else
        nmodes == length(wg.modes) ||
            throw(ArgumentError("nmodes not equal to number of existing modes in wg"))
    end

    k₀ = 2π * f / c₀ # free-space wavenumber [rad/m]
    rootϵ = mysqrt(ϵᵣ * complex(1.0, -tanδ))
    k = k₀ * rootϵ # wavenumber in dielectric
    k² = k^2
    ηnorm = inv(rootϵ) # η/η₀
    η = η₀ * ηnorm # intrinsic impedance in dielectric
    Zs = Zsurface(f, σ, Rq, :normal) # surface impedance ([Ω/□])
    Rs = real(Zs) # surface resistance ([Ω/□])

    # Update γ and Z
    mnold = (0, 0)
    local γte, γtm # So we can reuse old values within the loop below
    for (q, mode) in pairs(wg.modes)
        (; m, n, p, kco) = mode
        if method == :auto
            # For TE10 or TE01, use Lomakin 2018
            if (p, m, n) == (TE, 1, 0)
                (γ, Z) = rwgte10gz(a, b, f; ϵᵣ, tanδ, σ, Rq)
                Z /= η₀
            elseif (p, m, n) == (TE, 0, 1)
                (γ, Z) = rwgte10gz(b, a, f; ϵᵣ, tanδ, σ, Rq)
                Z /= η₀
            else
                # For all other modes use Collin's variational approach:
                (m, n) ≠ mnold && ((; γte, γtm) = collin_gammas(a, b, m, n, k, Zs / η))
                γ = isTE(p) ? γte : γtm 
                β = -im * γ
                Z = isTE(p) ? (k / β * ηnorm) : (β / k * ηnorm)
            end
        elseif method == :ploss
            γ = mysqrt(kco^2 - k²)
            Rsnorm = Rs / η
            params = (; a, b, m, n, p, kco)
            αc = alphaploss(params, k, Rsnorm)
            γ += αc
            β = -im * γ
            Z = isTE(p) ? (k / β * ηnorm) : (β / k * ηnorm)
        else
            throw(ArgumentError("Unknown method: $method"))
        end
        wg.modes[q] = update(mode; f, γ, Z)
        mnold = (m, n)
    end
    return wg
end


"""
    alphaploss(params, k, Rsnorm)

Compute the attenuation constant of a rectangular waveguide using the classical power loss method.

## Input Arguments
- `params`: A NamedTuple or other destructurable object containing the following fields:
  - `a`, `b`: waveguide dimensions in x and y directions [m].
  - `m`, `n`, `p::TETM`: Mode indices and mode type.
  - `kco`: Cutoff wavenumber of mode [rad/m].
- `k`: Wavenumber in the dielectric filling the waveguide.
- `Rsnorm`: Ratio of surface resistance of the waveguide metal walls to the intrinsic impedance of the 
  dielectric filling the waveguide.

  ## Return Value
- `α`:  For frequencies above cutoff, the attenuation constant [np/m] calculated via the classical 
  power loss method.  For frequencies below cutoff, zero is returned.

## References:
1. R. E. Collin, **Field Theory of Guided Waves** 2cd Edition, IEEE Press, 1991, Table 5.2, p. 351.
2. J. Uher, J. Bornemann and U. Rosenberg, **Waveguide Components for Antenna Feed
   Systems: Theory and CAD**, Artech House, 1991, p. 114.
"""
function alphaploss(params, k::Number, Rsnorm::Real)
    (; a, b, m, n, p, kco) = params
    iszero(Rsnorm) && return 0.0
    ratio = real((kco / k)^2)
    ratio ≥ 1 && return 0.0
    boa = b / a
    α = 2 * real(Rsnorm) / (b * sqrt(1 - ratio))
    m²b² = (m * b)^2
    n²a² = (n * a)^2
    if p == TE
        ϵ₀ₙ = iszero(n) ? 1 : 2
        α *= (
            (1 + boa) * ratio +
            boa * (ϵ₀ₙ / 2 - ratio) * (m²b² / boa + n²a²) / (m²b² + n²a²)
        )
    else
        n²a³ = n²a² * a
        α *= (m²b² * b + n²a³) / (m²b² * a + n²a³)
    end
    return α
end

"""
    alphaploss(rwg::RWG, modeindex, k, Rsnorm)

Compute the attenuation constant of a rectangular waveguide using the classical power loss method.

## Input Arguments
- `rwg`: An `RWG` instance with `modes` instantiated.
- `modeindex`: The index into `rwg.modes` for the mode of interest.
- `k`: Wavenumber in the dielectric filling the waveguide.
- `Rsnorm`: Ratio of surface resistance of the waveguide metal walls to the intrinsic impedance of the 
  dielectric filling the waveguide.

## Return Value
- `α`:  For frequencies above cutoff, the attenuation constant [np/m] calculated via the classical 
  power loss method.  For frequencies below cutoff, zero.

## References:
1. R. E. Collin, **Field Theory of Guided Waves** 2cd Edition, IEEE Press, 1991, Table 5.2, p. 351.
2. J. Uher, J. Bornemann and U. Rosenberg, **Waveguide Components for Antenna Feed
   Systems: Theory and CAD**, Artech House, 1991, p. 114.
"""
function alphaploss(rwg::RWG, modeindex::Integer, k, Rsnorm)
    (; a, b) = rwg
    (; m, n, p, kco) = rwg.modes[modeindex]
    params = (; a, b, m, n, p, kco)
    return alphaploss(params, k, Rsnorm)
end



"""
    rwgkco(a, b, m, n)

Compute the cutoff frequency for a rectangular waveguide.

## Input Variables
- `a`, `b`: The dimensions of the waveguide in the x and y directions, respectively.
- `m`, `n`: Mode indices in the x and y directions, respectively.

# Return value:
`kco`: Cutoff wavenumber in inverse units of `a` and `b`.
"""
rwgkco(a::Real, b::Real, m::Integer, n::Integer) = π * hypot(m / a, n / b)

"""
    collin_gamma(a, b, m, n, p, k, Zsnorm)

Compute variational estimate of the complex attenuation constant for a rectangular waveguide mode, 
using the variational method described on pages 350-354 of the reference.

## Arguments
- `a`, `b`: Waveguide dimensions along x and y, respectively [m].
- `m`, `n`: Nonnegative integer mode numbers along x and y, respectively, not both zero, and neither zero if `p == TM`.
- `p::TETM`: Mode type (TE or TM).
- `k`: The wavenumber in the dielectric medium filling the waveguide.
- `Zsnorm`: The surface impedance of the metal waveguide walls normalized to the intrinsic impedance of the dielectric medium.

## Return Value
`γ`: A variational estimate of the complex attenuation constant [Np/m] of the specified mode.
"""
function collin_gamma(a::Real, b::Real, m::Integer, n::Integer, p::TETM, k, Zsnorm)
    zeroindex = any(iszero, (m, n))
    isTM(p) && zeroindex && throw(ArgumentError("Zero mode index not allowed for TM mode"))

    (; γte, γtm) = collin_gammas(a::Real, b::Real, m::Integer, n::Integer, k, Zsnorm)
    if isTE(p)
        return γte
    else
        return γtm
    end
end

"""
    collin_gammas(a, b, m, n, k, Zsnorm)

Compute variational estimate of the complex attenuation constant for a rectangular waveguide mode, 
using the variational method described on pages 350-354 of the reference.

## Arguments
- `a`, `b`: Waveguide dimensions along x and y, respectively [m].
- `m`, `n`: Nonnegative integer mode numbers along x and y, respectively, not both zero, and neither zero if `p == TM`.
- `k`: The wavenumber in the dielectric medium filling the waveguide.
- `Zsnorm`: The surface impedance of the metal waveguide walls normalized to the intrinsic impedance of the dielectric medium.

## Return Value
`(; γte, γtm)`:  TE and TM complex attenuation constants [Np/m].
"""
function collin_gammas(a::Real, b::Real, m::Integer, n::Integer, k, Zsnorm)
    (m, n) ≠ (0, 0) || throw(ArgumentError("Both mode indices are zero"))
    zeroindex = any(iszero, (m, n))
    kco = rwgkco(a, b, m, n)

    # Compute quantities needed in objective function
    ϵ₀ₘ = iszero(m) ? 1 : 2
    ϵ₀ₙ = iszero(n) ? 1 : 2
    k² = complex(k^2)
    kco² = kco^2
    β₀² = k² - kco²
    β₀ = sqrt(β₀²)
    mπoa = m * π / a
    mπoa² = mπoa^2
    nπob = n * π / b
    nπob² = nπob^2
    # Constants in Eq. (66a) of reference (multiplied by a to make all terms unitless):
    f11 = im * k * ((a^2 * b) / (ϵ₀ₘ * ϵ₀ₙ)) / Zsnorm
    term11b = a * kco² * (2a / ϵ₀ₘ + 2b / ϵ₀ₙ) + a * β₀² / kco² * (b * nπob² + a * mπoa²)
    A12 = k² * (a / kco² * mπoa * nπob * (b - a))
    A21 = a * β₀² / kco² * mπoa * nπob * (b - a)
    term22b = k² * (a / kco² * (a * nπob² + b * mπoa²))

    # For TE₀ₙ or TEₘ₀ modes, solve (66a) in closed form:
    if zeroindex
        γ² = -term11b / f11 - β₀²
        γte = mysqrt(γ²)
        return (; γte, γtm = NaN)
    end

    parameters = (; a, β = β₀, β² = β₀², f11, term11b, A12, A21, term22b)

    # Initial brute search along a line in the complex Δγa plain to locate two minima:
    αpl = alphaploss((; a, b, p = TE, m, n, kco), k, real(Zsnorm))
    if real(k) > kco
        nsamples = 500
        αa_max = 3 * αpl * a
        αas = range(0, αa_max, length = nsamples)
        i1 = i2 = 0 # Initialize the locations of minima in αas
        dmag2_im1 = dmag2_im2 = -Inf # Initialize previous two results
        for (i, αa) in pairs(αas)
            Δβa = -αa
            γ = αa / a + im * (β₀ + Δβa / a)
            determ = det(collinmatrix(γ, parameters))
            dmag2 = abs2(determ)
            if dmag2_im2 > dmag2_im1 < dmag2
                if iszero(i1)
                    i1 = i - 1
                else
                    i2 = i - 1
                    break
                end
            end
            dmag2_im1, dmag2_im2 = dmag2, dmag2_im1 # Update prev. results
        end
        any(iszero, (i1, i2)) && error("Got a zero: i1 = $i1, i2 = $i2")
    else
        # At or below cutoff.  Can only find a single solution here (typ. TM)
        i1 = i2 = 1
    end
    # Now polish the minima
    αa_spread = αas[i2] - αas[i1]
    iszero(αa_spread) && (αa_spread = 1e-4) # At or below cutoff
    local data1, data2
    for i in (i1, i2)
        αa = αas[i]
        Δβa = -αa
        Δγa₀ = complex(αa, Δβa)
        x0 = [real(Δγa₀), imag(Δγa₀)]
        Δ = 0.2 * αa_spread
        xu = x0 + [Δ, Δ]
        xl = x0 - [Δ, Δ]

        rhobeg = 0.9 * Δ
        rhoend = 1e-10 * rhobeg
        x, info = bobyqa(x0; xl, xu, rhobeg, rhoend) do x
            collin_gamma_objective(x, parameters).objective
        end
        (; γ, A) = collin_gamma_objective(x, parameters)
        A /= norm(A)
        U, S, V = svd(A)
        if i == i1
            data1 = (; x, γ, A, U, S, V, info)
        else
            data2 = (; x, γ, A, U, S, V, info)
        end
    end
    norm(data1.A * data1.V[:, end]) < 1e-2 || @warn "Large residual for data1" data1
    norm(data2.A * data2.V[:, end]) < 1e-2 || @warn "Large residual for data2" data2
    γtm, γte = real(data1.γ) < real(data2.γ) ? (data1.γ, data2.γ) : (data2.γ, data1.γ)
    return (; γte, γtm)
end

#penalty_fun(x, x0 = 0.9, ϵ = 0.01) = log1p(exp((x - x0) / ϵ)) * ϵ
function penalty_fun(x; a = 100, n = 2, x0 = 0.9, ϵ = 0.01)
    xⁿ = x^n
    if xⁿ > 5
        result = a * (xⁿ - x0)
    else
        result = a * log1p(exp((xⁿ - x0) / ϵ)) * ϵ
    end
    return result
end

"""
    collinmatrix(γ, parameters)

Compute the 2x2 matrix from Eq. (66) on page 353 of Collin (after each row multiplied by `a` to render unitless).

# Arguments
- `γ`: The complex attenuation constant [Np/m] from which the matrix is to be computed.
- `parameters`: A named tuple consisting of at least the following fields:
  - `β²`: The square of `β`.
  - `f11`, `term11b`, `A12`, `A21`, `term22b`: Constants needed to evaluate the matrix
    elements defined in Collin Eq. (66) on page 353.  These have all been multiplied by 
    waveguide dimension `a` to render them unitless.

## Return Value
- `detabs2`:  The magnitude squared of the determinant of the matrix equation in Eq. (66), 
  after each row has been multiplied by `a` to make the matrix entries unitless.
"""
function collinmatrix(γ, parameters)
    (; β², f11, term11b, A12, A21, term22b) = parameters
    term11a = im * (γ^2 + β²) * f11
    A11 = term11a + term11b
    A22 = term11a + term22b
    A = @SMatrix [A11 A12; A21 A22]
    return A
end

"""
    collin_gamma_objective(Δγari, parameters)

Objective function for `collin_gamma`.

# Arguments
- `Δγari`: A 2-vector containing real and imaginary parts of `Δγa`. The latter is defined
  such that the current estimate for `γ` is obtained as `im*β + Δγa / a`.  Note that
  `Δγa` is unitless, and is also very small for a waveguide with good conductor walls.
- `parameters`: A named tuple consisting of at least the following fields:
  - `a`: Waveguide x dimension [m].
  - `p::TETM`: Desired mode type
  - `β`: Phase constant of the degenerate modes assuming PEC walls [rad/m].
  - `β²`: The square of `β`.
  - `f11`, `term11b`, `A12`, `A21`, `term22b`: Constants needed to evaluate the matrix
    elements defined in Collin Eq. (66) on page 353.  These have all been multiplied by 
    waveguide dimension `a` to render them unitless.

## Return Value
- # TBC
  
"""
function collin_gamma_objective(Δγari, parameters)
    (; a, β) = parameters
    Δγa = complex(Δγari[1], Δγari[2])
    γ = im * β + Δγa / a
    A = collinmatrix(γ, parameters)
    Anorm = A / norm(A)
    determ = det(Anorm)
    temag = abs(Anorm[1, 2])
    tmmag = abs(Anorm[1, 1])
    #penalty = penalty_fun(ratio)
    penalty = 0.0
    detobjective = abs2(determ) / norm(A)
    objective = detobjective + penalty
    return (; objective, Δγa, γ, A, penalty)
end

"""
    rwgte10gz(a, b, f; ϵᵣ, tanδ, σ, Rq) -> (γ, Z)

Accurately compute γ and Z (prop. constant and wave impedance) of TE10 mode in rectangular guide with rough walls.

## Required Positional Arguments
- `a`,  `b`: The waveguide x and y dimensions, `Unitful` quantities with dimensions of length.
- `f`: Frequency, a `Unitful` quantity with dimensions of inverse time.
## Optional Keyword Arguments
- `ϵᵣ`, `tanδ`: Dielectric constant and loss tangent, respectively, of material filling the waveguide. Default
  values are `1.0` and `0.0`, respectively.
- `σ`: The bulk conductivity of the waveguide metal walls, a Unitful quantity with the same dimensions as `u"S/m"`.
  Default value is `Inf*u"S/m"`.
- `Rq`: The RMS surface roughness of the waveguide walls, a `Unitful` quantity with dimension of length. Default
  value is `0.0u"m"`

## Return Values
- `γ`: Complex propagation constant [neper/m]
- `Z`: Complex wave impedance [Ω]

## References 
- [1] Lomakin, Konstantin, Gerald Gold, and Klaus Helmreich. "Analytical waveguide model precisely
  predicting loss and delay including surface roughness." IEEE Transactions on Microwave Theory
  and Techniques 66, no. 6 (2018): 2649-2662.
- [2] Gold, Gerald, and Klaus Helmreich. "A physical surface roughness model and its applications."
  IEEE Transactions on Microwave Theory and Techniques 65, no. 10 (2017): 3720-3732.
- [3] Lomakin, Konstantin, Gerald Gold, and Klaus Helmreich. "Transmission line model for rectangular 
  waveguides accurately incorporating loss effects." In 2017 IEEE 21st Workshop on Signal and Power 
  Integrity (SPI), pp. 1-4. IEEE, 2017.
"""
function rwgte10gz(
        a::Unitful.Quantity{<:Real, Unitful.dimension(u"m")},
        b::Unitful.Quantity{<:Real, Unitful.dimension(u"m")},
        f::Unitful.Quantity{<:Real, Unitful.dimension(u"Hz")};
        ϵᵣ::Real = 1.0,
        tanδ::Real = 0.0,
        σ::Unitful.Quantity{<:Real, Unitful.dimension(u"S/m")} = Inf * u"S/m",
        Rq::Unitful.Quantity{<:Real, Unitful.dimension(u"m")} = 0.0u"m")
    am, bm, Rqm = map(x -> Unitful.ustrip(u"m", x), (a, b, Rq))
    fhz = Unitful.ustrip(u"Hz", f)
    σspm = Unitful.ustrip(u"S/m", σ)
    return rwgte10gz(am, bm, fhz; ϵᵣ, tanδ, σ = σspm, Rq = Rqm)
end

"""
    rwgte10gz(wgspec, f; kwargs...) -> (γ, Z)

This method accepts a rectangular waveguide "specification string" instead of a and b dimensions.

## Required Positional Inputs
- `wgspec::AbstractString`: A standard (EIA, RCSC, or IEC) abbreviation for a rectangular waveguide size.
  Examples include `"WR650"` (EIA), `"WG7"` (RCSC), and `"R22"` (IEC).

## Optional Keyword Arguments
- `ϵᵣ`, `tanδ`: Dielectric constant and loss tangent, respectively, of material filling the waveguide. Default
  values are `1.0` and `0.0`, respectively.
- `σ`: The bulk conductivity of the waveguide metal walls, a Unitful quantity with the same dimensions as `u"S/m"`.
  Default value is `Inf*u"S/m"`.
- `Rq`: The RMS surface roughness of the waveguide walls, a `Unitful` quantity with dimension of length. Default
  value is `0.0u"m"`

## Return Values
- `γ`: Complex propagation constant [neper/m]
- `Z`: Complex wave impedance [Ω]
"""
function rwgte10gz(
        wgspec::AbstractString, f::Unitful.Quantity{<:Real, Unitful.dimension(u"Hz")}; kwargs...)
    (a, b) = lookup_rwg(wgspec)
    return rwgte10gz(a, b, f; kwargs...)
end

"""
    rwgte10gz(a::Real, b::Real, f::Real; ϵᵣ::Real=1.0, tanδ::Real=0.0, σ::Real=Inf, Rq::Real=0.0) -> (γ, Z)
    rwgte10gz(a::Real, b::Real, f::Real; epsr::Real=1.0, tandel::Real=0.0, sigma::Real=Inf, Rq::Real=0.0) -> (γ, Z)

This method accepts SI values as pure numbers (without attached units). 
## Required Positional Arguments
- `a`,  `b`: The waveguide x and y dimensions [m].
- `f`: Frequency [Hz]
## Optional Keyword Arguments
- `ϵᵣ` or `epsr`: Dielectric constant of material filling the waveguide. Default value is `1.0`.
- `tanδ` or `tandel`: Loss tangent of material filling the waveguide. Default value is `0.0`. 
- `σ` or `sigma`: The bulk conductivity of the waveguide metal walls [S/m]. Default value is `Inf`.
- `Rq`: The RMS surface roughness of the waveguide walls [m]. Default value is `0.0`.

## Return Value
The tuple `(γ, Z)` where:
- `γ`: Complex propagation constant [neper/m]
- `Z`: Complex wave impedance [Ω]
"""
function rwgte10gz(
        a::Real,
        b::Real,
        f::Real;
        ϵᵣ::Real = 1.0,
        epsr::Real = 1.0,
        tanδ::Real = 0.0,
        tandel::Real = 0.0,
        σ::Real = Inf,
        sigma::Real = Inf,
        Rq::Real = 0.0)
    ϵᵣ = max(ϵᵣ, epsr)
    tanδ = max(tanδ, tandel)
    σ = min(σ, sigma)

    isinf(σ) &&
        !iszero(Rq) &&
        throw(ArgumentError("Rq must be zero for infinite conductivity"))

    w, h = a, b # Use nomenclature from the reference
    ω = 2π * f
    Zs = Zsurface(f, σ, Rq)
    σeff = effective_conductivity(f, σ, Rq)

    Lₒ′ = μ₀ # Eq. (14) of [3]
    C′ = ϵ₀ * ϵᵣ  # Eq. (15) of [3]
    Lₒ′′ = μ₀ * (w / π)^2 # Eq. (16) of [3]
    G′ = ω * C′ * tanδ # Eq. (35) of [3]

    if isinf(σ)
        R′ = R′′ = Lᵢ′ = Lᵢ′′ = 0.0
    else
        δc = inv(real(Zs) * σeff) # Eq. (37) of [2], same as Eq. (35) of [1]
        δm = inv(imag(Zs) * σ) # Eq. (37) of [2], same as Eq. (36) of [1]
        R′ = 2 / (σeff * δc * h) # Eq. (37) of [1]
        R′′ = 2w * (w + 2h) / (π^2 * h * σeff * δc) # Eq. (38) of [1]
        Lᵢ′ = 2 / (σ * δm * h * ω) # Eq. (39) of [1]
        Lᵢ′′ = 2w * (w + 2h) / (π^2 * h * σ * δm * ω) # Eq. (40) of [1]
    end

    L′ = Lₒ′ + Lᵢ′ # Eq. (5) of [1]
    L′′ = Lₒ′′ + Lᵢ′′ # Eq. (6) of [1]

    # Eqs. (12) and (13) of [1]:
    Zckt = complex(R′, ω * L′)
    Yckt = inv(complex(R′′, ω * L′′)) + complex(G′, ω * C′)
    γ = mysqrt(Zckt * Yckt)
    Z = sqrt(Zckt / Yckt)
    return (γ, Z)
end

"""
    rwggz(p, m, n, a, b, f; ϵᵣ=1.0, tanδ=0.0, σ=Inf, Rq=0.0) -> (γ, Z)
    rwggz(p, m, n, a, b, f; epsr=1.0, tandel=0.0, sigma=Inf, Rq=0.0) -> (γ, Z)

Accurately compute γ and Z (prop. constant and wave impedance) of a mode in rectangular guide with rough walls.

This method accepts SI values as pure numbers (without attached units). 
## Required Positional Arguments
- `p::TETM`: The mode type.
- `m::Int`, `n::Int`: Mode indices in the x and y directions, resp.
- `a::Real`,  `b::Real`: The waveguide x and y dimensions [m].
- `f`: Frequency [Hz].
## Optional Keyword Arguments
- `ϵᵣ` or `epsr`: Dielectric constant of material filling the waveguide. Default value is `1.0`.
- `tanδ` or `tandel`: Loss tangent of material filling the waveguide. Default value is `0.0`. 
- `σ` or `sigma`: The bulk conductivity of the waveguide metal walls [S/m]. Default value is `Inf`.
- `Rq`: The RMS surface roughness of the waveguide walls [m]. Default value is `0.0`.

## Return Value
The tuple `(γ, Z)` where:
- `γ`: Complex propagation constant [neper/m]
- `Z`: Complex wave impedance [Ω]

## References:
- [1] Yeap, Kim Ho et al., "Attenuation in circular and rectangular waveguides." 
  Electromagnetics 37, no. 3 (2017): 171-184.
"""
function rwggz(
        p::TETM,
        m::Int,
        n::Int,
        a::Real,
        b::Real,
        f::Real;
        ϵᵣ::Real = 1.0,
        epsr::Real = 1.0,
        tanδ::Real = 0.0,
        tandel::Real = 0.0,
        σ::Real = Inf,
        sigma::Real = Inf,
        Rq::Real = 0.0)
    ϵᵣ = max(ϵᵣ, epsr)
    tanδ = max(tanδ, tandel)
    σ = min(σ, sigma)

    isinf(σ) &&
        !iszero(Rq) &&
        throw(ArgumentError("Rq must be zero for infinite conductivity"))

    ω = 2π * f
    k₀ = 2π * f / c₀ # free-space wavenumber [rad/m]
    rootϵ = mysqrt(ϵᵣ * complex(1.0, -tanδ))
    k = k₀ * rootϵ # wavenumber in dielectric
    k² = k^2
    ηnorm = inv(rootϵ) # η/η₀
    η = η₀ * ηnorm # intrinsic impedance in dielectric
    Zs = Zsurface(f, σ, Rq)
    params = (; m, n, k, a, b)
    Δs0 = @SVector [0.0, 0.0, 0.0, 0.0]
    prob = snls.SimpleNonLinearProblem(nlsovekdxakdyb, Δs0, params)
    sol = snls.solve(prob, Snls.SimpleNewtonRaphson)
    snls.SciMLBase.successful_retcode(sol) || error("Failed nonlinear solve")
    Δs = sol.u
    Δkxa = complex(Δs[1], Δs[2])
    Δkyb = complex(Δs[3], Δs[4])
    kx = (m * π + Δkxa) / a
    ky = (n * π + Δkyb) / b
    kρ² = kx^2 + ky^2
    kz = sqrt(k^2 - kρ²)
    γ = im * kz
    if isTE(p)
        Z = im * k * η / γ
    else
        Z = η * γ / (im * k)
    end

    return (γ, Z)
end

"""
    nlsolvedkxadkyb(Δs, params) -> residual::SVector{4, Float64}

Function defining a NonLinearProblem whose solution provides rigorously correct kx and ky values for a rectangular waveguide.

## Positional Arguments
- `Δs`: A 4-vector containing the unitless values `[real(Δkxa), imag(Δkxa), real(Δkyb), imag(Δkyb)].  Here `Δkxa` is the departure
  of `kx` from its nominal value (for PEC walls) of `m*π/a`.  Similarly, `Δkyb` is the departure of `ky` from its nominal value
  of `n*π/b`.  The entries in `Δs` are to be adjusted to produce a zero return vector.
- `params`:  A named tuple containing the fixed parameters of the waveguide problem.  These are
  - `m` and `n`: Mode indices in the x and y directions, respectively.
  - `k`: The (complex) wavenumber in the dielectric filling the waveguide [radians/m].
  - `a` and  `b`:  The waveguide (inner) dimensions in the x and y directions, respectively [m].
  - `ηwn`: Waveguide wall surface impedance normalized to the intrinsic impedance of the dielectric filling the waveguide.

## Return Value
- `residual`: A static 4-vector containing the real and imaginary parts of the two expressions that must be zeroed for a valid 
  solution.

## Reference: 

"""
function nlsolvedkxadkyb(Δs::AbstractVector, params)
    Δkxa = complex(Δs[1], Δs[2])
    Δkyb = complex(Δs[3], Δs[4])
    (; m, n, k, a, b, ηwn) = params
    kxa = m * π + Δkxa
    kx = kxa / a
    kyb = n * π + Δkyb
    ky = kyb / b
    kρ² = kx^2 + ky^2
    sux, cux = sincos((kxa + m * π) / 2)
    suy, cuy = sincos((kyb + n * π) / 2)
    kz = sqrt(k^2 - kρ²)

    eq1 = (im * k * ky * suy / kρ² + ηwn * cuy) * (im * ηwn * k * ky * cuy / kρ² - suy) -
          (kz * kx / kρ²)^2 * ηwn * cuy * suy
    eq2 = (im * k * kx * sux / kρ² + ηwn * cux) * (im * ηwn * k * kx * cux / kρ² - sux) -
          (kz * ky / kρ²)^2 * ηwn * cux * sux
    return @SVector [real(eq1), imag(eq1), real(eq2), imag(eq2)]
end

"""
    nlsolvekxakyb(rikakbs, params) -> residual::SVector{4, Float64}

Function defining a NonLinearProblem whose zero solves the determinantal equations for the
kx and ky values of a rectangular waveguide.

## Positional Arguments
- `rikakbs`: A 4-vector containing the unitless values `[real(kxa), imag(kxa), real(kyb), imag(kyb)].  Here `kxa` is the product
  of `kx` and `a` (the waveguide x dimension) Similarly, `kyb` is the product of `ky` and `b` (the waveguide y dimension).
  The entries in `rikakbs` are to be adjusted to produce a zero return vector.
- `params`:  A named tuple containing the fixed parameters of the waveguide problem.  These are
  - `m` and `n`: Mode indices in the x and y directions, respectively.
  - `k`: The (complex) wavenumber in the dielectric filling the waveguide [radians/m].
  - `a` and  `b`:  The waveguide (inner) dimensions in the x and y directions, respectively [m].
  - `ηwn`: Waveguide wall surface impedance normalized to the intrinsic impedance of the dielectric filling the waveguide.

## Return Value
- `residual`: A static 4-vector containing the real and imaginary parts of the two expressions that must be zeroed for a valid 
  solution.

## Reference: 
[1] Yeap, Kim Ho, et al. "Attenuation in circular and rectangular waveguides." **Electromagnetics**
37.3 (2017): 171-184.  See Eqs. (6a) and (6b).
[2] P. Simon, "Notes on Attenuation in Circular and Rectangular Waveguides", Eqs (7a) and (7b).

"""
function nlsolvekxakyb(rikakbs::AbstractVector, params)
    kxa = complex(rikakbs[1], rikakbs[2])
    kyb = complex(rikakbs[3], rikakbs[4])
    (; m, n, k, a, b, ηwn) = params
    kx = kxa / a
    ky = kyb / b
    kρ² = kx^2 + ky^2
    sux, cux = sincos((kxa + m * π) / 2)
    suy, cuy = sincos((kyb + n * π) / 2)
    kz = sqrt(k^2 - kρ²)

    eq1 = (im * k * ky * suy / kρ² + ηwn * cuy) * (im * ηwn * k * ky * cuy / kρ² - suy) -
          (kz * kx / kρ²)^2 * ηwn * cuy * suy
    eq2 = (im * k * kx * sux / kρ² + ηwn * cux) * (im * ηwn * k * kx * cux / kρ² - sux) -
          (kz * ky / kρ²)^2 * ηwn * cux * sux
    #eq1 *= 1e8
    #eq2 *= 1e8
    return @SVector [real(eq1), imag(eq1), real(eq2), imag(eq2)]
end

"""
    e0h0mat(kxa, kyb, params)

Compute the dispersion matrix `A` whose null space defines the possible z-component field coefficients.

## Input Arguments
- `kxa`, `kyb`: The (unitless) product of proposed x and y modal wavenumbers times the x and y dimensions of the
  waveguide, resp.
- `params`: A struct or named tuple containing the following parameters of the waveguide problem:
  - `m`, `n`: Mode numbers in the x and y directions, resp.
  - `a`, `b`: Waveguide x and y dimensions, resp. [m]
  - `k`: Wavenumber [radians/m] in the dielectric filling the waveguide.
  - `ηwn`: The wall surface impedance normalized to the intrinsic impedance of the dielectric filling the waveguide.

## Return Value
- `A`: `ComplexF64` matrix of dimensions 4×2 computed using the provided inputs, whose nullspace spans
  the possible modal vectors `[E₀, η*H₀]`.  Here `E₀` is the coefficient of the z-component of electric field,
  `H₀` is the coefficient of the z-component of magnetic field, and η is the intrinsic impedance of the
  dielectric filling the waveguide.
"""
function e0h0mat(kxa, kyb, params)
    # Set up matrix Eq. (11) of notes
    (; m, n, k, a, b, ηwn) = params
    kx = kxa / a
    ky = kyb / b
    kρ² = kx^2 + ky^2
    ux = (m * π + kxa) / 2
    uy = (n * π + kyb) / 2
    sux, cux = sincos(ux)
    suy, cuy = sincos(uy)
    kz = sqrt(k^2 - kρ²)
    hyz = ky * kz / kρ²
    hx0 = kx * k / kρ²
    hxz = kx * kz / kρ²
    hy0 = ky * k / kρ²

    mat = zeros(ComplexF64, 4, 2)
    mat[1, 1] = -im * ηwn * hx0 * cux + sux
    mat[2, 1] = im * hyz * sux
    mat[3, 1] = im * ηwn * hy0 * cuy - suy
    mat[4, 1] = im * hxz * suy

    mat[1, 2] = im * ηwn * hyz * cux
    mat[2, 2] = im * hx0 * sux + ηwn * cux
    mat[3, 2] = im * ηwn * hxz * cuy
    mat[4, 2] = -im * hy0 * suy - ηwn * cuy

    return mat
end

"""
    compute_As_kxakybs(FGHz, rwg, mode_params; iterate=false)

    Compute dispersion matrices for a rectangular waveguide with lossy walls using 4 choices of (δx, δy)

## Input Arguments
- `FGHz`: The desired analysis frequency in GHz.
- `rwg::RWG`: An `RWG` instance with at least fields `a`, `b`, `σ`, and `Rq` set.
- `mode_params`: Named tuple or struct with fields `m::Integer`, `n::Integer`, `p::TETM` defining the desired mode.
- `iterate`: If `false` (default), then the `δx` and `δy` values computed from [1] are used directly to compute kx and ky and thus γ.
  If true, then these values are used as starting values in a root-finding algorithm implementing the numerical search described in
  [2].

## Return value
`(As, kxakybs)`, where

- `As`: A 4-vector containing four complex matrices of size (4,2) containing the coefficients of E0 and η*H0 in the modal 
  dispersion relation.  `As[i]` contains the matrix obtained by choosing (δx, δy) to be
  - `i==1`: (δx1, δy1)
  - `i==2`: (δx2, δy1)
  - `i==1`: (δx1, δy2)
  - `i==1`: (δx2, δy2)

- `kxakybs`: A 4-vector containing tuples `(kx*a, ky*b)` computed for the same four cases as `As`.

## References
- [1] Yeap, Kim Ho, et al. "Attenuation in circular and rectangular waveguides." Electromagnetics 37.3 (2017): 171-184.
- [2] Yeap, Kim Ho, et al. "Attenuation in rectangular waveguides with finite conductivity walls." Radioengineering 20.2 (2011).
"""
function compute_As_kxakybs(FGHz::Real, rwg::RWG, mode_params; iterate = false)
    (; a, b, σ, Rq) = rwg
    As = [zeros(ComplexF64, 4, 2) for _ in 1:4]
    kxakybs = [(zero(ComplexF64), zero(ComplexF64)) for _ in 1:4]
    (; m, n, p) = mode_params
    fhz = FGHz * 1e9
    setup_modes!(rwg, fhz, 10)
    k = zero(ComplexF64) + 2π * fhz / c₀
    ηwn = Zsurface(fhz, σ, Rq) / η₀
    if isTE(p)
        acm = a * 100
        bcm = b * 100
        ΔTx = im * (2.62e-6 / acm^2)
        ΔTy = im * (5.2e-7 / bcm^2)
    else
        ΔTx = ΔTy = zero(ComplexF64)
    end
    if m > 0
        ka = k * a
        δx1 = im * ηwn * ka / (m * π) + ΔTx
        δx2 = im * ηwn * (m * π) / ka + ΔTx
    else
        δx1 = δx2 = ΔTx
    end

    if n > 0
        kb = k * b
        δy1 = im * ηwn * kb / (n * π) + ΔTy
        δy2 = im * ηwn * (n * π) / kb + ΔTy
    else
        δy1 = δy2 = ΔTy
    end
    params = (; a, b, m, n, k, ηwn)
    modei = findfirst(1:length(rwg.modes)) do i
        mode = rwg.modes[i]
        mode.p == p && mode.m == m && mode.n == n
    end
    isnothing(modei) && error("Mode not found in rwg.modes")

    for (case, δs) in enumerate(((δx1, δy1), (δx2, δy1), (δx1, δy2), (δx2, δy2)))
        δx, δy = δs
        Δkxa = 2 * δx
        Δkyb = 2 * δy
        Δs0 = @SVector[real(Δkxa), imag(Δkxa), real(Δkyb), imag(Δkyb)]
        if iterate
            prob = snls.NonlinearProblem(nlsolvedkxadkyb, Δs0, params)
            sol = snls.solve(
                prob, snls.SimpleNewtonRaphson(autodiff = snls.AutoFiniteDiff()))
            snls.SciMLBase.successful_retcode(sol) ||
                error("Failed nonlinear solve for case $case")
            Δs = sol.u
        else
            Δs = Δs0
        end
        Δkxa = complex(Δs[1], Δs[2])
        Δkyb = complex(Δs[3], Δs[4])
        kxa = m * π + Δkxa
        kyb = n * π + Δkyb
        kxakybs[case] = (kxa, kyb)
        A = WaveguideModes.e0h0mat(kxa, kyb, params)
        As[case] .= A
    end
    return As, kxakybs
end

"""
    rwg_modes(a, b, f; nmodes=10, ϵᵣ=1.0, tanδ=0.0, σ=Inf*u"S/m", Rq=0.0u"m") -> modedata

Compute cutoff frequency, guide wavelength, and attenuation constant for first few modes of a rectangular waveguide.

Note: This function is intended for programmatic use.  For interactive use, see `rwg_modetable`.  


## Required Positional Arguments
- `a`, `b`: The waveguide inner dimensions as any `Unitful` length quantity, e.g. `b=0.8128u"cm"`, `b=320u"mil"`, 
  or `b=0.320u"inch"`.  `a` and `b` must be expressed in the same length units.
- `f`: The frequency as any `Unitful` frequency quantity, e.g. `f=3u"GHz"`,  `f=3e9u"Hz"` or `f=3000u"MHz"`.

## Optional Keyword Arguments
- `nmodes`: The number of modes to be treated.  Defaults to `10`.
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
  indices in the x and y directions, resp., `fco` is the cutoff frequency expressed in the same units as
  input argument `f` (but without attached units), `λg` is the guide wavelength expressed in the same 
  units as `a` and `b` (but without attached units), and `α` is the mode attenuation constant in units 
  of dB/unitlength (without attached units), where unitlength is one unit of length in the same units 
  as `a` and `b`. The modes are listed in order of increasing cutoff frequency. For cutoff modes (other
  than TE10 or TE01 which are calculated more accurately) the value of `λg` will be `Inf`.
"""
function rwg_modes(
        a::Unitful.Quantity{<:Real, Unitful.dimension(u"m")},
        b::Unitful.Quantity{<:Real, Unitful.dimension(u"m")},
        f::Unitful.Quantity{<:Real, Unitful.dimension(u"Hz")};
        nmodes::Int = 10,
        ϵᵣ::Real = 1.0,
        epsr::Real = 1.0,
        tanδ::Real = 0.0,
        tandel::Real = 0.0,
        σ::Unitful.Quantity{<:Real, Unitful.dimension(u"S/m")} = Inf * u"S/m",
        sigma::Unitful.Quantity{<:Real, Unitful.dimension(u"S/m")} = Inf * u"S/m",
        Rq::Unitful.Quantity{<:Real, Unitful.dimension(u"m")} = 0.0u"m")
    ϵᵣ = max(ϵᵣ, epsr)
    tanδ = max(tanδ, tandel)
    σ = min(σ, sigma)

    # Verify inputs
    (a > 0u"m" && b > 0u"m") || throw(ArgumentError("a and b must be positive"))
    length_unit = Unitful.unit(a)
    length_unit == Unitful.unit(b) ||
        throw(ArgumentError("a and b must have the same units"))

    f > 0u"Hz" || throw(ArgumentError("f must be positive"))
    freq_unit = Unitful.unit(f)
    nmodes > 0 || throw(ArgumentError("nmodes must be positive"))
    ϵᵣ ≥ 1 || throw(ArgumentError("ϵᵣ must be greater than or equal to 1"))
    tanδ ≥ 0 || throw(ArgumentError("tanδ must be nonnegative"))
    σ > 0u"S/m" || throw(ArgumentError("σ must be posiive"))
    Rq ≥ 0u"m" || throw(ArgumentError("Rq must be nonnegative"))

    # Convert to Float64 in standard MKS units and remove units. 
    am = Unitful.ustrip(Float64, u"m", a)
    bm = Unitful.ustrip(Float64, u"m", b)
    fHz = Unitful.ustrip(Float64, u"Hz", f)
    σSpm = Unitful.ustrip(Float64, u"S/m", σ)
    Rqm = Unitful.ustrip(Float64, u"m", Rq)

    # Create list of modes to consider:
    rlattice = 4 * sqrt(nmodes / (am * bm * π))
    mmax = ceil(Int, rlattice * am)
    nmax = ceil(Int, rlattice * bm)
    modes = [RWGMode(; p, m, n, a = am, b = bm)
             for p in (TE, TM), m in 0:mmax, n in 0:nmax
             if
             (p == TE && (m, n) ≠ (0, 0)) || (p == TM && m > 0 && n > 0)]
    sort!(modes; by = x -> x.kco)
    deleteat!(modes, (nmodes + 1):length(modes))

    wg = RWG(; a = am, b = bm, ϵᵣ, tanδ, σ = σSpm, Rq = Rqm, modes)
    setup_modes!(wg, fHz)

    modedata = Matrix{Any}(undef, nmodes, 6)

    np2dB = 20 * log10(exp(1)) # Convert neper to dB
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
    lookup_rwg(wgspec::AbstractString) -> (a, b)

Find the dimensions of a rectangular waveguide given its "spec" (a string using standard waveguide nomenclature).

## Input Argument
- `wgspec::AbstractString`: A standard (EIA, RCSC, or IEC) abbreviation for a rectangular waveguide size.
  Examples include `"WR650"` (EIA), `"WG7"` (RCSC), and `"R22"` (IEC).

## Return Value
- `(a, b)`: Waveguide inner dimensions in the x and y directions, respectively, expressed as `Unitful` quantities
  with units of `u"inch"`.
"""
function lookup_rwg(wgspec::AbstractString)
    if wgspec == "WR2300" || wgspec == "WG0.0" || wgspec == "R3"
        a, b = 23.0u"inch", 11.5u"inch"
    elseif wgspec == "WR2100" || wgspec == "WG0" || wgspec == "R4"
        a, b = 21u"inch", 10.5u"inch"
    elseif wgspec == "WR1800" || wgspec == "WG1" || wgspec == "R5"
        a, b = 18.0u"inch", 9.0u"inch"
    elseif wgspec == "WR1500" || wgspec == "WG2" || wgspec == "R6"
        a, b = 15.0u"inch", 7.5u"inch"
    elseif wgspec == "WR1150" || wgspec == "WG3" || wgspec == "R8"
        a, b = 11.5u"inch", 5.75u"inch"
    elseif wgspec == "WR975" || wgspec == "WG4" || wgspec == "R9"
        a, b = 9.75u"inch", 4.875u"inch"
    elseif wgspec == "WR770" || wgspec == "WG5" || wgspec == "R12"
        a, b = 7.7u"inch", 3.85u"inch"
    elseif wgspec == "WR650" || wgspec == "WG6" || wgspec == "R14"
        a, b = 6.5u"inch", 3.25u"inch"
    elseif wgspec == "WR510" || wgspec == "WG7" || wgspec == "R18"
        a, b = 5.1u"inch", 2.55u"inch"
    elseif wgspec == "WR430" || wgspec == "WG8" || wgspec == "R22"
        a, b = 4.3u"inch", 2.15u"inch"
    elseif wgspec == "WG9"
        a, b = 3.5u"inch", 1.75u"inch"
    elseif wgspec == "WR340" || wgspec == "WG9A" || wgspec == "R26"
        a, b = 3.4u"inch", 1.7u"inch"
    elseif wgspec == "WR284" || wgspec == "WG10" || wgspec == "R32"
        a, b = 2.84u"inch", 1.34u"inch"
    elseif wgspec == "WG11"
        a, b = 2.372u"inch", 1.122u"inch"
    elseif wgspec == "WR229" || wgspec == "WG11A" || wgspec == "R40"
        a, b = 2.29u"inch", 1.145u"inch"
    elseif wgspec == "WR187" || wgspec == "WG12" || wgspec == "R48"
        a, b = 1.872u"inch", 0.872u"inch"
    elseif wgspec == "WR159" || wgspec == "WG13" || wgspec == "R58"
        a, b = 1.59u"inch", 0.795u"inch"
    elseif wgspec == "WR137" || wgspec == "WG14" || wgspec == "R70"
        a, b = 1.372u"inch", 0.622u"inch"
    elseif wgspec == "WR112" || wgspec == "WG15" || wgspec == "R84"
        a, b = 1.122u"inch", 0.497u"inch"
    elseif wgspec == "WR102"
        a, b = 1.02u"inch", 0.51u"inch"
    elseif wgspec == "WR90" || wgspec == "WG16" || wgspec == "R100"
        a, b = 0.9u"inch", 0.4u"inch"
    elseif wgspec == "WR75" || wgspec == "WG17" || wgspec == "R120"
        a, b = 0.75u"inch", 0.375u"inch"
    elseif wgspec == "WR62" || wgspec == "WG18" || wgspec == "R140"
        a, b = 0.622u"inch", 0.311u"inch"
    elseif wgspec == "WR51" || wgspec == "WG19" || wgspec == "R180"
        a, b = 0.51u"inch", 0.255u"inch"
    elseif wgspec == "WR42" || wgspec == "WG20" || wgspec == "R220"
        a, b = 0.42u"inch", 0.17u"inch"
    elseif wgspec == "WR34" || wgspec == "WG21" || wgspec == "R260"
        a, b = 0.34u"inch", 0.17u"inch"
    elseif wgspec == "WR28" || wgspec == "WG22" || wgspec == "R320"
        a, b = 0.28u"inch", 0.14u"inch"
    elseif wgspec == "WR22" || wgspec == "WG23" || wgspec == "R400"
        a, b = 0.224u"inch", 0.112u"inch"
    elseif wgspec == "WR19" || wgspec == "WG24" || wgspec == "R500"
        a, b = 0.188u"inch", 0.094u"inch"
    elseif wgspec == "WR15" || wgspec == "WG25" || wgspec == "R620"
        a, b = 0.148u"inch", 0.074u"inch"
    elseif wgspec == "WR12" || wgspec == "WG26" || wgspec == "R740"
        a, b = 0.122u"inch", 0.061u"inch"
    elseif wgspec == "WR10" || wgspec == "WG27" || wgspec == "R900"
        a, b = 0.1u"inch", 0.05u"inch"
    elseif wgspec == "WR8" || wgspec == "WG28" || wgspec == "R1200"
        a, b = 0.08u"inch", 0.04u"inch"
    elseif wgspec == "WR6" || wgspec == "WG29" || wgspec == "R1400"
        a, b = 0.065u"inch", 0.0325u"inch"
    elseif wgspec == "WR7" || wgspec == "WG29" || wgspec == "R1400"
        a, b = 0.065u"inch", 0.0325u"inch"
    elseif wgspec == "WR5" || wgspec == "WG30" || wgspec == "R1800"
        a, b = 0.051u"inch", 0.0255u"inch"
    elseif wgspec == "WR4" || wgspec == "WG31" || wgspec == "R2200"
        a, b = 0.043u"inch", 0.0215u"inch"
    elseif wgspec == "WR3" || wgspec == "WG32" || wgspec == "R2600"
        a, b = 0.034u"inch", 0.017u"inch"
    elseif wgspec == "WR2"
        a, b = 0.020u"inch", 0.010u"inch"
    elseif wgspec == "WR1"
        a, b = 0.010u"inch", 0.0050u"inch"
    else
        throw(ArgumentError("Unknown waveguide specification: $wgspec"))
    end
    return (a, b)
end

"""
    rwg_modetable(wgspec, frequency; kwargs...)

Pretty-print a table of rectangular waveguide mode properties: cutoff frequency, guide wavelength, and attenuation constant.

Guide wavelength and attenuation for the TE₁₀ and TE₀₁ modes are accurately computed at frequencies below, at, or 
above cutoff, for both smooth and rough imperfectly conducting surfaces, using the Gradient method as detailed 
in the paper: "Analytical waveguide model precisely predicting loss and delay including surface roughness," 
**IEEE Trans. Microwave Theory Tech.**, vol. 66, no. 6 (2018): pp. 2649-2662. For other modes the standard 
perturbational formulas are used in conjunction with the Gradient method to determine effective surface impedance.  
For cutoff modes (other than TE₁₀ and TE₀₁) the value of guide wavelength will be printed out as `Inf`.

Note: This function is intended for interactive use.  For programmatic use, see `rwg_modes`.  
The table is printed to the user's console using the `PrettyTables` package.  

## Required Positional Arguments
- `wgspec::AbstractString`: A standard (EIA, RCSC, or IEC) abbreviation for a rectangular waveguide size.
  Examples include `"WR650"` (EIA), `"WG7"` (RCSC), and `"R22"` (IEC).
- `frequency`: The frequency expressed as a `Unitful` quantity with dimensions of inverse time.  Examples: `2.35u"GHz"`,
  `2350u"MHz"`, `2.35e9u"Hz"`, etc.

## Optional Keyword Arguments
- `length_unit`:: A `Unitful:Units` object with dimension of length. Default value is `u"inch"`. The table will
  print values in the "guide wavelength" column using this unit; similarly the attenuation column will use
  units of dB/`length_unit`.
- `ϵᵣ` or `epsr`: Dielectric constant of material filling the waveguide. Default value is `1.0`.
- `tanδ` or `tandel`: Loss tangent of material filling the waveguide. Default value is `0.0`. 
- `σ` or `sigma`: The bulk conductivity of the waveguide metal walls expressed as a `Unitful` quantity
  with dimensions the same as those of `u"S/m"`.  Default value is `Inf*u"S/m"`.
- `Rq`: The RMS surface roughness of the waveguide walls expressed as a `Unitful` quantity having dimension
  of length. Default value is `0.0u"m"`.
- `col_fmts::String = ["%s", "%i", "%i", "%#.8g", "%8.4f", "%8.4f"]`: A vector of C-style format strings
  used to format the six columns of the table.
"""
function rwg_modetable(
        wgspec::AbstractString,
        f::Unitful.Quantity{<:Real, Unitful.dimension(u"Hz")};
        length_unit::Unitful.FreeUnits{Tl, Unitful.dimension(u"m")} = u"inch",
        kwargs...) where {Tl}
    (a, b) = lookup_rwg(wgspec)
    return rwg_modetable(length_unit(a), length_unit(b), f; wgname = wgspec, kwargs...)
end

"""
    rwg_modetable(a, b, frequency; kwargs...)

An alternative method for `rwg_modetable` allowing explicit specification of the waveguide dimensions.
Here, instead of `wgspec`, the first two arguments are the x and y dimensions of the waveguide, expressed
as `Unitful` quantities having dimension of length.  Both must be expressed in the same units, e.g. `u"mm"`
or `u"inch"`.  The default value of `length_unit` for this method is whatever length unit `a` and `b` use.
"""
function rwg_modetable(
        a::Unitful.Quantity{<:Real, Unitful.dimension(u"m")},
        b::Unitful.Quantity{<:Real, Unitful.dimension(u"m")},
        f::Unitful.Quantity{<:Real, Unitful.dimension(u"Hz")};
        nmodes::Int = 10,
        ϵᵣ::Real = 1.0,
        epsr::Real = 1.0,
        tanδ::Real = 0.0,
        tandel::Real = 0.0,
        σ::Unitful.Quantity{<:Real, Unitful.dimension(u"S/m")} = Inf * u"S/m",
        sigma::Unitful.Quantity{<:Real, Unitful.dimension(u"S/m")} = Inf * u"S/m",
        Rq::Unitful.Quantity{<:Real, Unitful.dimension(u"m")} = 0.0u"m",
        col_fmts = ["%s", "%i", "%i", "%#.8g", "%8.4f", "%8.4f"],
        wgname::AbstractString = "RWG")
    ϵᵣ = max(ϵᵣ, epsr)
    tanδ = max(tanδ, tandel)
    σ = min(σ, sigma)

    length_unit = Unitful.unit(a)
    freq_unit = Unitful.unit(f)

    modedata = rwg_modes(a, b, f; nmodes, ϵᵣ, tanδ, σ, Rq)

    title = wgname * ": a = $a, b = $b, f = $f"
    !isinf(σ) && (title *= ", σ = $σ")
    !iszero(Rq) && (title *= ", Rq = $Rq")
    !isone(ϵᵣ) && (title *= ", ϵᵣ = $ϵᵣ")
    !iszero(tanδ) && (title *= ", tanδ = $tanδ")
    maxtitlelen = 60
    linebreak = is_html_environment() ? "<br>\n" : "\n"
    linelen = 0
    title2 = ""
    for phrase in eachsplit(title, ',')
        if linelen + length(phrase) > maxtitlelen
            title2 = string(title2, linebreak)
            linelen = 0
        end
        linelen += length(phrase)
        title2 *= string(phrase, ",")
    end
    title = title2[begin:(end - 1)]

    column_labels = ["Type", "m", "n", "Cutoff Freq.", "Guide Wavelength", "Attenuation"]
    column_units = ["", "", "", "[$freq_unit]", "[$length_unit]", "[dB/$length_unit]"]
    header = (column_labels, column_units)

    if is_html_environment()
        p = pretty_table(
            modedata;
            backend = Val(:html),
            title,
            title_alignment = :c,
            header,
            formatters = ft_printf(col_fmts, 1:6),
            tf = tf_html_default
        )
        #Main.IJulia.display("text/html", p)
    else
        println()
        p = pretty_table(
            modedata;
            backend = Val(:text),
            title,
            title_alignment = :c,
            title_same_width_as_table = true,
            title_autowrap = true,
            header,
            formatters = ft_printf(col_fmts, 1:6),
            tf = tf_unicode_rounded
        )
    end

    return nothing
end
