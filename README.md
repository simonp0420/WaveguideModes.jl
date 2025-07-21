[![CI](https://github.com/simonp0420/WaveguideModes.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/simonp0420/WaveguideModes.jl/actions/workflows/CI.yml)
# WaveguideModes

Contains useful data structures and utilty programs for microwave engineers who work with metallic waveguides.  This 
first release is limited to uniformly filled rectangular waveguides.  Next up will be circular guides, possibly
followed by ridged guides and partially filled guides.

## Rectangular Waveguides
The most useful function to call from the Julia REPL is `rwg_modetable`, which prints out a nicely formatted table
of cutoff frequencies, guide wavelengths, and attenuation constants for rectangular waveguides.  The waveguide size
can be specified either by explicitly passing the `a` and `b` dimensions (using any 
[`Unitful`](https://painterqubits.github.io/Unitful.jl/stable/) length quantity), or by passing a string containing
a "waveguide spec" (a standard (EIA, RCSC, or IEC) abbreviation for a rectangular waveguide size). Examples include 
`"WR650"` (EIA), `"WG7"` (RCSC), and `"R22"` (IEC). 

Here is a short video showing the use of the `rwg_modetable` function:

https://github.com/user-attachments/assets/f59602ac-3497-4d69-953c-eaff624073a7

Here is the full documentation for `rwg_modetable`:

    rwg_modetable(wgspec, frequency; kwargs...)

Pretty-print a table of rectangular waveguide mode properties: cutoff frequency, guide wavelength, and attenuation constant.

Guide wavelength and attenuation for the TE₁₀ and TE₀₁ modes are accurately computed at frequencies below, at, or 
above cutoff, for both smooth and rough imperfectly conducting surfaces, using the Gradient method as detailed 
in the paper: "Analytical waveguide model precisely predicting loss and delay including surface roughness," 
**IEEE Trans. Microwave Theory Tech.**, vol. 66, no. 6 (2018): pp. 2649-2662. For other modes the standard 
perturbational formulas are used in conjunction with the Gradient method to determine effective surface impedance.  
For cutoff modes (other than TE₁₀ and TE₀₁) the value of guide wavelength will be printed out as `Inf`.

Note:This function is intended for interactive use.  For programmatic use, see `rwg_modes`.  
The table is printed to the user's console using the `PrettyTables` package.  

## Required Positional Arguments
- `wgspec::AbstractString`: A standard (EIA, RCSC, or IEC) abbreviation for a rectangular waveguide size.
  Examples include `"WR650"` (EIA), `"WG7"` (RCSC), and `"R22"` (IEC).
- `frequency`: The frequency expressed as a [`Unitful`](https://painterqubits.github.io/Unitful.jl/stable/)
  quantity with dimensions of inverse time.  Examples: `2.35u"GHz"`, `2350u"MHz"`, `2.35e9u"Hz"`, etc.

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

    rwg_modetable(a, b, frequency; kwargs...)

An alternative method for `rwg_modetable` allowing explicit specification of the waveguide dimensions.
Here, instead of `wgspec`, the first two arguments are the x and y dimensions of the waveguide, expressed
as [`Unitful`](https://painterqubits.github.io/Unitful.jl/stable/) quantities having dimension of length.
Both must be expressed in the same units, e.g. `u"mm"` or `u"inch"`.  The default value of `length_unit`
for this method is whatever length unit `a` and `b` use.
