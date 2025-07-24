#=
# Introduction
[`WaveguideModes`](https://github.com/simonp0420/WaveguideModes.jl) is a Julia package intended to 
assist microwave engineers and others who work with metallic waveguides.  

## For Microwave Engineers/Designers
Convenient, user-oriented functions `rwg_modetable` and `cwg_modetable` (the latter coming "soon"), 
are called from the Julia REPL to
produce pretty-printed tables of modal cutoff frequency, guide wavelength, and attenuation constant for guides whose 
walls may be lossy and non-smooth. Here's an example:
=#
using WaveguideModes
rwg_modetable("WR62", 12u"GHz"; σ = 58u"MS/m", Rq = 20e-6u"inch", length_unit=u"cm")
#=
In the case of the rectangular waveguide TE₁₀ and TE₀₁ modes, highly accurate 
formulas from [lomakin2018](@cite) are used. For the other rectangular waveguide modes and all circular waveguide 
modes, standard perturbational formulas [collin2007foundations](@cite) are used to incorporate effects of metal wall
losses into above-cutoff attenuation constants, using an effective surface resistance from [grujic2022](@cite) to 
accurately account for finite surface roughness. See the Tutorial section of the documentation for more details.


## For Julia Programmers
`WaveguideModes` implements the following hierarchies of types:

```mermaid
flowchart LR
    Any --> Waveguide
    Waveguide --> MetallicWaveguide
    MetallicWaveguide --> HomogeneousMetallicWaveguide
    HomogeneousMetallicWaveguide --> RWG
    HomogeneousMetallicWaveguide --> CWG
```

```mermaid
flowchart LR
    Any --> Mode
    Mode --> MetalPipeMode
    MetalPipeMode --> HomogeneousMetalPipeMode
    HomogeneousMetalPipeMode --> RWGMode
    HomogeneousMetalPipeMode --> CWGMode
```

`RWG` and `RWGMode` are concrete types for representing rectangular homogeneously-filled, metallic waveguides
and the modes they support.  `CWG` and `CWGMode` are the corresponding types for circular guides.
These types along with their constructors and utility functions such as `setup_modes!` provide a convenient software
"scaffolding" for those wishing to implement algorithms dealing with waveguide modes, e.g. mode matching analysis.

# Package Installation
You can obtain WaveguideModes using Julia's Pkg REPL-mode (hitting `]` as the first character of the command prompt):

```julia
(v1.10) pkg> add WaveguideModes
```

or with `using Pkg; Pkg.add("WaveguideModes")`.
=#

