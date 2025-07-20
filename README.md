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


https://github.com/user-attachments/assets/e4fe63cc-c65f-4063-a304-27fe43542c95

