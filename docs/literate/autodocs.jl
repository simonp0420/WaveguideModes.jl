# # Public API Reference

# ## Publicly Exported Functions

# ```@docs
# lookup_rwg
# rwg_modes
# rwg_modetable(::AbstractString, ::Unitful.Quantity{<:Real, Unitful.dimension(u"Hz")})
# rwg_modetable(::Unitful.Quantity{<:Real, Unitful.dimension(u"m")},
#               ::Unitful.Quantity{<:Real, Unitful.dimension(u"m")},
#               ::Unitful.Quantity{<:Real, Unitful.dimension(u"Hz")})
# rwgte10gz
# ```

# ## Publicly Exported Types

# ```@docs
# RWG
# RWGMode
# ```

# # Non-Exported API (Subject to Change)

# ## Private Functions

# ```@autodocs
# Modules = [WaveguideModes]
# Order   = [:function]
# Private = true
# Public = false
# ```

# ## Private Types

# ```@autodocs
# Modules = [WaveguideModes]
# Order   = [:type]
# Private = true
# Public = false
# ```

