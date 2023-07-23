using Documenter, HydroTools

CI = get(ENV, "CI", nothing) == "true"

# Logging.disable_logging(Logging.Warn)

# Make the docs, without running the tests again
# We need to explicitly add all the extensions here
makedocs(
  modules=[
    HydroTools
  ],
  format=Documenter.HTML(
    prettyurls=CI,
  ),
  pages=[
    "Introduction" => "index.md",
    "Meteorology" => "Meteorology.md",
    "Hydrology" => "Hydrology.md",
    "Potential Evapotranspiration models" => "PET.md",
    "Extreme Climate indexes" => "ExtremeClimate.md"
  ],
  sitename="HydroTools.jl",
  strict=false,
  clean=false,
)

# Enable logging to console again
# Logging.disable_logging(Logging.BelowMinLevel)

deploydocs(
  repo="github.com/CUG-hydro/HydroTools.jl.git", 
)
