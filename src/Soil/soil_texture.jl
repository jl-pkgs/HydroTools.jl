# TODO: 这里使用枚举型更合适

module USDA

# https://developers.google.com/earth-engine/datasets/catalog/OpenLandMap_SOL_SOL_TEXTURE-CLASS_USDA-TT_M_v02
@enum SoilTexture begin
  # USDA_UNDEFINED  # NA  
  CLAY            = 1 # CL
  SILTY_CLAY      = 2 # SICL
  SANDY_CLAY      = 3 # SACL
  CLAY_LOAM       = 4 # CLLO
  SILTY_CLAY_LOAM = 5 # SICLLO
  SANDY_CLAY_LOAM = 6 # SACLLO
  LOAM            = 7 # LO
  SILTY_LOAM      = 8 # SILO
  SANDY_LOAM      = 9 # SALO
  SILT            = 10 # SI
  LOAMY_SAND      = 11 # LOSA
  SAND            = 12 # SA
end

Base.to_index(i::SoilTexture) = Int(i)

# https://github.com/NigelVanNieuwenhuizen/USDA-Soil-Texture-Calculator/blob/9cbf7ee58d384074de26e9a50fa4d1872e9a37f8/USDASoilTextureCalculator.py#L261C17-L286C46

"""
    soil_texture(sand::Float64, silt::Float64)

# Arguments
- `sand`: percentage of sand, %
- `silt`: percentage of silt, %

# Example usage:
```
sand = 60.0
silt = 20.0
println("Soil Texture: ", soil_texture(sand, silt))
```
"""
function soil_texture(sand::Float64, silt::Float64)
  clay = 100 - sand - silt

  if sand <= 45 && silt <= 40 && clay >= 40
    CLAY
  elseif sand <= 65 && sand >= 45 && silt <= 20 && clay >= 35 && clay <= 55
    SANDY_CLAY
  elseif sand <= 20 && silt >= 40 && silt <= 60 && clay >= 40 && clay <= 60
    SILTY_CLAY
  elseif sand >= 45 && sand <= 80 && silt <= 28 && clay >= 20 && clay <= 35
    SANDY_CLAY_LOAM
  elseif sand >= 20 && sand <= 45 && silt >= 15 && silt <= 53 && clay >= 27 && clay <= 40
    CLAY_LOAM
  elseif sand <= 20 && silt >= 40 && silt <= 73 && clay >= 27 && clay <= 40
    SILTY_CLAY_LOAM
  elseif sand >= 43 && sand <= 85 && silt <= 50 && clay <= 20
    SANDY_LOAM
  elseif sand >= 23 && sand <= 52 && silt >= 28 && silt <= 50 && clay >= 7 && clay <= 27
    LOAM
  elseif sand <= 50 && silt >= 50 && silt <= 88 && clay <= 27
    SILT_LOAM
  elseif sand <= 20 && silt >= 80 && clay <= 12
    SILT
  elseif sand >= 70 && sand <= 90 && silt <= 30 && clay <= 15
    LOAMY_SAND
  elseif sand >= 85 && silt <= 15 && clay <= 10
    SAND
  else
    @error("Soil texture not found for sand: $sand, silt: $silt, clay: $clay")
  end
end

export soil_texture

end
