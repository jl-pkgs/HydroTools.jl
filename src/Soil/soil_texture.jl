# TODO: 这里使用枚举型更合适

# https://github.com/NigelVanNieuwenhuizen/USDA-Soil-Texture-Calculator/blob/9cbf7ee58d384074de26e9a50fa4d1872e9a37f8/USDASoilTextureCalculator.py#L261C17-L286C46

function soil_texture(sand::Float64, silt::Float64, clay::Float64)
  if sand <= 45 && silt <= 40 && clay >= 40
    "Clay"
  elseif sand <= 65 && sand >= 45 && silt <= 20 && clay >= 35 && clay <= 55
    "Sandy Clay"
  elseif sand <= 20 && silt >= 40 && silt <= 60 && clay >= 40 && clay <= 60
    "Silty Clay"
  elseif sand >= 45 && sand <= 80 && silt <= 28 && clay >= 20 && clay <= 35
    "Sandy Clay Loam"
  elseif sand >= 20 && sand <= 45 && silt >= 15 && silt <= 53 && clay >= 27 && clay <= 40
    "Clay Loam"
  elseif sand <= 20 && silt >= 40 && silt <= 73 && clay >= 27 && clay <= 40
    "Silty Clay Loam"
  elseif sand >= 43 && sand <= 85 && silt <= 50 && clay <= 20
    "Sandy Loam"
  elseif sand >= 23 && sand <= 52 && silt >= 28 && silt <= 50 && clay >= 7 && clay <= 27
    "Loam"
  elseif sand <= 50 && silt >= 50 && silt <= 88 && clay <= 27
    "Silt Loam"
  elseif sand <= 20 && silt >= 80 && clay <= 12
    "Silt"
  elseif sand >= 70 && sand <= 90 && silt <= 30 && clay <= 15
    "Loamy Sand"
  elseif sand >= 85 && silt <= 15 && clay <= 10
    "Sand"
  else
    "Not Available"
  end
end

# # Example usage:
# sand_percentage = 60.0
# silt_percentage = 20.0
# clay_percentage = 20.0

# println("Soil Texture: ", soil_texture(sand_percentage, silt_percentage, clay_percentage))
