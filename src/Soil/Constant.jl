export ρ_wat, ρ_ice, SILT, SAND, CLAY, θ_S


ρ_wat = 1000.0                       # Density of water (kg/m3)
ρ_ice = 917.0                        # Density of ice (kg/m3)

λ_fus = 0.3337e6                     # Heat of fusion for water at 0 C (J/kg)
tfrz  = K0                           # freezing Temperature (K)

"""
# Initialize soil texture variables

## Soil texture classes (Cosby et al. 1984)
%  1: sand
%  2: loamy sand
%  3: sandy loam
%  4: silty loam
%  5: loam
%  6: sandy clay loam
%  7  silty clay loam
%  8: clay loam
%  9: sandy clay
% 10: silty clay
% 11: clay

## References
- Cosby et al. 1984. Water Resources Research 20:682-690, Soil texture classes
- Clapp and Hornberger. 1978. Water Resources Research 14:601-604
"""
#  (Cosby et al. 1984. Water Resources Research 20:682-690)
SILT = [5.0, 12.0, 32.0, 70.0, 39.0, 15.0, 56.0, 34.0, 6.0, 47.0, 20.0] # Percent silt
SAND = [92.0, 82.0, 58.0, 17.0, 43.0, 58.0, 10.0, 32.0, 52.0, 6.0, 22.0] # Percent sand
CLAY = [3.0, 6.0, 10.0, 13.0, 18.0, 27.0, 34.0, 34.0, 42.0, 47.0, 58.0] # Percent clay

# Volumetric soil water content (%) at saturation (porosity)
# (Clapp and Hornberger. 1978. Water Resources Research 14:601-604)
θ_S = [0.395, 0.410, 0.435, 0.485, 0.451, 0.420, 0.477, 0.476, 0.426, 0.492, 0.482]


# 饱和质量
# cal_MS
