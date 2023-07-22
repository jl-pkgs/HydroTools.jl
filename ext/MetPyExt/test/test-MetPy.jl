using HydroTools

using PyCall

# end

# import metpy.calc as mpcalc
# from metpy.units import units

@time ca = pyimport("metpy.calc")
units = pyimport("metpy.units").units

Ta = [10, 20, 30.] * units.degC
RH = [80, 70, 75.] * units.percent

ca.heat_index(Ta, RH)


heat_index.(Ta.magnitude, RH.magnitude)
