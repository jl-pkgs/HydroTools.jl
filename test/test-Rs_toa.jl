using Test
using HydroTools
using Dates

@testset "Rs_toa Comparison" begin
  # 参数设置
  lat_deg = 30.0
  doy = 150
  hour_val = 15.0 # 下午3点
  
  # 构造时间
  # 注意：HydroTools 内部计算时通常不依赖年份的具体值，只依赖 DOY
  base_date = DateTime(2010) 
  time_inst = base_date + Day(doy) + Hour(hour_val)
  
  # 1. 瞬时值测试 (Instantaneous)
  # 基准值来自标准太阳辐射几何公式计算 (Sc * dr * coszen)
  # Benchmark: ~1001.97 W/m2
  expected_inst = 1001.97 
  ht_inst = cal_Rsi_toa_inst([time_inst]; lat=lat_deg)[1]
  
  # 比较结果 (允许 0.1% 的误差，涵盖公式系数的微小差异)
  @test isapprox(ht_inst, expected_inst, rtol=1e-3)
  
  # 2. 小时平均值测试 (Hourly Average)
  # 基准值来自该小时内的积分平均
  # Benchmark: ~894.94 W/m2
  expected_hourly = 894.94
  ht_hourly = cal_Rsi_toa_hour([time_inst]; lat=lat_deg)[1]
  
  # 比较结果
  @test isapprox(ht_hourly, expected_hourly, rtol=1e-3)
end
