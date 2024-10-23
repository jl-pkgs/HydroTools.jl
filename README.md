# HydroTools in Julia

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://jl-pkgs.github.io/HydroTools.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://jl-pkgs.github.io/HydroTools.jl/dev)
[![CI](https://github.com/jl-pkgs/HydroTools.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/jl-pkgs/HydroTools.jl/actions/workflows/CI.yml)
[![Codecov](https://codecov.io/gh/jl-pkgs/HydroTools.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/jl-pkgs/HydroTools.jl/tree/master)

> Translation of [hydroTools](https://github.com/CUG-hydro/hydroTools) in R

## Installation

```julia
] add https://github.com/jl-pkgs/HydroTools.jl
```

## Contents

- [x] 水汽转换函数：VPD, ea, es, q, RH
- [x] 极端降水指标、热浪指标、体感温度
- [x] 潜在蒸散发
- [x] 冠层辐射传输（n层）
- [x] 土壤水运动
- [x] 土壤热通量
- [x] 蒸散发与土壤热通量联合求解

- [x] 优化函数sceua

## TODO

温度求解存在较大误差，不确定误差来源`辐射`还是`参数`。

- [ ] 尝试采用ODE的方法解土壤温度

- [ ] 已知温度，反推`κ`, `cv`

## References

1. <https://github.com/CUG-hydro/hydroTools>
