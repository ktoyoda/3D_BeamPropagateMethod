include("parameter.jl")

# 計算条件###################
#計算レンジ
param_range = struct_range(xwidth = 100, ywidth = 100, zwidth = 100,trange = 50)
#計算ステップ
param_step = struct_step(( xstep = 0.1, ystep = 0.1, zstep = 0.05)
#材料情報
param_material = struct_materiarl((nb = 1.5, Δn0 = 0.01, τ = 0.1,α = 0)