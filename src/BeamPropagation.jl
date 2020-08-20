#計算に関係ないパッケージ
using Pkg

include("parameter.jl")

#計算に使うパッケージ
using LinearAlgebra
using SparseArrays
using Plots

# 計算条件###################
#計算レンジ
param_range = struct_range(xwidth = 100, ywidth = 100, zwidth = 100,trange = 50)
#計算ステップ
param_step = struct_step(xstep = 0.1, ystep = 0.1, zstep = 0.05)
#材料情報
param_material = struct_materiarl(nb = 1.5, Δn0 = 0.01, τ = 0.1,α = 0)
param_mode = struct_vortex_mode()

println("計算環境")
versioninfo()
Pkg.status()

println("計算条件")
@show param_range
@show param_step
@show param_material
@show param_mode

@time function returnPades(T::Int16,N::Int16)
end

#ADIのX方向差分
@time function calcStep1()
end

#ADIのY方向差分
@time function calcStep2()
end

@time function showMode()
end

@time function retN()
end

@time function boundary_set(u,k)
  
end


