#計算に関係ないパッケージ
using Pkg

include("parameter.jl")

#計算に使うパッケージ
using LinearAlgebra
using Plots

# 計算条件###################
#計算レンジ
param_range = struct_range(xwidth = 100, ywidth = 100, zwidth = 100,trange = 50)
#計算ステップ
param_step = struct_step(xstep = 0.1, ystep = 0.1, zstep = 0.05)
Nx = Int(param_range.xwidth / param_step.xstep)
Ny = Int(param_range.ywidth / param_step.ystep)
Nz = Int(param_range.zwidth / param_step.zstep)
Nt = Int(param_range.trange / param_step.tstep)
param_N = struct_N(Nx,Ny,Nz,Nt)

#材料情報
param_material = struct_materiarl(nb = 1.5, Δn0 = 0.01, τ = 0.1,α = 0)
param_mode = struct_vortex_mode(l=1, p=1)
param_beam = struct_beam_param()

#ビーム情報
param_beam = struct_beam

println("計算環境")
versioninfo()
Pkg.status()

println("計算条件")
@show param_range
@show param_step
@show param_N
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

@time function retN()
end

# 試験的。BPMのモード測定を使う
@time function showMode()
end


@time function boundary_set(u,k)
end

function initial_set(M::struct_vortex_mode,B::struct_beam,)
    
end
function main()
    #セル個数
    F0 = zeros(param_N.Nx,param_N.Ny)
    F0  = initial_set(param_mode,param_beam)
end
