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
param_beam = struct_beam()

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
@time function calcStep1(F_k11, F_kp12)

    for i in range(Nx)
        for j in range(Ny)
<<<<<<< HEAD
            continue
        end
    end
    boundary_set(u,k,F_kp12)
=======
            
        end
    end
    boundary_set(F_kp12)
>>>>>>> 553b1f734a59bd075a3c44afc74b4fec5b16b0d5
    return F_kp12
end

#ADIのY方向差分
@time function calcStep2(F_k12, F_kp21)
    
    for i in range(Nx)
        for j in range(Ny)
<<<<<<< HEAD
            continue
        end
    end
    boundary_set(u,k,F_kp21)
=======
            
        end
    end
    boundary_set(F_kp21)
>>>>>>> 553b1f734a59bd075a3c44afc74b4fec5b16b0d5
    return F_kp21

end

@time function retN()
end

# 試験的。BPMのモード測定を使う
@time function showMode()
end


@time function boundary_set(u,k,array)
end

function initial_set(M::struct_vortex_mode,B::struct_beam,)
    
end
function main()
    #セル個数
    F0 = zeros(param_N.Nx, param_N.Ny, param_N.Nz)
    F0[:,:,0]  = initial_set(param_mode,param_beam)
    F_result = zeros(param_N.Nx,param_N.Ny,param_N.Nz,param_N.Nt)
    
    #F_k_1stは現在のF_k
    #F_k_2ndは更新されたF_k
    F_k_1st = zeros(param_N.Nx,param_N.Ny)
    F_k_2nd = zeros(param_N.Nx,param_N.Ny)

    #左辺は更新されたF_k_1stが入る。
    #F_k_2ndは関数内部で毎回宣言した方がいいのか、
    #それとも一度宣言して引数として与えた方がいいのか（いまはこれ
    #は要検証。
    F_k_1st = calcStep1(F_k_1st,F_k_2nd)
    F_k_1st = calcStep1(F_k_1st,F_k_2nd)
    #push(F_k_1st )

<<<<<<< HEAD
=======

>>>>>>> 553b1f734a59bd075a3c44afc74b4fec5b16b0d5
end
