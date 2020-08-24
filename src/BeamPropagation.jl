#計算に関係ないパッケージ
using Pkg

include("parameter.jl")

#計算に使うパッケージ
using LinearAlgebra
using Plots
using KTOptical
using .Param
# 計算条件###################
#計算レンジ
range = Param.range(xwidth = 100um, ywidth = 100um, zwidth = 100um,trange = 10)
#計算ステップ
step = Param.step(xstep = 0.1um, ystep = 0.1um, zstep = 0.05um)
Nx = Int(param_range.xwidth / param_step.xstep)
Ny = Int(param_range.ywidth / param_step.ystep)
Nz = Int(param_range.zwidth / param_step.zstep)
Nt = Int(param_range.trange / param_step.tstep)
N = Param.N(Nx,Ny,Nz,Nt)

#材料情報
material = Param.materiarl(nb = 1.5, Δn0 = 0.01, τ = 0.1,α = 0)
mode = Param.vortex_mode(l=1, p=0)

#ビーム情報
beam = Param.beam(w = 20um, U0 = 100, wavelength = 1.06um)

println("計算環境")
versioninfo()
Pkg.status()

println("計算条件")
@show param_range
@show param_step
@show param_N
@show param_material
@show param_mode
@show param_beam

println("//////////////////////////////")
function returnPades(T::Int16,N::Int16)
end

#ADIのX方向差分
function calcStep1(F_k11, F_kp12)
    for i in range(1,length = Nx)
        for j in range(1,length = Ny)
            continue
        end
    end
    boundary_set(u,k,F_kp12)
    return F_kp12
end

#ADIのY方向差分
function calcStep1(F_k12, F_kp21)
    for i in range(1,length = Nx)
        for j in range(1,length = Ny)
            continue
        end
    end
    boundary_set(u,k,F_kp21)
    return F_kp21
end

function retN()
end

# 試験的。BPMのモード測定を使う
function showMode()
end


function boundary_set(u,k,array)
end

@time function initial_set(M, B ,F0)
    KTOptical.setParam(B.w, 0, B.wavelength)
    x = range(-param_range.xwidth/2, param_range.xwidth/2 ,step = param_step.xstep)
    y = range(-param_range.ywidth/2, param_range.ywidth/2 ,step = param_step.ystep)

    E = LG_I.(M.l, M.p,x,y')
    return E

end

function main()
    #セル個数
    F0 = zeros(param_N.Nx, param_N.Ny, param_N.Nz)
    E = initial_set(param_mode, param_beam, F0)
    gr()
    x = range(-param_range.xwidth/2, param_range.xwidth/2 ,step = param_step.xstep)

    F_result = zeros(Float64,(param_N.Nx,param_N.Ny,param_N.Nz));
    #F_k_1stは現在のF_k
    #F_k_2ndは更新されたF_k
    F_k_1st = zeros(param_N.Nx,param_N.Ny);
    F_k_2nd = zeros(param_N.Nx,param_N.Ny);
#    #左辺は更新されたF_k_1stが入る。
    #F_k_2ndは関数内部で毎回宣言した方がいいのか、
    #それとも一度宣言して引数として与えた方がいいのか（いまはこれ
    #は要検証。
    F_k_1st = calcStep1(F_k_1st,F_k_2nd);
    F_k_1st = calcStep1(F_k_1st,F_k_2nd);
    #push(F_k_1st )
end

@time main()