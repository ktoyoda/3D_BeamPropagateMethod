#計算に関係ないパッケージ
using Pkg
include("propagator.jl")
include("Params.jl")
include("CalcFunctions.jl")
#計算に使うパッケージ
using LinearAlgebra
using Plots
using KTOptical 
# 自作のパッケージを使うときにはGithubにアップロードして、
# GithubのHTTPSをadd $$HTTPS_URL$$する。
# プライベートでよい。

using BenchmarkTools
using JLD2
using Polynomials
um = Params.um

#using FileIO
# 計算条件###################
#計算レンジ
crange = Params.crange(x = 20um, y = 20um, z = 100um, t = 2)
#計算ステップ
steps = Params.steps(x = 0.2um, y = 0.2um, z = 0.5um, t = 0.1)
Nx = Int(floor(crange.x / steps.x))
Ny = Int(floor(crange.y / steps.y))
Nz = Int(floor(crange.z / steps.z))
Nt = Int(floor(crange.t / steps.t))
N = Params.N(Nx,Ny,Nz,Nt)

#材料情報
mtr = Params.material(nb = 1.5, Δn0 = -0.04, τ = 0.1,α = 0)

#ビーム情報
# ガウスモードならgauss_mode(0,0)
# LGモードならvortex_mode(1,0)
mode = vortex_mode(1,0)
beam = Params.beam(w = 1.06um *2 , U0 = 1, wavelength = 1.064um)

println("計算環境")
#versioninfo()
Pkg.status()

println("計算条件")
@show crange
@show steps
@show N
@show mtr
@show mode
@show beam

function main()
    #セル個数
    F0 = zeros(N.x, N.y, N.z)
#    @show size(F0)
    E = initial_set(mode, beam, F0)
    gr()
    #!!!!!!!!!!!!!!!!!!!!!!!!
    x = range(-crange.x/2, crange.x/2 ,length = N.x)
    y = range(-crange.y/2, crange.y/2 ,length = N.y)
    z = range(-crange.z/2, crange.z/2 ,length = N.z)

    F_result = zeros(ComplexF64,(N.x,N.y,N.z))
    Iintegral = zeros(ComplexF64,(N.x,N.y,N.z))
    Et = zeros(ComplexF64,(N.x,N.y,N.z))
    #F_k_1stは現在のF_k
    #F_k_2ndは更新されたF_k
    @show N.x * N.y, N.x* N.y
    F_k_1st = zeros(ComplexF64, N.x, N.y) #Zerosに型指定忘れないこと！！！
    F_k_half = zeros(ComplexF64, N.x, N.y)
    F_k_2nd = zeros(ComplexF64, N.x, N.y)
    matN = fill(mtr.nb,(N.x,N.y,N.z))
    @show size(F0)
    @show size(E)
    @show typeof(E)
    @show size(F_k_1st)
    F_k_1st = E
    #setNwaveguide!(matN, steps.x, steps.y, steps.z, 0, 0, mtr.nb, mtr.nb, 0.5)

    Nref = (2mtr.nb + (mtr.nb + mtr.Δn0)) /3

    for t in 1:N.t
        @show "Zmax:", N.z
        for k in 1:N.z
            if k%10 ==0
                @show "z", k,"/",N.z,"t",t, "/", N.t
            end
            # y 固定、　x方向移動
            calcStep1!(F_k_1st, F_k_half, k, matN, Nref)
            # x 固定、　y方向移動
            calcStep2!(F_k_half, F_k_2nd, k, matN, Nref)
            F_k_1st = F_k_2nd
            F_result[:,:,k] = F_k_2nd
        end
        #!! ブロードキャストしたくない 引数があった場合はどうすればよいのだろう
        matN = matN .+ IntensityTodN!.(Iintegral, Et, F_result, mtr, t*steps.t , step)
    #    @save "/savefile/F_"* string(t)*".jld2" F_k_2nd
        Ezx = abs.(F_result[Int(floor(N.y/2)-1),: , :])
        p1 = heatmap(Ezx,levels = 200)
        display(p1)        
    end
    # @show F_result

end








@time main()

#=
Ezx = abs.(F_result[Int(floor(N.y/2)-1),: , :])
#  @show Ezx
@show typeof(Ezx)
p1 = heatmap(Ezx,levels = 200)

@show maximum(abs.(F_result[:, :, :]))
@show maximum(Ezx)

display(p1)

Ezy = abs.(F_result[: , Int(floor(N.y/2)-1),:])
p3 = heatmap(-crange.x:steps.x:crange.x, -crange.z:steps.z:crange.z,
            Ezy,levels = 200)#,clim = (0,0))#,lim=(0.25,0))
@show size(F_result)
display(p3)

p2 = plot(Ezx[Int(floor(N.y/2)),:])
display(p2)
=#