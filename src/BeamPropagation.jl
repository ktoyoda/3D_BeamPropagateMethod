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
crange = Params.crange(x = 200um, y = 200um, z = 800um, t = 1)
#計算ステップ
#steps = Params.steps(x =1um, y = 1um, z = 0.5um, t = 1)
steps = Params.steps(x =5um, y = 5um, z = 5um, t = 1)

#steps = Params.steps(x = 1um, y = 1um, z = 1um, t = 0.1)
Nx = Int(floor(crange.x / steps.x))
Ny = Int(floor(crange.y / steps.y))
Nz = Int(floor(crange.z / steps.z))
Nt = Int(floor(crange.t / steps.t))
N = Params.N(Nx,Ny,Nz,Nt)
# 初期で作られる屈折率の強度依存分布と、伝搬でできる屈折率分布の比
ratio = 0.001
# smoothing range
smooth_range = 1

#材料情報
#Δn0は0.01以下にとらないと発散する。これを防ぐにはパでを使うしかない
mtr = Params.material(nb = 1.5, Δn0 = -0.01, τ = 0, α = 0, U = 1)

#ビーム情報
# ガウスモードならgauss_mode(0,0)
# LGモードならvortex_mode(1,0)
mode = vortex_mode(1,0)
beam = Params.beam(w = 25um , U0 = 1, wavelength = 0.532um)

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
@show ratio

@time function main()


    #セル個数
    F0 = zeros(N.x, N.y, N.z)
#    @show size(F0)
    E = initial_set(mode, beam, F0)
    gr()
    #!!!!!!!!!!!!!!!!!!!!!!!!
    x = range(-crange.x/2, crange.x/2 ,length = N.x)
    y = range(-crange.y/2, crange.y/2 ,length = N.y)
    z = range(0, crange.z ,length = N.z)
    @show x,length(x)
    F_result = zeros(ComplexF64,(N.x,N.y,N.z))
    Iintegral = zeros(Float64,(N.x,N.y,N.z))
    Et = zeros(Float64,(N.x,N.y,N.z))
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
    #ratioを大きくしすぎるとエラー！
    setNwaveguide!(matN, mode, ratio, steps.x, steps.y, steps.z, mtr.nb, mtr.Δn0, beam)
    matNzx = matN[Int(floor(N.y/2)-1),: , :]
    if mtr.Δn0 >0
        cl = (mtr.nb, mtr.nb+mtr.Δn0)
    else
        cl = (mtr.nb + mtr.Δn0 ,mtr.nb)
    end
    p = heatmap(z*10^3, x*10^6, matNzx, levels = 200, clim=cl)
    display(p)
    Nref = mtr.nb
    
        
    for t in 1:N.t
        @show "Zmax:", N.z
        F_k_1st = E
        F_k_half = zeros(ComplexF64, N.x, N.y)
        F_k_2nd = zeros(ComplexF64, N.x, N.y)
        for k in 1:N.z
            if k%100 ==0
                @show "z", k,"/",N.z,"t",t, "/", N.t
                @show findmax(abs2.(F_result))
                @show abs2.(F_result[N.x÷2,N.y÷2,N.z÷2])
                @show minimum(matN)
            end
            # y 固定、　x方向移動
            calcStep1!(F_k_1st, F_k_half, k, matN, Nref)
            # x 固定、　y方向移動
            calcStep2!(F_k_half, F_k_2nd, k, matN, Nref)
            F_k_1st = F_k_2nd
            F_result[:,:,k] = F_k_2nd
        end

        #!! ブロードキャストしたくない 引数があった場合はどうすればよいのだろう
        #!!!!→ スカラとして渡したい引数にRefを付ける！
        matN_renew = matN .+ IntensityTodN.(Iintegral, Et, F_result, Ref(mtr), Ref(t*steps.t) , Ref(steps),Ref(ratio))
    #    @show IntensityTodN.(Iintegral, Et, F_result, Ref(mtr), Ref(t*steps.t) , Ref(steps),Ref(ratio))
        matN = Smoothing(matN_renew,x,y,z,5)
        # 常にインプットされるのでEをF_k_1stにいれる

    #    @save "/savefile/F_"* string(t)*".jld2" F_k_2nd
        Ezx2 = view(F_result,:, Int(floor(N.y/2)-1), :)
        Ezx = abs.(Ezx2)
        @show size(Ezx)
        matNzx = view(matN,:, Int(floor(N.y/2)-1), :)
        p1 = heatmap(z*10^3, x*10^6, Ezx,levels = 200,clim=(0,3),xlabel = "Z [mm]", ylabel = "X [μm]" )
        strt = string(t)
        savefig(p1, strt*"_In.png")
        display(p1)
        if mtr.Δn0 >0
            cl = (mtr.nb, mtr.nb+mtr.Δn0)
        else
            cl = (mtr.nb + mtr.Δn0 ,mtr.nb)
        end
        p3 = heatmap(z*10^3, x*10^6, matNzx, levels = 200, clim=(1.40,1.5),xlabel = "Z [mm]", ylabel = "X [μm]")
        savefig(p3, strt*"_ri.png")
        display(p3)
        jldopen("xyz1105um_1s_"*strt*"_E.jld2","a+") do file1
            file1[strt*"/E"] = abs.(F_result)
        end
        jldopen("xyz1105um_1s_"*strt*"_N.jld2","a+") do file2
            file2[strt*"/N"] = matN
        end
    end
    return F_result
end

@time F_result = main()

