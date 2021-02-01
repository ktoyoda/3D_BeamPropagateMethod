#計算に関係ないパッケージ
using Pkg
include("propagator.jl")
include("Params.jl")
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
crange = Params.crange(x = 20um, y = 20um, z = 10um, t = 0.1)
#計算ステップ
steps = Params.steps(x = 1um, y = 1um, z = 1um, t = 0.1)
Nx = Int(crange.x / steps.x)
Ny = Int(crange.y / steps.y)
Nz = Int(crange.z / steps.z)
Nt = Int(crange.t / steps.t)
N = Params.N(Nx,Ny,Nz,Nt)

#材料情報
material = Params.materiarl(nb = 1.5, Δn0 = 0.2, τ = 0.1,α = 0)
mode = Params.vortex_mode(l=1, p=0)

#ビーム情報
beam = Params.beam(w = 2um, U0 = 100.0, wavelength = 1.064um)

println("計算環境")
#versioninfo()
Pkg.status()

println("計算条件")
@show crange
@show steps
@show N
@show material
@show mode
@show beam

#println("//////////////////////////////")
function returnPades(T::Integer,N::Integer)
    #Julia に型注釈は必要ないが、Padeの引数は整数型に限定したい。という説明

end

# 
# 2020/8/31
# 各Zに対して、calcstepを１，２の順で回していく。
# step 1 は、
# ADIの未知数Y方向(j) 定数X方向(i) 差分
# 
# Ax = B をつくってトーマスアルゴリズムJFで解く
# juliaの場合はB\Aでよい。(バックスラッシュ！）
# Aが N.ｙ＊N.ｙ(z=k+1の係数)
# すなわち、(x,y,z) = (i,j,k)において、
# i、k 固定で 
# A[J,J] の時、j=J
# A[J±1, J]の時、j=J±1
# である。 
# xがN.yサイズのベクトル
# これは各座標の波動方程式の解であるため、考慮する必要はない。
#
# BがN.yサイズのベクトル(z=kの係数)
# 
###### 各A,Bに入る屈折率項は座標情報が必要である。
# Aに入るのは()
# 
#F_k1 既知のビーム伝搬
#F_k_after 未知のビーム伝搬の格納先(F_k+1/2)
#k 今のZカウント
# calcStep1は
# 「未知数M個(Y方向)、バンド幅３の連立一次方程式をN個(X方向)個解く」
# 藪哲郎光導波路132Pより
#
# つまり、calcStep1の内部ではN.x回ループを回し、
# その都度Ax=Bを作成して解く。
# z = kに関しては
# Step1,Step2の組み合わせを複数回周回することになる。
# t = n に関しては、
# 上記で作った新しい屈折率マップをもとに、
# 再度回していく。 	\1/2 



function calcStep1!(F_k_before, F_k_after, k, matN)
    k0 = 2π / beam.wavelength
    Nref = 1.35
    # a,b,cは左辺用
    a = -1/(steps.y)^2
    b(i,j,k) = 1im*4*k0*matN[i, j, k]/steps.z + 2/steps.y^2 - k0^2(matN[i, j, k]^2 - Nref^2)/2
    c = 1/(steps.x)^2
    
    # d は 右辺用
    d(i,j,k) = 1im*4*k0*matN[i, j, k]/steps.z - 2/steps.x^2 + k0^2(matN[i, j, k]^2 - Nref^2)/2
    B = zeros(ComplexF64,N.y,1)

    # 各回 iを固定して、jxj　の行列を作って計算する。
    # よって、最初のループは固定方向のx(i)で、
    # 内部で、diagmで行列方向のN.y x N.y = jxjを作る。
    for i in 1:N.x
        # Ax=BのA (z = k+1 における係数行列)を作る。
        # mapでベクトルを作っておいて、diagmで対角行列にする。
        # bが位置によって値が異なるのでmapで対応。
        # diagmは配列を正方行列の対角に配置する。
        # pair（=>) でズレを表現できる。
        A = diagm(map(j -> b(i,j,k),1:N.y))
        A += diagm(1 => fill(c,N.y-1))
        A += diagm(-1 => fill(a,N.y-1))

        #透明境界条件(TBC) for A#######
        # 参考文献  左貝潤一 光導波路の電磁界解析 p.145
        #左端---------------------
        imKxL = -(1/steps.y)   * log(abs(F_k_before[i,2]/F_k_before[i,1]))
        reKxL = -(1im/steps.y) * log(abs(F_k_before[i,2]/F_k_before[i,1]*exp(imKxL*steps.y)))

        # reKxLおよびreKxRが負の時、
        # 左貝だと符号逆転させて、藪だと0にすると書いてある。前者のほうがそれっぽい気がする。
        if real(reKxL)<0
            reKxL = -1*reKxL      #左貝方式
            # reKxL = 0           #藪方式
        end
        ###!!!ηLの計算違う？左貝144p。どっちが正解？１
        #ηL = exp(1im* reKxL * steps.y - imKxL*steps.y)
        ηL = exp(1im* (reKxL+ imKxL) *steps.y)

        ###!!! ここA[1,1]か？ ← 1,1でいい
        ###!!! 左貝(7.25b)と(7,37a)比較して 差分をとる。
        A[1,1] -= ηL/(steps.y)^2

        #右端---------------------
        imKxR = (1/steps.y) * log(abs(F_k_before[i,N.y]/F_k_before[i,N.y-1]))
        reKxR = (1im/steps.y) * log(abs(F_k_before[i,N.y-1]/F_k_before[i,N.y]*exp(-imKxR*steps.y)))

        if real(reKxR)<0
            reKxR = -reKxR      #左貝方式
            # reKxR = 0           #藪方式
        end
        ηR = exp(1im* (reKxR + imKxR) * steps.y)
        A[N.y,N.y] -= ηR/(steps.y)^2

        #########################

        
        # 透明境界条件(TBC) for B#######
        colBL = (2-ηL)/steps.y^2 - (matN[i,1,k]^2-Nref^2)*k0^2 + (4im*Nref*k0)/steps.z
        colBR = (2-ηR)/steps.y^2 - (matN[i,N.y,k]^2-Nref^2)*k0^2 + (4im*Nref*k0)/steps.z
        colC = -1/(steps.x)^2

        #@show F_k_before

        for j in 1:N.y
            if i == 1
                B[j] = -conj(colBL)*F_k_before[1,j] - colC*F_k_before[2,j]
            elseif i == N.x
                B[j] = -conj(colBR)*F_k_before[N.x,j] - colC*F_k_before[N.x-1,j]
            else
                B[j] = c*F_k_before[i,j] + d(i,j,k)*F_k_before[i-1,j] + d(i,j,k)* F_k_before[i+1,j]
            end
        end

        
        #########################

        F_k_after[i,:] = B\A
        #F_k_after[]を使って、新しい屈折率マップを作る。

    end
    return F_k_after
end

#ADIの未知数X方向 定数Y方向 差分
function calcStep2!(F_k_before, F_k_after,k, matN)
    #正確には、k+1とKの屈折率の平均を取るべきだと思うが、
    #今のところはk+1を抜き出す   
    k0 = 2π / beam.wavelength
    Nref = 1.35
    # a,b,cは左辺用
    a = -1/(steps.x)^2
    b(i,j,k) = 1im*4*k0*matN[i, j, k]/steps.z + 2/steps.x^2 - k0^2(matN[i, j, k]^2 - Nref^2)/2
    c = 1/(steps.y)^2
    
    # d は 右辺用
    d(i,j,k) = 1im*4*k0*matN[i, j, k]/steps.z - 2/steps.y^2 + k0^2(matN[i, j, k]^2 - Nref^2)/2
    B = zeros(ComplexF64,N.x,1)

    for j in 1:N.y
        # Ax=BのA (z = k+1 における係数行列)を作る。
        # mapでベクトルを作っておいて、diagmで対角行列にする。
        # bが位置によって値が異なるのでmapで対応。
        # diagmは配列を正方行列の対角に配置する。
        # pair（=>) でズレを表現できる。
        A = diagm(map(i -> b(i,j,k),1:N.x))
        A += diagm(1 => fill(c,N.x-1))
        A += diagm(-1 => fill(a,N.x-1))

        #透明境界条件(TBC) for A#######
        # 参考文献  左貝潤一 光導波路の電磁界解析 p.145
        #端---------------------
        imKxL = -(1/steps.x)   * log(abs(F_k_before[2,j]/F_k_before[1,j]))
        reKxL = -(1im/steps.x) * log(abs(F_k_before[2,j]/F_k_before[1,j]*exp(imKxL*steps.x)))

        # reKxLおよびreKxRが負の時、
        # 左貝だと符号逆転させて、藪だと0にすると書いてある。前者のほうがそれっぽい気がする。
        if real(reKxL)<0
            reKxL = -1*reKxL      #左貝方式
            # reKxL = 0           #藪方式
        end
        ###!!!ηLの計算違う？左貝144p。どっちが正解？１
        #ηL = exp(1im* reKxL * steps.y - imKxL*steps.x)
        ηL = exp(1im* (reKxL+ imKxL) *steps.x)
        
        ###!!! ここA[1,1]か？ ← 1,1でいい
        ###!!! 左貝(7.25b)と(7,37a)比較して 差分をとる。
        A[1,1] -= ηL/(steps.x)^2

        #右端---------------------
        imKxR = (1/steps.x) * log(abs(F_k_before[N.x, j]/F_k_before[N.x-1, j]))
        reKxR = (1im/steps.x) * log(abs(F_k_before[N.x-1, j]/F_k_before[N.x, j]*exp(-imKxR*steps.y)))

        if real(reKxR)<0
            reKxR = -reKxR      #左貝方式
            # reKxR = 0           #藪方式
        end
        ηR = exp(1im* (reKxR + imKxR) * steps.y)
        A[N.x, N.x] -= ηR/(steps.x)^2

        #########################

        
        # 透明境界条件(TBC) for B#######
        colBL = (2-ηL)/steps.x^2 - (matN[1, j, k]^2-Nref^2)*k0^2 + (4im*Nref*k0)/steps.z
        colBR = (2-ηR)/steps.x^2 - (matN[N.x, j, k]^2-Nref^2)*k0^2 + (4im*Nref*k0)/steps.z
        colC = -1/(steps.y)^2

        #@show F_k_before

        for i in 1:N.x
            if j == 1
                B[i] = -conj(colBL)*F_k_before[i,1] - colC*F_k_before[i,2]
            elseif j == N.y
                B[i] = -conj(colBR)*F_k_before[i,N.y] - colC*F_k_before[i,N.y-1]
            else
                B[i] = c*F_k_before[i, j] + d(i, j, k)*F_k_before[i,j-1] + d(i,j,k)* F_k_before[i,j+1]
            end
        end

        #########################

        F_k_after[:,j] = B\A
        #F_k_after[]を使って、新しい屈折率マップを作る。

    end
end


function retN()
    return 1
end

# 試験的。BPMのモード測定を使う
function showMode()
end

# 電界の初期条件
function initial_set(Mode, Beamparam ,F0)
    KTOptical.setParam(Beamparam.w, 0, Beamparam.wavelength)

    x = range(-crange.x/2, crange.x/2 - steps.x, length = N.x)
    y = range(-crange.y/2, crange.y/2 - steps.y ,length = N.y)
    E = LG_E.(Mode.l, Mode.p, x, y')
    #plot()
    return E

end

function renewN!(matN, E3d)
end

function main()
    #セル個数
    F0 = zeros(N.x, N.y, N.z)
#    @show size(F0)
    E = initial_set(mode, beam, F0)
    gr()
    #!!!!!!!!!!!!!!!!!!!!!!!!
    x = range(-crange.x/2, crange.x/2 ,step = steps.x)
    y = range(-crange.y/2, crange.y/2 ,step = steps.y)
    z = range(-crange.z/2, crange.z/2 ,step = steps.z)

    F_result = zeros(ComplexF64,(N.x,N.y,N.z));
    #F_k_1stは現在のF_k
    #F_k_2ndは更新されたF_k
    @show N.x * N.y, N.x* N.y
    F_k_1st = zeros(ComplexF64, N.x, N.y) #Zerosに型指定忘れないこと！！！
    F_k_2nd = zeros(ComplexF64, N.x, N.y)
    matN = zeros(ComplexF64,(N.x,N.y,N.z))
    @show size(F0)
    @show size(E)
    @show size(F_k_1st[:,:,1])
    F_k_1st[:,:,1] = E
    setNwaveguide!(matN, steps.x, steps.y, steps.z, 5um, 0, material.nb, material.nb + material.Δn0, 0.5)
    #@show matN
    
    #左辺は更新されたF_k_1stが入る。
    # S
    #初期条件を F_K_1st に入れる。


    for t in 1:N.t
        @show "Zmax:", N.z
        for k in 1:N.z
            # x 固定、　y方向移動
            @show t,k
            calcStep1!(F_k_1st, F_k_2nd, k, matN)
            # y 固定、　x方向移動
            calcStep2!(F_k_2nd, F_k_1st, k, matN)
            F_result[:,:,k] = F_k_1st
        end
    #    @save "/savefile/F_"* string(t)*".jld2" F_k_2nd
    end
    # @show F_result
    println("calc done")
    println("plotting.....") 
    Ezx = abs.(F_result[2, :, :]).^2
    @show Ezx
    p1 = contourf(x, z, Ezx )
    display(p1)
    #//////////////////////////////////////////////////////////////////////


end

@time main()

