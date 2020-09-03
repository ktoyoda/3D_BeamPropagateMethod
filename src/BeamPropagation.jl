#計算に関係ないパッケージ
using Pkg

include("parameter.jl")
include("propagator.jl")
#計算に使うパッケージ
using LinearAlgebra
using Plots
using KTOptical
using .Param
using BenchmarkTools
using JLD2

#using FileIO
# 計算条件###################
#計算レンジ
um = Param.um
crange = Param.crange(x = 50um, y = 50um, z = 100um, t = 10)
#計算ステップ
step = Param.step(x = 0.25um, y = 0.25um, z = 0.25um)
Nx = Int(crange.x / step.x)
Ny = Int(crange.y / step.y)
Nz = Int(crange.z / step.z)
Nt = Int(crange.t / step.t)
N = Param.N(Nx,Ny,Nz,Nt)

#材料情報
material = Param.materiarl(nb = 1.5, Δn0 = 0.01, τ = 0.1,α = 0)
mode = Param.vortex_mode(l=1, p=0)

#ビーム情報
beam = Param.beam(w = 20um, U0 = 100.0, wavelength = 1.064um)

println("計算環境")
#versioninfo()
Pkg.status()

println("計算条件")
@show crange
@show step
@show N
@show material
@show mode
@show beam

println("//////////////////////////////")
function returnPades(T::Integer,N::Integer)
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
# 各A,Bに入る屈折率項は座標情報が必要である。
# Aに入るのは()
# 
##F_k11 既知のビーム伝搬
#F_kp12 未知のビーム伝搬の格納先(F_k+1/2)
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
# 再度回していく。

function calcStep1!(F_k11, F_kp12,k,matN)
    k0 = 2π / beam.wavelength
    a = -1/(step.y)^2
    b(i,j,k) = 1im*4*k0*matN[i, j, k]/step.z + 2/step.y^2 - k0^2(matN[i, j, k]-1.35^2)/2
    c = 1/(step.x)^2
    d(i,j,k) = 1im*4*k0*matN[i, j, k]/step.z - 2/step.x^2 + k0^2(matN[i, j, k]-1.35^2)/2
    B = zeros(ComplexF64,N.y,1)

    for i in 1:N.x
        # Ax=BのA, z = k+1を作る。
        # mapでベクトルを作っておいて、diagmで対角行列にする。
        # bが位置によって値が異なるのでmapで対応。

        A = diagm(map(j -> b(i,j,k),1:N.y))
        A += diagm( 1 => map(j -> b(i,j,k),1:N.y-1))
        A += diagm(-1 => map(j -> b(i,j,k),1:N.y-1))

        #透明境界条件(TBC) for A#######
        #左端---------------------
        imKxL = -(1/step.y)   * log(abs(F_k11[i,2]/F_k11[i,1]))
        reKxL = -(1im/step.y) * log(abs(F_k11[i,2]/F_k11[i,1]*exp(imKxL*step.y)))
        if real(reKxL)<0
            reKxL = -1*reKxL
        end
        ηL = exp(1im* reKxL * step.y - imKxL*step.y)

        A[1,1] -= ηL/(step.y)^2

        #右端---------------------
        imKxR = (1/step.y)   * log(abs(F_k11[i,N.y-2]/F_k11[i,N.y-1]))
        reKxR = (1im/step.y) * log(abs(F_k11[i,N.y-2]/F_k11[i,N.y-1]*exp(-imKxR*step.y)))
        if real(reKxR)<0
            reKxR = -reKxR
        end
        ηR = exp(1im* reKxR * step.y - imKxR* step.y)
        A[N.y,N.y] -= ηR/(step.y)^2
        #########################
        for j in 2:N.y-1
            B[i] = c*F_k11[N.y * j + i] + d(i,j,k)*F_k11[N.y*j + i-1] + d(i,j,k)* F_k11[N.y*j + i+1]
        end

        # 透明境界条件(TBC) for B#######
        # 参考文献 ●●● p.xxx
        colBL = (2-ηL)/step.y^2 - (matN[i,1,k]-1.35^2)*k0^2 + (4im*matN[i,1,k]*k0)/step.z
        colBR = (2-ηR)/step.y^2 - (matN[i,N.y,k]-1.35^2)*k0^2 + (4im*matN[i,N.y,k]*k0)/step.z
        colC = -1/(step.x)^2
        B[1] = -conj(colBL)*F_k11[i,1] - colC*F_k11[i,2]
        B[N.x] = -conj(colBR)*F_k11[i,N.y] - colC*F_k11[i,N.y-1]
        #########################
        F_kp12[i,:] = B\A



    end
    return F_kp12
end

#ADIの未知数X方向 定数Y方向 差分
function calcStep2!(F_k12, F_kp21,k, matN)
    #正確には、k+1とKの屈折率の平均を取るべきだと思うが、
    #今のところはk+1を抜き出す   
    return F_kp21
end

function retN()
    return 1
end

# 試験的。BPMのモード測定を使う
function showMode()
end

# 電界の初期条件
@time function initial_set(M, B ,F0)
    KTOptical.setParam(B.w, 0, B.wavelength)
    x = range(-crange.x/2, crange.x/2 ,step = step.x)
    y = range(-crange.y/2, crange.y/2 ,step = step.y)

    E = LG_I.(M.l, M.p,x,y')
    return E

end

# calcstep
function renewN!(matN, E)
    if size(matN) == size(E)
        println("matN and E have different sizes.[renewE] ")
        return 0
    end

    for matN 

end

function main()
    #セル個数
    F0 = zeros(N.x, N.y, N.z)
    E = initial_set(mode, beam, F0)
    gr()
    x = range(-crange.x/2, crange.x/2 ,step = step.x)

    F_result = zeros(Float64,(N.x,N.y,N.z));
    #F_k_1stは現在のF_k
    #F_k_2ndは更新されたF_k
    @show N.x * N.y, N.x* N.y
    F_k_1st = zeros(N.x, N.y)
    F_k_2nd = zeros(N.x, N.y)
    matN = zeros(ComplexF64,(N.x,N.y,N.z))
    typeof(matN)
    setNwaveguide!(matN, step.x, step.y, step.z, 4um, 0, material.nb, material.nb + material.Δn0, 0.5)
    #左辺は更新されたF_k_1stが入る。
    # 

    for t in 1:N.t
        for k in 1:N.z-1
            calcStep1!(F_k_1st, F_k_2nd, k, matN)
            calcStep2!(F_k_2nd, F_k_1st, k, matN)
        end
    #    @save "/savefile/F_"* string(t)*".jld2" F_k_2nd
    end
    "done"
end

@time main()