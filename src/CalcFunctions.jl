#計算に使うパッケージ
using LinearAlgebra
using Plots
using KTOptical 
# 自作のパッケージを使うときにはGithubにアップロードして、
# GithubのHTTPSをadd HTTPS_URLする。
# プライベートでよい。

using BenchmarkTools
using JLD2
using Polynomials

function returnPades(T::Integer,N::Integer)

end

# 各Zに対して、calcstepを１，２の順で回していく。
# step 1 は、
# ADIの未知数Y方向(j) 定数X方向(i) 差分
# 
# Ax = B をつくってトーマスアルゴリズムJFで解く
# juliaの場合はA\Bでよい。(バックスラッシュ）

###### 各A,Bに入る屈折率項は座標情報が必要である。
# これはmatNとしてある。非線形の場合は、tのスイープで更新する。
# ビーム伝搬法をダグラスガンの三次元
# calcStep1では
# F_k_before 既知のビーム伝搬
# F_k_half  未知のビーム伝搬の格納先(F_k+1/2)
# calcStep2では
# F_k_half 既知のビーム伝搬
# F_k_after  未知のビーム伝搬の格納先(F_k+1)
# k 今のZカウント
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

# テキスト中nr と nについて
# 

#ADIの未知数X方向 定数Y方向 差分
function calcStep1!(F_k_before, F_k_half, k, matN, Nref)
    #正確には、k+1とKの屈折率の平均を取るべきだと思うが、
    #今のところはk+1を抜き出す   
    k0 = 2π / beam.wavelength
    # a,b,cは左辺用
    ax = -1/(steps.x)^2
    b(i,j,k) = 1im*4k0*Nref/steps.z + 2/steps.x^2 - k0^2(matN[i, j, k]^2 - Nref^2)/2
    cx = -1/(steps.x)^2
    
    #Bをつくる。 d は 右辺用 F_k_2ndの係数
    ay = +1/(steps.y)^2
    d(i,j,k) = 1im*4*k0*Nref/steps.z - 2/steps.y^2 + k0^2(matN[i, j, k]^2 - Nref^2)/2
    cy = +1/(steps.y)^2
    B = zeros(ComplexF64,N.x,1)

    for j in 1:N.y
        #! jに関して透明境界条件必要じゃないんか？
        #! 左貝146Pにφ₀ʳ≒φ₁ʳ*ηLで近似できると書いてあった。
        # Ax=BのA (z = k+1 における係数行列)を作る。
        # mapでベクトルを作っておいて、diagmで対角行列にする。
        # bが位置によって値が異なるのでmapで対応。
        # diagmは配列を正方行列の対角に配置する。
        # pair（=>) でズレを表現できる。
        A = diagm(map(i -> b(i,j,k),1:N.x))
        A += diagm(1 => fill(cx,N.x-1))
        A += diagm(-1 => fill(ax,N.x-1))

        #透明境界条件(TBC) for A#######
        # 参考文献  左貝潤一 光導波路の電磁界解析 p.145
        #左端 ηLを作成する。---------------------
        #左貝(7.29)
        #=
        imKxL = -(1/steps.x)   * log(F_k_before[2,j]/F_k_before[1,j])
        reKxL = -(1im/steps.x) * log(F_k_before[2,j]/F_k_before[1,j]*exp(imKxL*steps.x))
        KxL = reKxL + imKxL
        =#
        KxL = -(1/(1im*steps.x))   * log(F_k_before[1,j]/F_k_before[2,j])

        if real(KxL) < 0
            #KxL = -reKxL + imKxL
            KxL = -real(KxL) + 1im*imag(KxL)
        end

        ηL = exp(1im*KxL*(-steps.x))
        A[1,1] += ηL*ax

        #右端 ηRを作成する。---------------------
        # 左貝(7.35)
        #=
        imKxR = (1/steps.x) * log(F_k_before[N.x, j]/F_k_before[N.x-1, j])
        reKxR = (1im/steps.x) * log(F_k_before[N.x, j]/F_k_before[N.x-1, j]*exp(-imKxR*steps.y))
        KxR = reKxR + 1im* imKxR
        =#
        KxR = (-1/(1im*steps.x)) * log(F_k_before[N.x, j]/F_k_before[N.x-1, j])
        if real(KxR) > 0
            KxR = -real(KxR) + 1im* imag(KxR)      #左貝方式
        end

        ηR = exp( 1im* KxR * steps.x)
        A[N.x, N.x] += ηR*cx

        #########################

        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
         # 透明境界条件(TBC)を適用したBの係数#######
         d_b(i,j,k) = - 2/steps.y^2 + (matN[i, j, k]^2-Nref^2)*k0^2 + (4im*Nref*k0)/steps.z
         d_ac = 1/(steps.y)^2

        #@show F_k_before

        for i in 1:N.x
            if j == 1
                K_temp = (-1/(1im*steps.y))   * log(F_k_before[i,1] / F_k_before[i,2])
                F0 = F_k_before[i,1] * exp(1im * K_temp*(-steps.y))
                B[i] = d_b(i,j,k) * F_k_before[i,j] + d_ac * F_k_before[i,j+1] +d_ac*F0 #+ ηU？
            elseif j == Ny
                K_temp = (-1/(1im*steps.y))   * log(F_k_before[i,Ny] / F_k_before[i,Ny-1])
                F0 = F_k_before[i,Ny] * exp(1im * K_temp*(-steps.y))
                B[i] = d_b(i,j,k) * F_k_before[i,j] + d_ac * F_k_before[i,j-1] +d_ac * F0#+ ηB？
            else
                B[i] = d_b(i,j,k) * F_k_before[i,j] + d_ac * F_k_before[i,j-1] + d_ac*F_k_before[i,j+1]    
            end
        end
        #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        #########################
        F_k_half[:,j] = A \ B
        #F_k_after[]を使って、新しい屈折率マップを作る。

    end
end

function calcStep2!(F_k_half, F_k_next, k, matN, Nref)

    k0 = 2π / beam.wavelength
    # Aをつくる。 a,b,cは左辺用
    ay = -1/(steps.y)^2
    b(i,j,k) = 1im*4*k0*Nref/steps.z + 2/steps.y^2 - k0^2(matN[i, j, k]^2 - Nref^2)/2
    cy = -1/(steps.y)^2
    
    # Bをつくる。 d は 右辺用。 F_k_halfの係数
    ax = 1/(steps.x)^2
    d(i,j,k) = 1im*4*k0*Nref/steps.z - 2/steps.x^2 + k0^2(matN[i, j, k]^2 - Nref^2)/2
    cx = 1/(steps.x)^2
    B = zeros(ComplexF64,N.y,1)
    # yをつくる。

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
        A += diagm(1 => fill(cy,N.y-1))
        A += diagm(-1 => fill(ay,N.y-1))
        
        #透明境界条件(TBC) for A#######--------------------------------------------------
        # 参考文献  左貝潤一 光導波路の電磁界解析 p.145
        #左端---------------------
        # 藪107p
        #= ```
        KxL = -1 / (1im*steps.y) * log(F_k_half[i,1]/F_k_half[i,2])
        if real(KxL)<0
            KxL = imag(KxL)
        end

        ηL = exp(1im*KxL*(-steps.y))
        A[1,1] += ηL*ax
        ```
        =# 
        #透明境界条件(TBC) for A#######
        # 参考文献  左貝潤一 光導波路の電磁界解析 p.145
        #上端---------------------
        #左貝(7.29)
        #=
        imKyU = -(1/steps.y)   * log(F_k_half[i,2]/F_k_half[i,1])
        reKyU = -(1im/steps.y) * log(F_k_half[i,2]/F_k_half[i,1]*exp(imKyU*steps.y))
        KyU = reKyU + 1im* imKyU
        =#
        KyU = -(1/(1im*steps.y))   * log(F_k_half[i,1]/F_k_half[i,2])
        if real(KyU) < 0
        #    KyU = -reKyU + 1im * imKyU
           KyU = -real(KyU) +1im * imag(KyU)
        end

        ηU = exp(1im * KyU * (-steps.y))
        A[1,1] += ηU*ay
        ##!!!!!


        #下端---------------------
        #左貝7.35
        KyB = (-1/(1im*steps.y))   * log(F_k_half[i,N.y] / F_k_half[i,N.y-1])

        if real(KyB) > 0
            KyB = -real(KyB) + 1im* imag(KyB)           #藪方式
        end
        
        ηB = exp(1im* KyB * steps.y)
        A[N.y,N.y] += ηB*cy

        #########################------------------------------------------------------
         # 透明境界条件(TBC)を適用したBの係数#######
         d_b(i,j,k) = - 2/steps.x^2 + (matN[i, j, k]^2-Nref^2)*k0^2 + (4im*Nref*k0)/steps.z
         d_ac = 1/(steps.x)^2

        for j in 1:N.y
            if i == 1
                K_temp = (-1/(1im*steps.x))   * log(F_k_half[1,j] / F_k_half[2,j])
                F0 = F_k_half[1,j] * exp(1im * K_temp*(-steps.x))
                B[j] = d_b(i,j,k) * F_k_half[i,j] + d_ac * F_k_half[i+1,j] + d_ac*F0

            elseif i == N.x
                K_temp = (-1/(1im*steps.x))   * log(F_k_half[Nx,j] / F_k_half[Nx-1,j])
                F0 = F_k_half[Nx,j] * exp(1im * K_temp*(-steps.x))
                B[j] = d_b(i,j,k) * F_k_half[i,j] + d_ac * F_k_half[i-1,j] + d_ac*F0 

            else
                B[j] = d_b(i,j,k) * F_k_half[i,j] + d_ac * F_k_half[i-1,j] + d_ac*F_k_half[i+1,j]    
            end
        end

        #########################

        F_k_next[i,:] = A \ B
    end
end

function retN()
    return 1
end

# 試験的。BPMのモード測定を使う
function showMode()
end

# 多重ディスパッチでvortexかgaussかを切り替える
# 電界の初期条件
function initial_set(Mode::vortex_mode, Beamparam ,F0)
    KTOptical.setParam(Beamparam.w, 0, Beamparam.wavelength)

    x = range(-crange.x/2, crange.x/2 - steps.x, length = N.x)
    y = range(-crange.y/2, crange.y/2 - steps.y ,length = N.y)

    #ここをLG_E HG_Eののなんかラッパ的な奴にできない？
    E = LG_E.(Mode.l, Mode.p, x, y')
    #plot()
    return E
end

# 電界の初期条件
function initial_set(Mode::gauss_mode, Beamparam ,F0)
    KTOptical.setParam(Beamparam.w, 0, Beamparam.wavelength)

    x = range(-crange.x/2, crange.x/2 - steps.x, length = N.x)
    y = range(-crange.y/2, crange.y/2 - steps.y ,length = N.y)

    #ここをLG_E HG_Eののなんかラッパ的な奴にできない？
    E = HG_E.(Mode.m, Mode.n, x, y')
    #plot()
    return E
end


function renewN!(matN, E3d)
end
