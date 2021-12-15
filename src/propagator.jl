include("Params.jl")

"""
'IntensityToN(matN, E, baseN::Float64, deltaN::Float64, t::Float64, τ::Float64, Uᵨ::Float64)'
Yarivの論文
A. S. Kewitsch and A. Yariv, "Self-focusing and self-trapping of optical beams upon photopolymerization," Opt. Lett. 21(1), 24 (2008).
をベースに""屈折率変化""を返す。
E はその時点の電場強度
Iintegral はここまでの電場強度二乗の和(外部で保持しているのを記録) 
baseN は初期屈折率
deltaN は屈折率変化量
t は 現在の時間
stept は時間ステップ
τ は モノマーのラジカル寿命。積分の遅延時間。金析出であれば長めにとる？ 
Uᵩ は 反応開始閾値
"""
function IntensityTodN!(Iintegral, Et, E, mat, t_now::Float64, step)
    xy_square = step.x * step.y
    I = abs2(E)
    if I  > mat.U # 面積強度が臨海強度を超えたなら
        if !(Et == 0)
            Et = t_now # 臨海強度を超えた時の時間を記録
        end  
    end
    # 積算
    if Et + mat.τ > t_now
        Iintegral += I
    end

    return mat.Δn0*(1-exp(-1/mat.U*Iintegral))
end

function setNwaveguide!(N::Array{Float64,3},xstep, ystep, zstep, diameter, angle,baseN , propN, starting_separate)
    Nx = size(N,1)
    Ny = size(N,2)
    Nz = size(N,3)
    xrange_max = xstep * Nx / 2
    yrange_max = ystep * Ny / 2
    zrange_max = zstep * Nz / 2
    xrange = range(-xrange_max,xrange_max,length = Nx)
    yrange = range(-yrange_max,yrange_max,length = Ny)
    zrange = range(-zrange_max,zrange_max,length = Nz)
    fill!(N,baseN)
    for (i,x_pos) in enumerate(xrange)
        for (j,y_pos) in enumerate(yrange)
            for (k,z_pos) in enumerate(zrange)
                #セパレート前
#                if k < Nz * starting_separate

                    if x_pos^2 + y_pos^2 < (diameter/2)^2 
                        N[i, j, k] = propN
                    else 
                        N[i, j, k] = baseN
                    end
                    #斜め方向にして正射影を取ればいいと思うんだけど、
                    #いまはまっすぐの導波路を作ろう。
                    #if (x_pos-angle*x_pos)
                #end
            end
        end
    end
end