include("Params.jl")

using KTOptical

"""
'IntensityToN(Iintegral, baseN::Float64, deltaN::Float64, t::Float64, τ::Float64, Uᵨ::Float64)'
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
function IntensityTodN(Iintegral, Et, E, mat, t_now::Float64, step,ratio)
    
    xy_square = step.x * step.y
    I = abs2(E)
    if I /xy_square > mat.U/xy_square # 面積強度が臨海強度を超えたなら
        if Et == 0
            Et = t_now # 臨海強度を超えた時の時間を記録
        end  
    end
    # 積算
    if Et + mat.τ < t_now
        Iintegral += I*step.t
    end
    if Iintegral == NaN
        Iintegral = 0
    end
    if mat.Δn0*(1-ratio)*(1-exp((-1/mat.U)*Iintegral)) == NaN
        @show mat.Δn0
        @show ratio
        @show mat.U
        @show Iintegral
    end
    return mat.Δn0*(1-ratio)*(1-exp((-1/mat.U)*Iintegral))
end

function Smoothing(N::Array{Float64,3}, x, y, z, smooting::Int64)
    temp_N = copy(N)
    lenx = length(x)
    leny = length(y)
    lenz = length(z)
    calcrange = smooting ÷ 2

    for (i, x_) in enumerate(x)
        for (j, y_) in enumerate(y)
            for (k, z_) in enumerate(z)
                edge_lx, edge_ly, edge_lz = i-calcrange, j-calcrange, k-calcrange 
                edge_rx, edge_ry, edge_rz = i+calcrange,j+calcrange,k+calcrange
#                @show i,j,k,calcrange
                if 1 > i-calcrange
                    edge_lx = 1
                end
                if lenx <= calcrange + i 
                    edge_rx = lenx
                end 

                if 1 > j- calcrange 
                    edge_ly = 1
                end

                if leny < calcrange + j 
                    edge_ry = leny
                end 
                
                if 1 > k- calcrange 
                    edge_lz = 1
                end
                if lenz <= calcrange + k 
                    edge_rz = lenx
                end 
                
                temp_N[i,j,k] = sum(N[edge_lx: edge_rx, edge_ly: edge_ry, edge_lz: edge_rz]) / 
                        length(N[edge_lx: edge_rx, edge_ly: edge_ry, edge_lz: edge_rz])     
            end
        end
    end
    return temp_N
end

# 一様な光ファイバ
function setNwaveguide!(N::Array{Float64,3},xstep, ystep, zstep, diameter, angle,baseN , deltaN, starting_separate)
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
                        N[i, j, k] = baseN + deltaN
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


# 第二引数にモード(mode)を取る。モードのIntensityに相当した屈折率分布を与える。ピークの屈折率をratio * I(peak) にする
function setNwaveguide!(N::Array{Float64,3},mode, ratio, xstep, ystep, zstep, baseN , deltaN, Beam)


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

    Intensity = KTOptical.I.(Ref(mode),xrange,yrange')

    max_I = maximum(Intensity)
    Intensity = Intensity / max_I * Beam.U0
    @show maximum(Intensity)
    for (i,x_pos) in enumerate(xrange)
        for (j,y_pos) in enumerate(yrange)
            for (k,z_pos) in enumerate(zrange)
                #セパレート前
#                if k < Nz * starting_separate
                    N[i, j, k] = Intensity[i,j] * ratio *  deltaN + baseN
                    
                    #斜め方向にして正射影を取ればいいと思うんだけど、
                    #いまはまっすぐの導波路を作ろう。
                    #if (x_pos-angle*x_pos)
                #end
            end
        end
    end
end