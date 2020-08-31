

#ビームの集光をテーパー型の屈折率分布として表現する。
function setNfocus(N::Array{ComplexF64,2},NA,n,separatingpoint,)
    
end

function setNwaveguide!(N::Array{ComplexF64,3},xstep, ystep, zstep, diameter, angle,baseN , propN, starting_separate)
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
                if k < Nz * starting_separate
                    if x_pos^2 + y_pos^2 < diameter
                        N[i, j, k] = propN
                    end
                end
                if k >= Nz * starting_separate
                    if x_pos^2 + y_pos^2 < diameter
                        N[i, j, k] = propN
                    end
                    #斜め方向にして正射影を取ればいいと思うんだけど、
                    #いまはまっすぐの導波路を作ろう。
                    #if (x_pos-angle*x_pos)
                end
            end

        end
    end


end