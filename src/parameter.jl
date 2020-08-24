#コード分割のため、
#パラメーターはすべてここで定義

# struct の初期化のため、Parameters をインクルードする。

#negative value の累乗は.0を付ける
m = 1
mm = 10^(-3.0m)
um = 10^(-6.0m)
nm = 10^(-9.0m)

# 実験条件####################
module Param
    using Parameters
    @with_kw struct range
        # セルの距離(um)
        xwidth::Float64 = 100um
        ywidth::Float64 = 100um
        zwidth::Float64 = 100um
        trange::Float64 = 10
    end
    @with_kw struct step
        # セルのステップ（um
        xstep::Float64 = 0.1um
        ystep::Float64 = 0.1um
        zstep::Float64 = 0.1um
        tstep::Float64 = 0.1
    end

    @with_kw struct N
        Nx::Int32
        Ny::Int32
        Nz::Int32
        Nt::Int32
    end
    @with_kw struct materiarl
        nb::Float64 =1.5
        Δn0::Float64 = 0.03
        τ::Float64 = 0.1
        α::Float64 = 0
    end
    @with_kw struct beam
        w::Float64 = 3um
        U0::Float64 = 100
        wavelength::Float64 = 1.06um
    end
    @with_kw struct gauss_mode
        m::Int16 = 1
        n::Int16 = 1
    end
    @with_kw struct vortex_mode
        l::Int16 = 1
        p::Int16 = 1
    end
end