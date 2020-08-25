#コード分割のため、
#パラメーターはすべてここで定義

# struct の初期化のため、Parameters をインクルードする。

#negative value の累乗は.0を付ける

# 実験条件####################
module Param
   export  m
   export  mm
   export  um
   export  nm 


    m = 1
    mm = 10^(-3.0m)
    um = 10^(-6.0m)
    nm = 10^(-9.0m)


    using Parameters
    @with_kw struct crange{AbstractFloat}
        # セルの距離(um)
        x::AbstractFloat = 100um
        y::AbstractFloat = 100um
        z::AbstractFloat = 100um
        t::AbstractFloat = 10
    end
    @with_kw struct step{AbstractFloat}
        # セルのステップ（um
        x::AbstractFloat = 0.1um
        y::AbstractFloat = 0.1um
        z::AbstractFloat = 0.1um
        t::AbstractFloat = 0.1
    end

    @with_kw struct N{Integer}
        x::Integer
        y::Integer
        z::Integer
        t::Integer
    end
    @with_kw struct materiarl{AbstractFloat,Number}
        nb::AbstractFloat =1.5
        Δn0::AbstractFloat = 0.03
        τ::AbstractFloat = 0.1
        α::Number = 0
    end
    @with_kw struct beam{AbstructFloat}
        w::AbstractFloat = 3um
        U0::AbstractFloat = 100.0
        wavelength::AbstractFloat = 1.06um
    end
    @with_kw struct gauss_mode{Inteteger}
        m::Integer = 1
        n::Integer = 1
    end
    @with_kw struct vortex_mode{Integer}
        l::Integer = 1
        p::Integer = 1
    end
end