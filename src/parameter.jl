#コード分割のため、
#パラメーターはすべてここで定義

# struct の初期化のため、Parameters をインクルードする。
using Parameters

#negative value の累乗は.0を付ける
m = 1
mm = 10^(-3.0m)
um = 10^(-6.0m)
nm = 10^(-9.0m)

# 実験条件####################
@with_kw struct struct_range
    # セルの距離(um)
    xwidth::Int16 = 100um
    ywidth::Int16 = 100um
    zwidth::Int16 = 100um
    trange::Int16 = 50um
end
@with_kw struct struct_step
    # セルのステップ（um
    xstep::Float16 = 0.1um
    ystep::Float16 = 0.1um
    zstep::Float16 = 0.1um
    tstep::Float16 = 0.1um
end

struct struct_N
    Nx::Int
    Ny::Int
    Nz::Int
    Nt::Int
end
@with_kw struct struct_materiarl
    nb::Float16 =1.5
    Δn0::Float16 = 0.03
    τ::Float16 = 0.1
    α::Float16 = 0
end
@with_kw struct struct_beam_param
    w::Float16 = 3um
    U0::Float16 = 100
    wavelength::Float16 = 1.06um
end
@with_kw struct struct_gauss_mode
    m::Int16 = 1
    n::Int16 = 1
end
@with_kw struct struct_vortex_mode
    l::Int16 = 1
    p::Int16 = 1
end