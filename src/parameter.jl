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
    xwidth::Int16 = 100
    ywidth::Int16 = 100
    zwidth::Int16 = 100
    trange::Int16 = 50
end
@with_kw struct struct_step
    xstep::Float16 = 0.1
    ystep::Float16 = 0.1
    zstep::Float16 = 0.1
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