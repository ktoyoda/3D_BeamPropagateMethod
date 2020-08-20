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
    xwidth::Float32 = 100um
    ywidth::Float32 = 100um
    zwidth::Float32 = 100um
    trange::Float32 = 10
end
@with_kw struct struct_step
    # セルのステップ（um
    xstep::Float32 = 0.1um
    ystep::Float32 = 0.1um
    zstep::Float32 = 0.1um
    tstep::Float32 = 0.1
end

struct struct_N
    Nx::Int32
    Ny::Int32
    Nz::Int32
    Nt::Int32
end
@with_kw struct struct_materiarl
    nb::Float16 =1.5
    Δn0::Float16 = 0.03
    τ::Float16 = 0.1
    α::Float16 = 0
end
@with_kw struct struct_beam
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