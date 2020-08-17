#using Unitful

m = 1
mm = 10^-3m
um = 10^-6m
nm = 10^-9m


# 実験条件####################
const a = 10nm 
const np = 1.53 # 粒子の屈折率
const nb = 1.32 # 媒体の屈折率
const ρ0 = 0.2  # 摂動がないときの粒子密度(単位確認)
const ε0 = 2.4 # 誘電率(単位確認)
const l = 1 # とぽろじかるちゃーじ
const w = 5
const λ = 1.06
# 計算条件###################
#計算レンジ
const xwidth = 100
const ywidth = 100
const zwidth = 3
const trange = 50
#計算ステップ
const xstep = 0.1
const ystep = 0.1
const zstep = 0.05
const D = 0
##/////
const φ0 = 1
const kbT = 1.38064852 * 10^ -23
const Vp = 4π*a^3/3
const m = np/nb
const α = 3Vp*ε0*nb^2*((m^2-1)/(m^2+2))
##################
x = -xwidth/2:xstep:xwidth/2
y = -ywidth/2:ystep:ywidth/2
const xn = length(x)
const yn = length(y)

D = 0
Φ = deg2rad(D)
rot_axis(x,y,Φ) = x*cos(Φ) - y*sin(Φ)
p = rot_axis.(x',  y , Φ)
q = rot_axis.(y , -x', Φ)
r = .√(p.^2+q.^2)
φ = @. (1/√(factorial(abs(l))))*((√2*r/w)^abs(l))*exp(-r^2/w^2)
φ2 = abs(φ0)^2
σ = ((128π^5*a^2*nb^4)/3)*(a/λ)^4*((m^2-1)/(m^2+2))^2
γ = σ*ρ0*exp(α*φ2/(4kbT))

@show σ
#@show dimension(σ)
