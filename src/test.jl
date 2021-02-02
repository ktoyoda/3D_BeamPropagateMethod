x = 5
using LinearAlgebra

A = diagm(fill(2,x))
A += diagm(1 => fill(3,x-1))
A += diagm(-1 => fill(4,x-1))

@show log(2.8)