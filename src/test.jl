using LinearAlgebra
A = [3 4; 5 6]
B = [1;1]

@show A
@show B

A[2,:] = B
@show A
@show A\B