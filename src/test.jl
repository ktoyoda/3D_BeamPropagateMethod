
a = zeros(Int32,5,5)
b = ones(Int32,5,1)
@show a
a[5,:] = b
@show a