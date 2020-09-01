using LinearAlgebra
using Distributed
function test(n)
    b(i) = i^2
    A = diagm(pmap(x -> b(x),1:n))
    A +=diagm(1 => pmap(x -> b(x),1:n-1))
    A +=diagm(-1 => pmap(x -> b(x),1:n-1))
end 


@time test(5000)