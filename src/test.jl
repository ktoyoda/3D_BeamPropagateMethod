# 1億まで足し続けるプログラム
println("Julia 非関数")
sum = 0
@time for i in 1:100000000
    global sum += i
end

println("Julia 関数")
function sumtest()
    sum = 0
    for i in 1:100000000
        sum += i
    end

    @show sum

end

@time sumtest()