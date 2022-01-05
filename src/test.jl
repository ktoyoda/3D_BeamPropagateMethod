
a = [1 -1 2;2 -2 3; 3 -1 0]

@show a

b = map(a) do x
    if x < 0
        return -x
    else
        return x
    end
end

@show b