#problem 2.1
function silnia_rek(N=1)
    if N < 2
        return 1
    else
        return silnia_rek(N-1) * N
    end
end

#problem 2.2
function silnia_iter(N=1)
    result = 1
    for i in 1:N
        if N > 1
            result = i * result;
        else
            return 1
        end
    end
    println(result)
    return result
end

silnia_iter(20)
println(silnia_rek(BigInt))