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

# silnia_iter(20)
# println(silnia_rek(BigInt))

#problem 2.3
function liczba_parz(num)
    return (num % 2 == 0) ? true : false
end

#problem 2.4
function liczba_pierw(N)
    if N == 1
        return false
    elseif N == 2
        return true 
    end

    is_prime = true 
    for i in 3:2:(N-1)
        if N % i == 0
            is_prime = false
        end
    end
    return is_prime
end

# liczba_pierw(19)

#problem 2.5
function str_rev(str)
    stack = []
    for letter in str
        append!(stack, letter)
    end
    for letter in str
        print(pop!(stack))
    end
end

#str_rev("Kinia")

#problem 2.6
function palindrom(str)
    str = lowercase(str)
    stack = []
    new_str = []
    for letter in str
        append!(stack, letter)
    end
    for letter in str
        append!(new_str, pop!(stack))
    end

    is_pal = true
    for i in 1:lastindex(str)
        if str[i] != new_str[i] 
            is_pal = false
        end
    end

    return is_pal
end

palindrom("ADA")

#problem 2.7
function sierpinski_area(N)
    if N == 0
        return 1
    end
    return (3^N)*(sqrt(3)/4)*((2/(3^(0.25)))*(0.5^N))^2 
end

# sierpinski_area(4)

#problem 2.8
function newton(a, x0, err)
    xk = x0
    xk1 = a
    xk_1(xk) = 0.5*(xk+a/xk)
    i = 1
    while true
        xk1=xk_1(xk)
        xk=xk1
        i+=1

        if (xk^2 - a <= err)
            println("Znaleziono $xk po $i iteracjach")
            return xk
        end
    end
    return xk1
end

# newton(1010, 1, 0.0001)

#problem 2.9
using Plots
function kzbieznosc()
    N = 1000
    xp = LinRange(-1,2,N)
    yp = LinRange(-1,1,N)

    Mk = zeros(Int, N, N)
    for i in 1:N
        for j in 1:N
            zn = 0 + 0*im
            k = 1
            p = xp[i] + im * yp[j]
            while abs(zn) < 2 && k < 200
                zn = zn^2 + p 
                k+=1
            end
            Mk[i, j] = k
        end
    end
    heatmap(xp, yp, Mk)
    # surface(xp, yp, Mk, color=:viridis)
end

kzbieznosc()

#problem 2.10

#problem 2.11