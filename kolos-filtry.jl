# filtr HP, suma
function rozwiazanie(;
    order::Int = 30,
    fp::Float64 = 192.0,
    f0::Float64 = 49.92,
    z::Vector{Int} = [1, 10, 1, 19, 21],
)
    F0 = f0/fp
    M = order

    δ(n) = (n==0) ? 1 : 0
    w = [0.5 + 0.5cospi(2n/(2M+1)) for n in -M÷2:M÷2]
    h = [(n!=0) ? δ(n) - 2F0*sinc(2*F0*n) : 1-2F0 for n in -M÷2:M÷2]
    return sum(h .* w)
end
# rozwiazanie()

function rozwiazanie(;
    order::Int = 22,
    fp::Float64 = 197.0,
    f0::Float64 = 33.49,
    z::Vector{Int} = [12, 8, 17, 2, 13, 11],
)
    M=order
    F0=f0/fp
    δ(n)= n==0 ? 1 : 0
    w = [1-abs(n)/(M+1) for n in -M÷2:M÷2]
    h = [n!=0 ? δ(n)-2F0*sinc(2*F0*n) : 1-2*F0 for n in -M÷2:M÷2]
    return sum([(h.*w)[i] for i in z])
end
# rozwiazanie()

function rozwiazanie(;
    order::Int = 98,
    fp::Float64 = 194.0,
    f1::Float64 = 17.46,
    f2::Float64 = 67.9,
    z::Vector{Int} = [16, 29, 20, 8],
)
    M=order
    F1=f1/fp
    F2=f2/fp
    w = [0.54 + 0.46*cospi(2n/(M+1)) for n in -M÷2:M÷2]
    h = [n!=0 ? 2F2*sinc(2*F2*n)-2F1*sinc(2*F1*n) : 2*(F2-F1) for n in -M÷2:M÷2]
    return sum([(h.*w)[i] for i in z])
end
rozwiazanie()