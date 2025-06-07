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