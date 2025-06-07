# fala piloksztalnta, srednia
function rozwiazanie(;
    fp::Float64 = 363.32,
    t1::Float64 = -6.95,
    N::Int = 76,
)
    g(t) = -2*(t-floor(t+1/2))
    x = range(start=t1, step=1/fp, length=N)
    y = [2.4*g(2.4*t-2.7) for t in x]
    return sum(y)/length(y)
end
# rozwiazanie()

# fala piloksztaltna, srednia
function rozwiazanie(;
    fp::Float64 = 355.84,
    t1::Float64 = 2.97,
    N::Int = 198,
)
    g(t) = 2*(t-floor(t+1/2)) #piloksztaltny
    x = range(start=t1, step=1/fp, length=N)
    y = [1.7*g(4.8*t-2.4) for t in x]
    return sum(y)/length(y)
end
# rozwiazanie()

# fala trojkatna
function rozwiazanie(;
    fp::Float64 = 348.34,
    t1::Float64 = 7.9,
    N::Int = 178,
)
    g(t) = 4*abs(t - floor(t+3/4)+1/4) - 1
    x = range(start=t1, step=1/fp, length=N)
    y = [3.5*g(2t-2.6) for t in x]
    return sum(y)/N
end
# rozwiazanie()