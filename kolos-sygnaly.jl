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

# fala piloksztaltna
function rozwiazanie(;
    fp::Float64 = 435.66,
    t1::Float64 = -4.54,
    N::Int = 905,
)
    g(t)=-2*(t - floor(t+1/2))
    x=range(start=t1,step=1/fp,length=N)
    y=[5.1*g(1.2*t-3.3) for t in x]
    return sum(y)/N
end
# rozwiazanie()

function rozwiazanie(;
    fp::Float64 = 468.7,
    t1::Float64 = -5.99,
    N::Int = 733,
)
    g(t) = sign(sin(2π*t)) # bipolarna fala prostokatna git?
    x=range(start=t1, step=1/fp, length=N)
    y=[1.6*g(4.3*t-1.1) for t in x]
    return sum(y)/N
end
# rozwiazanie()
          
function rozwiazanie(;
    fp::Float64 = 361.14,
    t1::Float64 = -4.16,
    N::Int = 162,
)
    g(t) = 2(t-floor(t+1/2))
    x=range(start=t1, step=1/fp, length=N)
    y=[1.4*g(t-3.2) for t in x]
    return sum(y)/N
end
# rozwiazanie()

function rozwiazanie(;
    fp::Float64 = 471.46,
    t1::Float64 = 0.28,
    N::Int = 922,
)
    g(t) = sign(sin(2π*t))
    x=range(start=t1, step=1/fp, length=N)
    y=[3*g(2.5*t-1) for t in x]
    return sum(y)/N
end
# rozwiazanie()

function rozwiazanie(;
    fp::Float64 = 471.46,
    t1::Float64 = 0.28,
    N::Int = 922,
)
    g(t) = sign(sin(2π*t))
    x=range(start=t1, step=1/fp, length=N)
    y = [3*g(2.5*t-1) for t in x]
    return sum(y)/N
end
# rozwiazanie()

function rozwiazanie(;
    fp::Float64 = 434.88,
    t1::Float64 = 0.92,
    N::Int = 826,
)
    g(t)=-2*(t-floor(t+1/2))
    x=range(start=t1, step=1/fp, length=N)
    y=[3.9*g(0.3*t-2) for t in x]
    return sum(y)/N
end
# rozwiazanie()

function rozwiazanie(;
    fp::Float64 = 498.99,
    t1::Float64 = 9.38,
    N::Int = 504,
)
    g(t)=2*(t-floor(t+1/2))
    x=range(start=t1, step=1/fp, length=N)
    y=[3.4*g(3.6*t-3) for t in x]
    return sum(y)/N
end
# rozwiazanie()

function rozwiazanie(;
    fp::Float64 = 406.56,
    t1::Float64 = -2.67,
    N::Int = 643,
)
    g(t)=4*abs(t+1/4-floor(t+3/4))-1
    x=range(start=t1, step=1/fp, length=N)
    y=[0.6*g(4.7*t-4.7) for t in x]
    return sum(y)/N
end
# rozwiazanie()

function rozwiazanie(;
    fp::Float64 = 352.67,
    t1::Float64 = 0.09,
    N::Int = 237,
)
    g(t)=sign(sin(2π*t))
    x=range(start=t1, step=1/fp, length=N)
    y=[3.3*g(4.1*t-4.2) for t in x]
    return sum(y)/N
end
# rozwiazanie()