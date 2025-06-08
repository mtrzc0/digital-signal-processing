# dft, suma amplitud
function rozwiazanie(;
    fp::Int = 682,
    x::Vector{ComplexF64} = ComplexF64[-0.35 + 0.22im, -0.68 - 1.13im, 0.69 + 0.08im, 0.24 - 0.56im, 0.03 + 0.22im, -0.08 + 0.33im, -0.0 + 0.71im, -0.15 + 0.87im, 0.74 - 0.19im, 0.75 + 1.81im, 0.89 - 0.66im, 0.15 - 0.36im, 0.33 - 0.75im, 0.22 - 0.59im, 0.71 - 0.62im, -0.09 + 0.63im, -0.1 + 0.49im, -0.16 - 0.3im, 0.25 + 1.87im, -0.17 + 0.11im, -1.15 - 0.24im, 0.55 + 0.76im, -0.26 + 0.18im, -0.55 - 0.18im, 1.34 - 0.02im, -0.56 - 0.26im, 1.23 + 0.46im, 0.17 + 0.34im, -0.7 - 0.65im, -0.91 - 0.81im, -0.29 - 0.3im],
    f::Vector{Int} = [-176, -154, 66, 110, 176, 198, 242],
)
    N=length(x)
    sum=0
    for fi in f
        X=0
        for n in 1:N-1
            X+=x[n+1]*cispi(-2*fi/fp*n)
        end
        sum+=abs(X/N)
    end
    return sum
end
# rozwiazanie()
     
function rozwiazanie(;
    fp::Int = 987,
    x::Vector{ComplexF64} = ComplexF64[0.8 - 0.62im, 0.68 + 0.19im, -0.06 + 0.05im, -0.23 + 0.4im, 0.17 - 0.23im, 0.43 + 0.15im, -0.93 - 0.33im, 0.81 + 0.06im, -1.37 - 0.23im, -0.78 + 0.4im, 0.13 - 0.62im, 1.91 - 0.61im, -0.63 + 0.37im, -1.69 - 0.31im, 0.05 - 0.87im, 0.87 + 0.07im, -0.62 + 0.44im, -0.78 - 0.88im, 0.51 - 0.92im, 1.11 + 0.63im, 0.73 - 0.0im, -0.39 + 0.8im, -1.32 - 0.67im, 0.43 + 0.0im, 0.1 - 0.51im, 1.2 + 0.01im, -0.17 - 0.22im, 0.15 + 0.41im, -0.17 + 0.47im, -0.67 + 0.81im, -0.13 - 0.06im, -0.23 + 0.31im, -0.01 + 0.57im, 1.42 - 0.27im, 0.74 - 0.14im, -0.1 + 1.28im, 0.59 - 0.48im, -0.42 - 1.0im, 0.04 - 0.1im, 0.09 + 0.78im, 1.38 - 0.14im, -0.61 + 0.15im, -1.36 - 0.02im, 0.5 - 1.12im, 0.65 - 0.75im, 0.63 - 0.87im, 0.32 + 0.89im],
    f::Vector{Int} = [-441, -420, -63, -42, -21, 210, 462],
)
    N=length(x)
    sum=0
    for fi in f
        X=0
        for n in 1:N-1
            X+=x[n+1]*cispi(-2*(fi/fp)*n)
        end
        sum += angle(X/N)
    end
    return sum
end
# rozwiazanie()

          
function rozwiazanie(;
    fp::Int = 893,
    x::Vector{ComplexF64} = ComplexF64[-0.05 + 0.15im, -0.51 - 0.06im, 0.46 + 0.01im, 0.9 + 0.31im, -0.56 + 0.38im, 0.87 + 0.61im, 0.53 + 0.4im, -1.63 + 1.23im, 0.31 + 0.77im, 1.43 + 0.06im, 0.43 + 0.66im, -1.15 + 0.41im, -0.47 - 0.39im, -1.86 - 0.44im, 1.38 - 0.15im, 0.18 + 0.08im, 0.1 - 1.8im, -0.9 - 0.57im, -0.04 + 1.66im, -0.97 + 1.17im, -0.51 - 1.03im, 0.15 - 0.8im, -0.09 - 0.5im, -0.9 + 1.25im, -0.2 + 1.12im, -0.24 + 1.15im, 1.38 + 0.03im, 0.53 - 0.76im, 0.46 + 0.58im, 0.59 + 0.89im, -1.09 + 0.51im, 0.08 + 0.89im, -0.39 - 0.49im, 0.53 + 0.71im, -0.59 - 0.26im, -0.42 - 0.3im, 0.14 - 0.52im, 0.85 - 0.05im, -0.21 - 0.58im, -0.1 + 0.32im, 0.13 - 1.0im, 1.26 + 0.05im, 0.01 + 0.26im, -0.2 - 1.52im, -1.42 + 0.38im, -0.21 - 1.1im, 0.35 - 0.54im],
    f::Vector{Int} = [133, 171, 323, 418],
)
    N=length(x)
    sum=0
    for fi in f
        X=0
        for n in 1:N-1
            X+=x[n+1]*cispi(-2*fi/fp*n)
        end
        sum+=angle(X/N)
    end
    return sum
end
rozwiazanie()