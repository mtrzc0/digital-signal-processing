# dft, suma amplitud
function rozwiazanie(;
    fp::Int = 682,
    x::Vector{ComplexF64} = ComplexF64[-0.35 + 0.22im, -0.68 - 1.13im, 0.69 + 0.08im, 0.24 - 0.56im, 0.03 + 0.22im, -0.08 + 0.33im, -0.0 + 0.71im, -0.15 + 0.87im, 0.74 - 0.19im, 0.75 + 1.81im, 0.89 - 0.66im, 0.15 - 0.36im, 0.33 - 0.75im, 0.22 - 0.59im, 0.71 - 0.62im, -0.09 + 0.63im, -0.1 + 0.49im, -0.16 - 0.3im, 0.25 + 1.87im, -0.17 + 0.11im, -1.15 - 0.24im, 0.55 + 0.76im, -0.26 + 0.18im, -0.55 - 0.18im, 1.34 - 0.02im, -0.56 - 0.26im, 1.23 + 0.46im, 0.17 + 0.34im, -0.7 - 0.65im, -0.91 - 0.81im, -0.29 - 0.3im],
    f::Vector{Int} = [-176, -154, 66, 110, 176, 198, 242],
)
    N = length(x)
    @show K = mod.(round.(Int, f ./ fp * N), N)
    dft(x) = [sum(x[n]*cispi(-2*(n-1)*(k-1)/length(x)) for n in 1:N) for k in K]
    @show X = dft(x)
    return sum(abs.(X)/N)
end
# rozwiazanie()