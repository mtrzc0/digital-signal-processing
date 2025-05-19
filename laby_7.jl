include("CPS.jl")

using Plots
using Random
using LinearAlgebra

#problem 7.3
function dtft_test()
    f = 1
    fs = 2
    n = 0:10
    x_sin = sin.(2*pi*f*n) # Sygnał sinus
    w_sin = LinRange(-pi, pi, 1000)
    X_sin = CPS.dtft(f, x_sin, fs)

    println(w_sin, X_sin)
end
# dtft_test()

#problem 7.4
function dft_test()
    x = [0, -3, 0, 3]
    X = CPS.dft(x)
    println(X)
    println(x)
end
# dft_test()

#problem 7.5
function idft_test()
    X = [0, im*6, 0, -im*6]
    x = CPS.idft(X)
    println(x)
end
idft_test()
