include("CPS.jl")

using Plots
using Random
using LinearAlgebra

#problem 5.1
function quantize_test()
    L=[1, 2, 3]
    f = CPS.quantize(L)

    println(f(2.6))
end
quantize_test()