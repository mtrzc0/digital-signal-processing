include("CPS.jl")

using Plots
using Random
using LinearAlgebra

#problem 4.6
function running_mean_test()
    x::Vector{Real}=1:10

    running_mean=CPS.running_mean(x, 1)
    println(running_mean)
    println(length(running_mean))
end
# running_mean_test()

#problem 4.7
function running_energy_test()
    x::Vector{Real}=1:10

    running_mean=CPS.running_energy(x, 1)
    println(running_mean)
    println(length(running_mean))
end
# running_energy_test()

#problem 4.8
function running_power_test()
    x::Vector{Real}=1:10

    running_mean=CPS.running_power(x, 1)
    println(running_mean)
    println(length(running_mean))
end
# running_power_test()

#problem 5.1
function quantize_test()
    L=[1:10]
    CPS.quantize(L)
end