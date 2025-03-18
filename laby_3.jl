using Plots
include("CPS.jl")

#problem 3.1
function vector256(fs, ts, A, f, phase)
    return [A*sin(2*pi * f * (ts + (t-1)/fs) + phase) for t in 1:256]
end

# vector256(1000, 0.25, 2, 25, pi/4)

#problem 3.2
function vectorN(fs, A, f, phase, ts, te)
    return [A*exp(im*2*pi * f * (ts + (t-1)/fs) + phase) for t in 1:((te-ts)*fs + 1)]
end

# vectorN(2048, 0.25, pi/2, pi, 5, 10)

#problem 3.7

function ramp_test()
    n=1000
    stop=10
    y::Vector{Real} = []
    t=0:(stop/n):stop
    for i in 1:(n+1)
        push!(y, CPS.ramp_wave(t[i]))
    end
    plot(t, y)
end

ramp_test()

#problem 3.8
function sawtooth_test()
    n=1000
    stop=10
    y::Vector{Real} = []
    t=0:(stop/n):stop
    for i in 1:(n+1)
        push!(y, CPS.sawtooth_wave(t[i]))
    end
    plot(t, y)
end

sawtooth_test()