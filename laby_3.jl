include("CPS.jl")

using Plots
using Random

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

#problem 3.3
function white_noise()
    return [(1/sqrt(2))*randn() for t in 1:1000]
end

function white_noise_test()
    t = [1:1000]
    plot(t, white_noise())
end
# white_noise_test()

#problem 3.4
function white_noise_complex()
    return [(1/sqrt(2))*(randn() + im* randn()) for t in 1:1000]
end

function white_noise_complex_test()
    t = [1:1000]
    noise = white_noise_complex() 
    ReNoise = real(noise)
    ImgNoise = imag(noise)
    plot(t, ReNoise, ImgNoise)
end
white_noise_complex_test()

#problem 3.5
function ci_rectangular_test()
    N=1000
    start=-1
    stop=1
    y::Vector{Real} = []
    t=start:(stop/N):stop
    for i in 1:(2*N+1)
        push!(y, CPS.ci_rectangular(t[i]))
    end
    plot(t, y)
end
# ci_rectangular_test()

#problem 3.6
function ci_triangle_test()
    N=1000
    start=-1
    stop=1
    y::Vector{Real} = []
    t=start:(stop/N):stop
    for i in 1:(2*N+1)
        push!(y, CPS.ci_triangle(t[i]))
    end
    plot(t, y)
end
# ci_triangle_test()

#problem 3.7
function ci_literka_M_test()
    N=1000
    start=-1
    stop=1
    y::Vector{Real} = []
    t=start:(stop/N):stop
    for i in 1:(2*N+1)
        push!(y, CPS.ci_literka_M(t[i]))
    end
    plot(t, y)
end
# ci_literka_M_test()

#problem 3.8
function ci_literka_U_test()
    N=1000
    start=-1
    stop=1
    y::Vector{Real} = []
    t=start:(stop/N):stop
    for i in 1:(2*N+1)
        push!(y, CPS.ci_literka_U(t[i]))
    end
    plot(t, y)
end
# ci_literka_U_test()

#problem 3.9
function ramp_test()
    N=1000
    start=-2
    stop=2
    y::Vector{Real} = []
    t=start:(stop/N):stop
    for i in 1:(2*N+1)
        push!(y, CPS.ramp_wave(t[i]))
    end
    plot(t, y)
end
# ramp_test()

#problem 3.10
function sawtooth_test()
    N=1000
    start=-2
    stop=2
    y::Vector{Real} = []
    t=start:(stop/N):stop
    for i in 1:(2*N+1)
        push!(y, CPS.sawtooth_wave(t[i]))
    end
    plot(t, y)
end
# sawtooth_test()

#problem 3.11
function triangular_wave_test()
    N=1000
    start=-2
    stop=2
    y::Vector{Real} = []
    t=start:(stop/N):stop
    for i in 1:(2*N+1)
        push!(y, CPS.triangular_wave(t[i]))
    end
    plot(t, y)
end
# triangular_wave_test()

#problem 3.12
function square_wave_test()
    N=1000
    start=-2
    stop=2
    y::Vector{Real} = []
    t=start:(stop/N):stop
    for i in 1:(N+1)
        push!(y, CPS.square_wave(t[i]))
    end
    plot(t, y)
end
# square_wave_test()

#problem 3.13
function pulse_wave_test()
    N=1000
    stop=1
    y::Vector{Real} = []
    t=0:(stop/N):stop
    for i in 1:(N+1)
        push!(y, CPS.pulse_wave(t[i]))
    end
    plot(t, y)
end
# pulse_wave_test()

#problem 3.14
function impulse_repeater_test()
    t = -10:0.1:10
    f = t -> exp.(t)
    plot(t, CPS.impulse_repeater(f, -2, 3))
end
# impulse_repeater_test()

#problem 3.15
function ramp_bl_test()
    N=1000
    start=-5
    stop=5
    y::Vector{Real} = []
    t=start:(stop/N):stop

    for i in 1:(2*N+1)
        push!(y, CPS.ramp_wave_bl(t[i]))
    end
    plot(t, y)
end
# ramp_bl_test()

function fseries_test()
    N = 1000
    T = 1
    nf = 5

    f(x) = CPS.ramp_wave(x)
    F = CPS.fseries(f, T, N, nf)

    t = -3:0.1:3
    plot(t, CPS.impulse_repeater(F, -2, 3))
end
fseries_test()