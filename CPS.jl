module CPS

using LinearAlgebra

const author = Dict(
    "index" => "419399",
    "name" => "Mateusz Trzeciak",
    "email" => "trzeciakmat@student.agh.edu.pl",
)

###############################################################################
#  Narzedzia                                                                  #
###############################################################################
function fseries(f::Function, T, N, nf)
    t = range(0, T; length=N)
    dt = T / N
    n = 1:nf

    a0 = sum(f.(t)) * dt / T
    an = N * (2 * sum(f.(t) .* cos.(n' .* (2π / T) .* t) .* dt, dims=1) / T)
    bn = N * (2 * sum(f.(t) .* sin.(n' .* (2π / T) .* t) .* dt, dims=1) / T)

    F(t) = a0 + sum((an .* cos.(n' .* (2π / T) * t) .* dt) + (bn .* sin.(n' .* (2π / T) * t) .* dt))

    return F
end

###############################################################################
# Parametry sygnalow                                                          #
###############################################################################
mean(x::AbstractVector)::Number = missing
function mean(x::AbstractVector)::Number
    return sum(x)/length(x)
end
peak2peak(x::AbstractVector)::Real = missing
function peak2peak(x::AbstractVector)::Real 
    return max(x)-min(x)
end
energy(x::AbstractVector)::Real = missing
function energy(x::AbstractVector)::Real
    return sum(x .* x)
end
power(x::AbstractVector)::Real = missing
function power(x::AbstractVector)::Real 
    return sum(x .* x)/length(x)
end
rms(x::AbstractVector)::Real = missing
function rms(x::AbstractVector)::Real 
    return sqrt(sum(x .* x)/length(x))
end

function running_mean(x::AbstractVector, M::Integer)::Vector
    x_start=M+1
    x_stop=length(x) - M
    index=x_start:x_stop

    running::Vector{Float64} = []
    for i in index
        append!(running, mean(x[(i-M):(i+M)]))
    end

    return running
end

function running_energy(x::AbstractVector, M::Integer)::Vector
    x_start=M+1
    x_stop=length(x) - M
    index=x_start:x_stop

    running::Vector{Float64} = []
    for i in index
        append!(running, energy(x[(i-M):(i+M)]))
    end

    return running
end

function running_power(x::AbstractVector, M::Integer)::Vector
    x_start=M+1
    x_stop=length(x) - M
    index=x_start:x_stop

    running::Vector{Float64} = []
    for i in index
        append!(running, power(x[(i-M):(i+M)]))
    end

    return running
end

###############################################################################
# Modelowanie sygnalow                                                        #
###############################################################################
ci_rectangular(t::Real; T::Real=1.0)::Real = missing
function ci_rectangular(t::Real; T::Real=1.0)::Real
    return (abs(t) != 0.5*T) ? (abs(t) < 0.5*T ? 1 : 0) : 0.5
end

ci_triangle(t::Real; T::Real=1.0)::Real = missing
function ci_triangle(t::Real; T::Real=1.0)::Real
    return (abs(t) <= 1) ? (T-abs(t)) : 0
end

ci_literka_M(t::Real; T=1.0)::Real = missing
function ci_literka_M(t::Real; T=1.0)
    return (abs(t) < 0.5*T) ? 0.2*(abs(10*t)) : 0
end

ci_literka_U(t::Real; T=1.0)::Real = missing
function ci_literka_U(t::Real; T=1.0)
    return (abs(t) < 0.5*T) ? 4*t^2 : 1
end

ramp_wave(t::Real)::Real = missing
function ramp_wave(t::Real)::Real 
    return t>0 ? t%1 : t%1 + 1
end

sawtooth_wave(t::Real)::Real = missing
function sawtooth_wave(t::Real)::Real 
    return t>0 ? 1-t%1 : (1-t)%1
end

triangular_wave(t::Real)::Real = missing
function triangular_wave(t::Real)
    return 4*abs(t - floor(t+ 3/4)+1/4)-1
end

square_wave(t::Real)::Real = missing
function square_wave(t::Real)
    return 1*(2*floor(t)-floor(2*t)+1)
end

pulse_wave(t::Real; D::Real=0.2)::Real = missing
function pulse_wave(t::Real, D::Real=0.2)::Real 
    return D + ((abs(t)< floor(t) + D) ? 1 : 0) 
end

function impulse_repeater(g::Function, t0::Real, t1::Real)::Function 
    return t -> g(mod(t-t0, t1-t0)+t0)
end

function ramp_wave_bl(t; A=1.0, T=1.0, band=20.0)
    N=1000
    f(x) = A .* CPS.ramp_wave(x)
    F = CPS.fseries(f, T, N, band)
    return F(t)
end

function sawtooth_wave_bl(t; A=1.0, T=1.0, band=20.0)
    N=1000
    f(x) = A .* CPS.sawtooth_wave(x)
    F = CPS.fseries(f, T, N, band)
    return F(t)
end

function triangular_wave_bl(t; A=1.0, T=1.0, band=20.0)
    N=1000
    f(x) = A .* CPS.triangular_wave(x)
    F = CPS.fseries(f, T, N, band)
    return F(t)
end

function square_wave_bl(t; A=1.0, T=1.0, band=20.0)
    N=1000
    f(x) = A .* CPS.square_wave(x)
    F = CPS.fseries(f, T, N, band)
    return F(t)
end

function pulse_wave_bl(t; D=0.2, A=1.0, T=1.0, band=20.0)
    N=1000
    f(x) = A .* CPS.pulse_wave(x, D)
    F = CPS.fseries(f, T, N, band)
    return F(t)
end

function impulse_repeater_bl(g::Function, t0::Real, t1::Real, band::Real)::Function
    N=1000
    T=t1-t0
    return CPS.impulse_repeater(fseries(g, T, N, band), t0, t1)
end

#TODO
function rand_siganl_bl(f1::Real, f2::Real)::Function
    B=f2-f1
end

kronecker(n::Integer)::Real = missing
function kronecker(n::Integer)::Real
    return (n==0) ? 1 : 0
end
heaviside(n::Integer)::Real = missing
function heaviside(n::Integer)::Real
    return (n>=1) ? 1 : 0
end

# Dyskretne okna czasowe
rect(N::Integer) = missing
function rect(N::Integer)
    return N>0 ? ones(N) : 0
end

triang(N::Integer) = missing
function triang(N::Integer)
end

hanning(N::Integer) = missing
hamming(N::Integer) = missing
blackman(N::Integer) = missing

function chebwin(N; Îą=-100)
    missing
end

function kaiser(N; Î˛=0, K=20)
    missing
end

###############################################################################
# Próbkowanie i kwantyzacja                                                   #
###############################################################################

quantize(L::AbstractVector)::Function = missing
function quantize(L::AbstractVector)::Function
    return x -> L[argmin(abs.(x .- L))]
end

SQNR(N::Integer)::Real = missing
function SQNR(N::Integer)::Real
    return 20*log10(2^N)
end

SNR(Psignal::Real, Pnoise::Real)::Real = missing
function SNR(Psignal::Real, Pnoise::Real)::Real
    return Psignal/Pnoise
end

function interpolate(
    m::AbstractVector,
    s::AbstractVector,
    kernel::Function = sinc
)::Function
    missing
end

###############################################################################
# Obliczanie dyskretnej transformacji Fouriera                                #
###############################################################################

function dtft(f::Real, signal::Vector, fs::Real)
    N = length(signal)
    w = zeros(ComplexF64,N)
    fd = f/fs
    for i in 0:(N-1)
        w[i+1] = exp(-im*2*pi*fd*i)
    end
    return w.*signal
end

function dft(x::Vector)::Vector
    N = length(x)
    w = zeros(ComplexF64,(N, N))
    for i in 0:(N-1)
        for j in 0:(N-1)
            w[i+1, j+1] = exp(-im*2*pi/N)^(i*j)
        end
    end
    return w*x
end

function idft(x::Vector)::Vector
    N = length(x)
    w = zeros(ComplexF64,(N, N))
    for i in 0:(N-1)
        for j in 0:(N-1)
            w[i+1, j+1] = exp(-im*2*pi/N)^(i*j)
        end
    end
    #bez dzielenia przez N
    return inv(w)*x
end

function goertzel(x::Vector, k::Integer)::Complex
    missing
end

function recursiv_dft(N::Integer)::Function
    missing
end

function exp_recursiv_dft(N::Integer, a::Number)::Function
    missing
end

function cos_recursiv_dft(N::Integer, a::Vector)::Function
    missing
end

function rdft(x::Vector)::Vector
    missing
end

function irdft(X::Vector, N::Integer)::Vector
    missing
end

function fft_radix2_dit_r(x::Vector)::Vector
    missing
end

function ifft_radix2_dif_r(x::Vector)::Vector
    missing
end

function fft(x::Vector)::Vector
    missing
end

function ifft(x::Vector)::Vector
    missing
end

function rfft(x::Vector)::Vector
    missing
end

function irfft(x::Vector, N::Integer)::Vector
    missing
end

###############################################################################
# Analiza częstotliwościowa sygnałów dyskretnych                              #
###############################################################################

fftfreq(N::Integer, fs::Real) = missing
rfftfreq(N::Integer, fs::Real) = missing

amplitude_spectrum(x::Vector, w::Vector = rect(length(x)))::Vector = missing
power_spectrum(x::Vector, w::Vector = rect(length(x)))::Vector = missing
psd(x::Vector, w::Vector = rect(length(x)), fs::Real = 1.0)::Vector = missing

function welch(x::Vector, w::Vector = rect(length(x)), L::Integer = 1, fs::Real = 1.0)::Vector
    missing
end

# Modelowanie sygnałów niestacjonarnych
chirp_lin(t, f0, f1, T, φ = 0) = missing
chirp_exp(t, f0, f1, T, φ = 0) = missing

function stft(x::Vector, w::Vector, L::Integer)::Matrix
    missing
end

function istft(
    X::AbstractMatrix{<:Complex},
    w::AbstractVector{<:Real},
    L::Integer = 0,
    N::Integer = length(w);
)::AbstractVector{<:Real}
    missing
end

"""
Inputs:
    * X - spectrogram
    * w - STFT window
    * L - STFT overlap
    * N - number of Griffin-Lim iterations

Outputs:
    * Y - estimated STFT representation (amplitude + phase)
    * loss_values - a vector with reconstruction loss values (optional, recommended for dev)
"""
function phase_retrieval(X, w, L, N)
    missing, missing
end

###############################################################################
# Systemy dyskretne                                                           #
###############################################################################

function lti_amp(f::Real; b::Vector, a::Vector)::Real
    missing
end

function lti_phase(f::Real; b::Vector, a::Vector)::Real
    missing
end

function conv(f::Vector, g::Vector)::Vector
    missing
end

function fast_conv(f::Vector, g::Vector)::Vector
    missing
end

function overlap_add(f::Vector, g::Vector, L::Integer)::Vector
    missing
end

function overlap_save(f::Vector, g::Vector, L::Integer)::Vector
    missing
end

function lti_filter(b::Vector, a::Vector, x::Vector)::Vector
    missing
end

function filtfilt(b::Vector, a::Vector, x::Vector)::Vector
    missing
end

end