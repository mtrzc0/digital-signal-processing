module CPS

const author = Dict(
    "index" => "419399",
    "name" => "Mateusz Trzeciak",
    "email" => "trzeciakmat@student.agh.edu.pl",
)



###############################################################################
# Parametry sygnalow                                                          #
###############################################################################
mean(x::AbstractVector)::Number = missing
peak2peak(x::AbstractVector)::Real = missing
energy(x::AbstractVector)::Real = missing
power(x::AbstractVector)::Real = missing
rms(x::AbstractVector)::Real = missing

function running_mean(x::AbstractVector, M::Integer)::Vector
    missing
end

function running_energy(x::AbstractVector, M::Integer)::Vector
    missing
end

function running_power(x::AbstractVector, M::Integer)::Vector
     missing
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
    return (abs(t) < 0.5*T) ? 1 - abs(t) : 0
end

ci_literka_U(t::Real; T=1.0)::Real = missing
function ci_literka_U(t::Real; T=1.0)
    return (abs(t) < 0.5*T) ? 4*t^2 : 0
end

ramp_wave(t::Real)::Real = missing
function ramp_wave(t::Real)::Real 
    return abs(2*(1-(t%1)))
end

sawtooth_wave(t::Real)::Real = missing
function sawtooth_wave(t::Real)::Real 
    return abs(2*(t%1))
end

triangular_wave(t::Real)::Real = missing
function triangular_wave(t::Real)
    return 4*(4*abs(t - floor(t+ 3/4)+1/4)-1)
end

square_wave(t::Real)::Real = missing
function square_wave(t::Real)
    return (4*floor(t)-2*floor(2*t)+1)
end
# pulse_wave(t::Real; D::Real=0.2)::Real = missing

function impuse_repeater(g::Function, t0::Real, t1::Real)::Function 
    missing
end

function ramp_wave_bl(t; A=1.0, T=1.0, band=20.0)
    missing
end

function sawtooth_wave_bl(t; A=1.0, T=1.0, band=20.0)
    missing
end

function triangular_wave_bl(t; A=1.0, T=1.0, band=20.0)
    missing
end

function square_wave_bl(t; A=1.0, T=1.0, band=20.0)
    missing
end

# function pulse_wave_bl(t; Ď=0.2, A=1.0, T=1.0, band=20.0)
#     missing
# end


function impuse_repeater_bl(g::Function, t0::Real, t1::Real, band::Real)::Function
    missing
end

function rand_siganl_bl(f1::Real, f2::Real)::Function
    missing
end



kronecker(n::Integer)::Real = missing
heaviside(n::Integer)::Real = missing

# Dyskretne okna czasowe
rect(N::Integer) = missing
triang(N::Integer) = missing
hanning(N::Integer) = missing
hamming(N::Integer) = missing
blackman(N::Integer) = missing

function chebwin(N; Îą=-100)
    missing
end

function kaiser(N; Î˛=0, K=20)
    missing
end

end