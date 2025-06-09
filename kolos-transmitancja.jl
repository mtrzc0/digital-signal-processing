# transmitancja, zera i bieguny
function rozwiazanie(;
    zz::Vector{ComplexF64} = ComplexF64[0.618237344495477 + 0.7859914667928534im, 0.618237344495477 - 0.7859914667928534im, 0.7880749527481846 + 0.6155792953397204im, 0.7880749527481846 - 0.6155792953397204im, 1.0 + 0.0im],
    pp::Vector{ComplexF64} = ComplexF64[0.4289087685769398 + 0.8645008129514766im, 0.4289087685769398 - 0.8645008129514766im, 0.1483245154470054 + 0.7829549048618939im, 0.1483245154470054 - 0.7829549048618939im, -0.3639175917140422 + 0.0im],
    k::Float64 = 0.14804660838937625,
    F::Vector{Float64} = [0.14, 0.17, 0.22],
)

    N = length(F)
    M = length(zz)
    K = length(pp)
    H = ones(ComplexF64, N)
    for n in 1:NF
        for m in 1:M
            H[n] *= (1-zz[m]/cispi(2*F[n]))
        end
        for k in 1:K
            H[n] /= 1-pp[k]/cispi(2*F[n])
        end
    end
    return sum(angle.(k*H))/N
end
# rozwiazanie()

# transmitancja, zera i bieguny
function rozwiazanie(;
    zz::Vector{ComplexF64} = ComplexF64[-1.0 + 0.0im, -1.0 + 0.0im, -1.0 + 0.0im, -1.0 + 0.0im, -1.0 + 0.0im, -1.0 + 0.0im],
    pp::Vector{ComplexF64} = ComplexF64[0.4829321112470466 + 0.8373793951972854im, 0.4829321112470466 - 0.8373793951972854im, 0.6259309407570826 + 0.6473274226098389im, 0.6259309407570826 - 0.6473274226098389im, 0.8075750782653377 + 0.2579163503294173im, 0.8075750782653377 - 0.2579163503294173im],
    k::Float64 = 0.0006201142564604331,
    F::Vector{Float64} = [0.11, 0.23, 0.44, 0.46, 0.47],
)
    NF = length(F)
    NM = length(zz)
    NK = length(pp)
    H = ones(ComplexF64, NF)

    for n in 1:NF
        for m in 1:NM
            H[n] *= 1-zz[m]/cispi(2*F[n])
        end
        for k in 1:NK
            H[n] /= 1-pp[k]/cispi(2*F[n])
        end
    end
    return sum(angle.(k*H))/NF
end
# rozwiazanie()

# transmitancja, zera i bieguny
function rozwiazanie(;
    zz::Vector{ComplexF64} = ComplexF64[0.6372786378076709 - 0.7706334652699683im, 0.1966960130534983 + 0.9804645217695834im, 0.6372786378076709 + 0.7706334652699683im, 0.1966960130534983 - 0.9804645217695834im, 0.44338445406409877 + 0.8963315379335266im, 0.44338445406409877 - 0.8963315379335266im],
    pp::Vector{ComplexF64} = ComplexF64[-0.36224162461493203 + 0.5983744930699401im, 0.7716708794401513 - 0.33916798101897755im, -0.36224162461493203 - 0.5983744930699401im, 0.7716708794401513 + 0.33916798101897755im, -0.0864628177533698 + 0.0im, 0.5102848542524746 + 0.0im],
    k::Float64 = 0.1517558578114602,
    F::Vector{Float64} = [0.02, 0.04, 0.38, 0.38, 0.4],
)
    M=length(zz)
    K=length(pp)
    N=length(F)

    h = ones(ComplexF64, N)
    for n in 1:N
        for m in 1:M
            h[n] *= 1-zz[m]/cispi(2*F[n])
        end
        for k in 1:K
            h[n] /= 1-pp[k]/cispi(2*F[n])
        end
    end 
    return sum(abs.(k*h)/N)
end
# rozwiazanie()

# transmitancja, sumy
function rozwiazanie(;
    b::Vector{Float64} = [1.7536549719840665e-7, 0.0, -1.05219298319044e-6, 0.0, 2.6304824579761e-6, 0.0, -3.507309943968133e-6, 0.0, 2.6304824579761e-6, 0.0, -1.05219298319044e-6, 0.0, 1.7536549719840665e-7],
    a::Vector{Float64} = [1.0, -5.505694192165391, 18.030017993408308, -40.24823881336225, 68.34799776211682, -90.4722151684874, 96.14190988223578, -81.76012331863264, 55.81795900812245, -29.70323127063752, 12.024170722518782, -3.3177893557296083, 0.5446010675601196],
    F::Vector{Float64} = [0.03, 0.21, 0.21, 0.33, 0.34, 0.47],
)
    M=length(b)
    K=length(a)
    N=length(F)
    h=zeros(ComplexF64, N)

    for n in 1:N
        licz=0
        mian=0
        for m in 1:M
            licz += b[m].*cispi(-2*F[n]*(m-1))
        end
        for k in 2:K
            mian += a[k].*cispi(-2*F[n]*(k-1))
        end
        h[n] = licz / (mian+1)
    end
    return sum(abs.(h))/N
end
# rozwiazanie()

function rozwiazanie(;
    b::Vector{Float64} = [0.6668547023844128, -1.720764471909391, 3.4806618241912046, -3.8658925116288825, 3.4806618241912037, -1.720764471909391, 0.6668547023844127],
    a::Vector{Float64} = [1.0, -2.2575274288819127, 3.948876910415752, -3.813445803637163, 2.951472093696976, -1.2364482229285918, 0.3946840490385092],
    F::Vector{Float64} = [0.14, 0.43, 0.48, 0.5],
)
    M=length(b)
    K=length(a)
    N=length(F)
    h=zeros(ComplexF64, N)

    for n in 1:N
        mian=0
        licz=0
        for m in 1:M
            licz+=b[m].*cispi(-2*F[n]*(m-1))
        end
        for k in 2:K
            mian+=a[k].*cispi(-2*F[n]*(k-1))
        end
        h[n] = licz/(mian +1)
    end
    return sum(abs.(h))/N
end
rozwiazanie()

function rozwiazanie(;
    b::Vector{Float64} = [5.566396272183253e-7, 0.0, -3.339837763309952e-6, 0.0, 8.34959440827488e-6, 0.0, -1.1132792544366507e-5, 0.0, 8.34959440827488e-6, 0.0, -3.339837763309952e-6, 0.0, 5.566396272183253e-7],
    a::Vector{Float64} = [1.0, -6.312470431235184, 22.280704254930217, -53.26748029453974, 95.22115428099832, -131.95980405673234, 145.1934458068807, -127.20440344811901, 88.48086156410804, -47.71128142721598, 19.236501494991767, -5.2531773469685685, 0.802333642696997],
    F::Vector{Float64} = [0.0, 0.06, 0.24, 0.34, 0.35, 0.41],
)
    M=length(b)
    K=length(a)
    N=length(F)
    h=zeros(ComplexF64, N)

    for n in 1:N
        licz=0
        mian=0
        for m in 1:M
            licz+=b[m].*cispi(-2*F[n]*m)
        end
        for k in 2:K
            mian+=a[k].*cispi(-2*F[n]*k)
        end
        h[n]=licz/(mian+1)
    end
    return sum(abs.(h))/N
end
# rozwiazanie()

function rozwiazanie(;
    zz::Vector{ComplexF64} = ComplexF64[1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im, -1.0 + 0.0im, -1.0 + 0.0im, -1.0 + 0.0im, -1.0 + 0.0im, -1.0 + 0.0im, -1.0 + 0.0im],
    pp::Vector{ComplexF64} = ComplexF64[0.543701277994591 - 0.8319925555042785im, 0.29457887606763117 + 0.9483370303390719im, 0.543701277994591 + 0.8319925555042785im, 0.29457887606763117 - 0.9483370303390719im, 0.508848478970115 - 0.8409843938995955im, 0.3274254690915373 + 0.9249844307416222im, 0.508848478970115 + 0.8409843938995955im, 0.3274254690915373 - 0.9249844307416222im, 0.45219457592657175 - 0.864724019694428im, 0.3858890222126597 + 0.8953357230268396im, 0.45219457592657175 + 0.864724019694428im, 0.3858890222126597 - 0.8953357230268396im],
    k::Float64 = 2.9839279959995376e-7,
    F::Vector{Float64} = [0.13, 0.22, 0.43, 0.48],
)
    M=length(zz)
    K=length(pp)
    N=length(F)
    h=ones(ComplexF64, N)

    for n in 1:N
        for m in 1:M
            h[n]*=1-zz[m]*cispi(-2*F[n])
        end
        for k in 1:K
            h[n]/=1-pp[k]*cispi(-2*F[n])
        end
    end
    return sum(angle.(k*h))/N
end
rozwiazanie()

function rozwiazanie(;
    zz::Vector{ComplexF64} = ComplexF64[0.5653138377631579 + 0.8248759087483948im, 0.5653138377631579 - 0.8248759087483948im, 0.9090424796534677 + 0.41670345593176317im, 0.9090424796534677 - 0.41670345593176317im],
    pp::Vector{ComplexF64} = ComplexF64[-0.27610745027671657 + 0.6519220973669482im, -0.27610745027671657 - 0.6519220973669482im, -0.1706544392027319 + 0.2080237764058478im, -0.1706544392027319 - 0.2080237764058478im],
    k::Float64 = 0.05804559733761661,
    F::Vector{Float64} = [0.01, 0.44, 0.46],
)
    M=length(zz)
    K=length(pp)
    N=length(F)
    h=ones(ComplexF64, N)

    for n in 1:N
        for m in 1:M
            h[n]*=1-zz[m]*cispi(-2*F[n])
        end
        for k in 1:K
            h[n]/=1-pp[k]*cispi(-2*F[n])
        end
    end
    return sum(abs.(k*h))/N
end
rozwiazanie()