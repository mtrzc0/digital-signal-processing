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
            licz+=b[m].*cispi(-2*F[n]*(m-1))
        end
        for k in 2:K
            mian+=a[k].*cispi(-2*F[n]*(k-1))
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
# rozwiazanie()

function rozwiazanie(;
    b::Vector{Float64} = [0.04314525874746985, -0.2588715524848191, 0.6471788812120477, -0.862905174949397, 0.6471788812120477, -0.2588715524848191, 0.04314525874746985],
    a::Vector{Float64} = [1.0, -0.38574424913684624, 1.3013655603245946, 0.13270353683808567, 0.5561248128464324, 0.17251600265356282, 0.160210615122618],
    F::Vector{Float64} = [0.15, 0.16, 0.24],
)
    M=length(b)
    K=length(a)
    N=length(F)
    h=zeros(ComplexF64, N)

    for n in 1:N
        licz=0
        mian=0
        for m in 1:M
           licz+=b[m].*cispi(-2*F[n]*(m-1)) 
        end
        for k in 2:K
            mian+=a[k].*cispi(-2*F[n]*(k-1))
        end
        h[n]=licz/(1+mian)
    end
    return sum(angle.(h))/N
end
# rozwiazanie()

function rozwiazanie(;
    b::Vector{Float64} = [0.009835403572646848, -1.875671518477884e-5, 0.010373051324321612, 0.010373051324321616, -1.8756715184778463e-5, 0.009835403572646846],
    a::Vector{Float64} = [1.0, -3.043987586489922, 3.967544846014289, -2.6908248277859954, 0.9430865718564958, -0.13543960723130033],
    F::Vector{Float64} = [0.01, 0.03, 0.11, 0.29],
)
    M=length(b)
    K=length(a)
    N=length(F)
    h=zeros(ComplexF64, N)

    for n in 1:N
        licz=0
        mian=0
        for m in 1:M
           licz+=b[m]*cispi(-2*F[n]*(m-1)) 
        end
        for k in 2:K
            mian+=a[k]*cispi(-2*F[n]*(k-1))
        end
        h[n]=licz/(1+mian)
    end
    return sum(abs.(h))/N
end
# rozwiazanie()

function rozwiazanie(;
    b::Vector{Float64} = [0.0036995705131426322, -0.006588701086142232, 0.005123154642098176, -2.0536741824576278e-19, -0.005123154642098175, 0.006588701086142232, -0.003699570513142633],
    a::Vector{Float64} = [1.0, -2.978210743260216, 5.773741495194624, -6.608511452910479, 5.494231273099657, -2.69585928272244, 0.8614467267882088],
    F::Vector{Float64} = [0.22, 0.23, 0.24, 0.31, 0.32, 0.36],
)
    M=length(b)
    K=length(a)
    N=length(F)
    h=zeros(ComplexF64, N)

    for n in 1:N
        licz=0
        mian=0
        for m in 1:M
            licz+=b[m]*cispi(-2*F[n]*(m-1))
        end
        for k in 2:K
            mian+=a[k]*cispi(-2*F[n]*(k-1))
        end
        h[n]=licz/(1+mian)
    end
    return sum(abs.(h))/N
end
# rozwiazanie()

function rozwiazanie(;
    b::Vector{Float64} = [0.09690588994826067, -0.4004081385768237, 0.8184157087376562, -1.0223679701595116, 0.8184157087376562, -0.40040813857682367, 0.09690588994826065],
    a::Vector{Float64} = [1.0, -0.6592300629638518, 1.1219500725455562, -0.5095075748876536, 0.28631368108070654, -0.06709412594593915, 0.009731927261286931],
    F::Vector{Float64} = [0.01, 0.08, 0.2, 0.31],
)
    M=length(b)
    K=length(a)
    N=length(F)
    h=zeros(ComplexF64, N)

    for n in 1:N
        licz=0
        mian=0
        for m in 1:M
            licz+=b[m]*cispi(-2*F[n]*(m-1))
        end
        for k in 2:K
            mian+=a[k]*cispi(-2*F[n]*(k-1))
        end
        h[n]=licz/(1+mian)
    end
    return sum(abs.(h))/N
end
# rozwiazanie()

function rozwiazanie(;
    zz::Vector{ComplexF64} = ComplexF64[-1.0 + 0.0im, -1.0 + 0.0im, -1.0 + 0.0im],
    pp::Vector{ComplexF64} = ComplexF64[0.4978837005263824 + 0.6901410648428219im, 0.4978837005263824 - 0.6901410648428219im, 0.662772422800638 + 0.0im],
    k::Float64 = 0.03070522291579536,
    F::Vector{Float64} = [0.08, 0.4, 0.43, 0.44],
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
# rozwiazanie()

function rozwiazanie(;
    b::Vector{Float64} = [6.238698354847942e-5, 0.0, -0.00024954793419391767, 0.0, 0.0003743219012908765, 0.0, -0.00024954793419391767, 0.0, 6.238698354847942e-5],
    a::Vector{Float64} = [1.0, -3.6330270926156896, 8.463618351718203, -12.591890016179114, 14.06407276350723, -11.131225695520087, 6.6138091257372205, -2.5089273284907283, 0.6105348075612239],
    F::Vector{Float64} = [0.02, 0.05, 0.12, 0.26, 0.39],
)
    M=length(b)
    K=length(a)
    N=length(F)
    h=zeros(ComplexF64, N)

    for n in 1:N
        licz=0
        mian=0
        for m in 1:M
            licz+=b[m]*cispi(-2*F[n]*(m-1))
        end
        for k in 2:K
            mian+=a[k]*cispi(-2*F[n]*(k-1))
        end
        h[n]=licz/(mian+1)
    end
    return sum(angle.(h))/N
end
# rozwiazanie()

function rozwiazanie(;
    zz::Vector{ComplexF64} = ComplexF64[1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im, 1.0 + 0.0im],
    pp::Vector{ComplexF64} = ComplexF64[0.44827249582137546 - 0.8556007217889913im, 0.44827249582137546 + 0.8556007217889913im, 0.17291160541979061 - 0.8495953651002048im, 0.17291160541979061 + 0.8495953651002048im, -0.4296758972867222 - 0.5095435898474966im, -0.4296758972867222 + 0.5095435898474966im],
    k::Float64 = 0.03839993839240795,
    F::Vector{Float64} = [0.0, 0.09, 0.14, 0.18, 0.48],
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
# rozwiazanie()