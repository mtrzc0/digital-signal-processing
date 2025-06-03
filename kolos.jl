# Z ciągłego sygnału f(t)∈Rf(t)∈R o paśmie ograniczonym od dołu i góry przez częstotliwość ∣B∣<12Δm∣B∣<2Δm1​, zostało pobrane 8080 próbek w równych odstępach czasu Δm=mn+1−mnΔm=mn+1​−mn​. Wartości sygnału oraz momenty w których zostały pobrane kolejne próbki, znajdują się odpowiednio w wektorze s∈R80s∈R80 oraz w wektorze m∈R80m∈R80, gdzie sn=f(mn)sn​=f(mn​). Na podstawie wektorów mm oraz ss, znajdź sygnał g(t)g(t), będący rekonstrukcją sygnału f(t)f(t) otrzymaną z wykorzystaniem wzoru interpolacyjnego Whittakera-Shannona. Jako rozwiązanie podaj sumę wartości sygnału g(t)g(t) dla momentów t∈R5t∈R5, to znaczy:

function rozwiazanie(;
    m::Vector{Float64} = [-2.1, -2.098, -2.096, -2.094, -2.092, -2.09, -2.088, -2.086, -2.084, -2.082, -2.08, -2.078, -2.076, -2.074, -2.072, -2.07, -2.068, -2.066, -2.064, -2.062, -2.06, -2.058, -2.056, -2.054, -2.052, -2.05, -2.048, -2.046, -2.044, -2.042, -2.04, -2.038, -2.036, -2.034, -2.032, -2.03, -2.028, -2.026, -2.024, -2.022, -2.02, -2.018, -2.016, -2.014, -2.012, -2.01, -2.008, -2.006, -2.004, -2.002, -2.0, -1.998, -1.996, -1.994, -1.992, -1.99, -1.988, -1.986, -1.984, -1.982, -1.98, -1.978, -1.976, -1.974, -1.972, -1.97, -1.968, -1.966, -1.964, -1.962, -1.96, -1.958, -1.956, -1.954, -1.952, -1.95, -1.948, -1.946, -1.944, -1.942],
    s::Vector{Float64} = [0.0551, 0.7162, 0.8468, 0.049, 0.8809, 0.0812, 0.1316, 0.638, 0.4675, 0.1182, 0.099, 0.0764, 0.1608, 0.9381, 0.4316, 0.9282, 0.4092, 0.5099, 0.6477, 0.4124, 0.5819, 0.2118, 0.1954, 0.0332, 0.0359, 0.7057, 0.2342, 0.9939, 0.5717, 0.009, 0.0632, 0.9706, 0.4979, 0.1286, 0.4441, 0.0787, 0.2207, 0.0379, 0.4855, 0.8344, 0.9141, 0.6686, 0.6008, 0.2538, 0.1748, 0.8121, 0.4465, 0.5145, 0.0834, 0.1011, 0.7905, 0.8337, 0.2869, 0.7453, 0.0794, 0.9527, 0.8243, 0.9327, 0.0539, 0.3925, 0.6367, 0.097, 0.2839, 0.0693, 0.2656, 0.0359, 0.7723, 0.71, 0.2734, 0.4897, 0.9531, 0.7275, 0.0807, 0.8415, 0.834, 0.5258, 0.7134, 0.5195, 0.4541, 0.0816],
    t::Vector{Float64} = [-1.9586, -1.9766, -1.9956, -1.9952, -2.0518],
)

    t_out = zeros(Float64, length(t))
    T = m[2] - m[1]
    for i in eachindex(t_out)
        for n in eachindex(s)
            t_out .+= sinc((t[i]-m[n])/T) .* s[n]
        end
    end
    return sum(t_out)
end
# rozwiazanie()

# Dany jest dyskretny system liniowy niezmienny w czasie, który jest opisany poniższym równaniem różnicowym.
# ∑m=0Mbmx[n−m]=∑k=0Kaky[n−k]
# m=0∑M​bm​x[n−m]=k=0∑K​ak​y[n−k] Współczynniki tego systemu znajdują się odpowiednio w wektorach b=[b0,b1,…,bM]∈RM+1b=[b0​,b1​,…,bM​]∈RM+1 i a=[a0,a1,…,aK]∈RK+1a=[a0​,a1​,…,aK​]∈RK+1. System ten został pobudzony impulsowym sygnałem x[n]∈Rx[n]∈R, którego niezerowe próbki znajdują się w wektorze x=[x[0],x[1],…,x[16]]T∈R17x=[x[0],x[1],…,x[16]]T∈R17. Znajdź sygnał y[n]∈Ry[n]∈R będący odpowiedą systemu na pobudzenie go sygnałem x[n]x[n]. Jako odpowiedź podaj moc sygnału utworzonego z pierwszych L=42L=42 próbek znalezionego sygnału y[n]y[n], to znaczy dla n=0,1,…,41n=0,1,…,41
function rozwiazanie(;
    b::Vector{Float64} = [0.17071186813902622, -0.6828474725561049, 1.0242712088341572, -0.6828474725561049, 0.17071186813902622],
    a::Vector{Float64} = [1.0, -0.9302156732075635, 1.0304520379086297, -0.22591463885332275, 0.2520337923516871],
    x::Vector{Float64} = [-0.77, -0.39, -0.81, 0.25, 0.78, -0.63, 0.24, 0.44, -0.05, 0.9, -0.58, -0.43, -0.1, 0.05, 0.05, -0.94, 0.42],
    L::Int = 42,
)
    N = length(x)
    M = length(b)
    K = length(a)
    y = zeros(Float64,L)
    for n in 1:L
        for m in 1:M
            if 0 < n-m+1 <= N
                y[n]+=b[m]*x[n-m+1]
            end
        end
        for k in 2:K
            if 0 < n-k+1 <= L
                y[n]-=a[k]*y[n-k+1]
            end
        end
    end

    return sum(y.*y)/L
end
rozwiazanie()

# Dany jest dyskretny system liniowy niezmienny w czasie, który jest opisany poniższym równaniem różnicowym.
# ∑m=0Mbmx[n−m]=∑k=0Kaky[n−k]
# m=0∑M​bm​x[n−m]=k=0∑K​ak​y[n−k] Współczynniki tego systemu znajdują się odpowiednio w wektorach b=[b0,b1,…,bM]∈RM+1b=[b0​,b1​,…,bM​]∈RM+1 i a=[a0,a1,…,aK]∈RK+1a=[a0​,a1​,…,aK​]∈RK+1. System ten został pobudzony impulsowym sygnałem x[n]∈Rx[n]∈R, którego niezerowe próbki znajdują się w wektorze x=[x[0],x[1],…,x[15]]T∈R16x=[x[0],x[1],…,x[15]]T∈R16. Znajdź sygnał y[n]∈Ry[n]∈R będący odpowiedą systemu na pobudzenie go sygnałem x[n]x[n]. Jako odpowiedź podaj wartość skuteczną sygnału utworzonego z pierwszych L=32L=32 próbek znalezionego sygnału y[n]y[n], to znaczy dla n=0,1,…,31n=0,1,…,31. 
function rozwiazanie(;
    b::Vector{Float64} = [0.47939668565494736, -2.396718455777035, 6.390951633328491, -10.884046460065196, 12.9616414851719, -10.884046460065196, 6.390951633328492, -2.3967184557770356, 0.4793966856549476],
    a::Vector{Float64} = [1.0, -4.1984253859082425, 9.378322690837937, -13.567509522160536, 13.900025177246398, -10.198938105424634, 5.334737255706614, -1.8376536321486605, 0.34743102338159215],
    x::Vector{Float64} = [0.74, -0.18, 0.19, 0.03, 0.79, 0.61, 0.32, -0.3, 0.11, -0.7, -0.94, 0.23, -0.95, 0.41, -0.26, -0.6],
    L::Int = 32,
)

    N = length(x)
    M = length(b)
    K = length(a)
    y = zeros(Float64, L)

    for n in 1:L
        for m in 1:M
            if 0 < n-m+1 <= N
                y[n] += b[m]*x[n-m+1]
            end
        end
        for k in 2:K
            if 0 < n-k+1 <= L
                y[n] -= a[k]*y[n-k+1]
            end
        end
    end

    return sqrt(sum(y .* y)/length(y))
end
# rozwiazanie()
          
# takie jak nizej xd
function rozwiazanie(;
    fp::Float64 = 363.32,
    t1::Float64 = -6.95,
    N::Int = 76,
)
    g(t) = -2*(t-floor(t+1/2))
    x = range(start=t1, step=1/fp, length=N)
    y = [2.4*g(2.4*t-2.7) for t in x]
    return sum(y)/length(y)
end
# rozwiazanie()

# Oblicz średnią dyskretnego sygnału x∈R198x∈R198. Dyskretny sygnał xx powstał w wyniki pobrania N=198N=198 próbek z ciągłego sygnału y(t)=1.7⋅g(4.8⋅t−2.4)y(t)=1.7⋅g(4.8⋅t−2.4) z szybkością fp=355.84fp​=355.84 próbek na sekundę. Pierwsza próbka x1=y(t1)x1​=y(t1​) została pobrana w chwili t1=2.97t1​=2.97. Funkcja g(t)g(t) zwraca wartości sygnału fali piłokształtnej o narastającym zboczu i następujących parametrach: amplituda 11, okres 11 sekunda, składowa stała 00, g(0)=0g(0)=0, oraz dgdt∣t=0=2dtdg​∣t=0​=2. 
function rozwiazanie(;
    fp::Float64 = 355.84,
    t1::Float64 = 2.97,
    N::Int = 198,
)
    g(t) = 2*(t-floor(t+1/2)) #piloksztaltny
    x = range(start=t1, step=1/fp, length=N)
    y = [1.7*g(4.8*t-2.4) for t in x]
    return sum(y)/length(y)
end
# rozwiazanie()
        
function rozwiazanie(;
    order::Int = 30,
    fp::Float64 = 192.0,
    f0::Float64 = 49.92,
    z::Vector{Int} = [1, 10, 1, 19, 21],
)
    F0 = f0/fp
    M = order

    δ(n) = (n==0) ? 1 : 0
    w = [0.5 + 0.5cospi(2n/(2M+1)) for n in -M÷2:M÷2]
    h = [(n!=0) ? δ(n) - 2F0*sinc(2*F0*n) : 1-2F0 for n in -M÷2:M÷2]
    return sum(h .* w)
end
# rozwiazanie()

# Dany jest idealny równomierny 99-bitowy kwantyzator q(x)q(x), którego najmniejszy poziom kwantyzacji ma wartość a=0.0017a=0.0017, natomist największy poziom kwantyzacji ma wartość b=1.0b=1.0. Oblicz sygnał błęd kwantyzacji tego przetwornika dla dyskretnego sygnału x∈R58x∈R58. Jako rozwiązanie podaj wartość skuteczną sygnału błędu kwantyzacji.
function rozwiazanie(;
    a::Float64 = 0.0017,
    b::Float64 = 1.0,
    x::Vector{Float64} = [0.32519, 0.66977, 0.22112, 0.86303, 0.84786, 0.78474, 0.74299, 0.89449, 0.40521, 0.67294, 0.12336, 0.8489, 0.89497, 0.992, 0.37029, 0.93429, 0.23715, 0.44712, 0.55996, 0.91066, 0.77214, 0.72759, 0.24131, 0.96841, 0.92991, 0.97906, 0.69307, 0.18012, 0.01652, 0.42954, 0.63751, 0.01805, 0.80662, 0.93477, 0.84842, 0.3272, 0.08417, 0.99283, 0.17723, 0.65494, 0.42265, 0.52204, 0.34743, 0.0591, 0.70169, 0.50724, 0.89066, 0.52285, 0.27484, 0.74983, 0.0017, 0.02743, 0.98037, 0.99545, 0.71153, 0.76461, 0.54003, 0.33978],
)
    N=9
    L=range(start=a,stop=b,length=2^N)
    quant(L)=x->L[argmin(abs.(x .- L))]
    q = quant(L)
    xq = q.(x)
    e = x - xq
    return sqrt(sum(e .* e)/length(e))
end
# rozwiazanie()

# Dany jest dyskretny system liniowy niezmienny w czasie, który jest opisany Z-transmitacją H(z)H(z). Transmitancja H(z)H(z) jest zdefiniowana poprzez z=[z1,z2,…,zM]∈CMz=[z1​,z2​,…,zM​]∈CM, p=[p0,p1,…,pN]∈CNp=[p0​,p1​,…,pN​]∈CN, k∈Rk∈R które są odpowiednio zerami transmitancji, biegnami transmitancji oraz współczynnikiem wzmocnienica systemu. Zbadaj czy system jest stabilny i jako odpowiedź podaj 1.01.0 jeżeli system jest stabilny lub −1.0−1.0 jeżeli system jest niestabilny.       
function rozwiazanie(;
    z::Vector{ComplexF64} = ComplexF64[0.7747988212229449 - 0.6322078666321191im, 0.6338425965825865 + 0.7734620629077061im, 0.7747988212229449 + 0.6322078666321191im, 0.6338425965825865 - 0.7734620629077061im, 0.7666641181763219 - 0.6420483859499398im, 0.6456791161713742 + 0.7636088520572251im, 0.7666641181763219 + 0.6420483859499398im, 0.6456791161713742 - 0.7636088520572251im, 0.7373439644119779 - 0.675517489148303im, 0.6833668826051585 + 0.7300751356939281im, 0.7373439644119779 + 0.675517489148303im, 0.6833668826051585 - 0.7300751356939281im],
    p::Vector{ComplexF64} = ComplexF64[0.6234220702867519 + 0.7791638624837406im, 0.7794189603309418 - 0.6237826096508589im, 0.6234220702867519 - 0.7791638624837406im, 3.1229915436017106 + 2.4993846879966686im, 0.6031422630334754 + 0.7797891758613019im, 0.7801747257138382 - 0.6077733503452839im, 0.6031422630334754 - 0.7797891758613019im, 0.7801747257138382 + 0.6077733503452839im, 0.4726111372588799 + 0.7293157695771523im, 0.7564132806349837 - 0.5074855000695155im, 0.4726111372588799 - 0.7293157695771523im, 0.7564132806349837 + 0.5074855000695155im],
    k::Float64 = 0.5890882678896954,
)
    for pole in p
        if abs(pole) > 1
            return -1
        end
    end
    for pole in p
        if abs(pole) == 1
            return 0
        end
    end
    return 1
end
# rozwiazanie()

# Dany jest dyskretny sygnał x∈C33x∈C33, którego próbki zostały probrane z ciągłego sygnału o ograniczonym paśmie, z szybkością fp=693fp​=693 próbek na sekunde. Oblicz 33-punktową dyskretną transformację Fouriera tego sygnału oraz znajdź wartości dyskretnego widma Fouriera tego sygnału dla następujących częstotliwość f=[−189,−168,−42,105,252,273]f=[−189,−168,−42,105,252,273]. Jako rozwiązanie podaj sumę faz tych składowych częstotliwościowych.         
function rozwiazanie(;
    fp::Int = 693,
    x::Vector{ComplexF64} = ComplexF64[-0.65 + 0.13im, -1.02 + 0.82im, -0.9 + 0.88im, -1.49 + 0.52im, 0.68 + 0.37im, -0.16 - 0.52im, 0.37 + 1.18im, 0.0 - 0.77im, 0.31 + 0.09im, -0.25 - 0.38im, -0.6 - 0.21im, -0.66 - 0.23im, -0.25 - 0.1im, 0.56 + 0.25im, -0.07 + 0.04im, -1.35 - 0.18im, 0.02 - 0.8im, 0.56 + 0.22im, 0.82 + 1.47im, 0.63 + 1.05im, -0.61 + 0.47im, -0.8 - 0.17im, -0.5 - 0.87im, 0.88 + 0.2im, -0.35 + 0.09im, -0.79 + 0.02im, 0.48 + 0.51im, 0.17 - 0.48im, 0.2 - 1.1im, 0.26 - 0.9im, 0.19 + 1.27im, -0.0 + 0.78im, -1.2 - 0.36im],
    f::Vector{Int} = [-189, -168, -42, 105, 252, 273],
) 
    N = length(x)
    dft(x) = [sum(x[n]*cispi(-2(n-1)*(k-1)/length(x)) for n in eachindex(x)) for k in eachindex(x)]/length(x)
    delta_f = fp /N
    out = 0
    F = zeros(Float64, length(x))
    X = dft(x)

    F[Int(floor(N/2)+2)] = -(floor((N-1)/2))*delta_f

    for i in 1:N-1
        if i+1 != floor(N/2)+2
            F[i+1] = F[i]+delta_f
        end
    end

    for i in 1:N
        if F[i] in f
            out += abs(X[i])
        end
    end
    println(out)
    return out
end
rozwiazanie()

# Oblicz odpowiedź impulsową h∈R51h∈R51 nierekursywnego filtru górnoprzepustowego rzędu 50 o liniowej charakterystyce fazowej. Filtr zaprojektuj tak aby przy częstotliwości próbkowania fp=192.0fp​=192.0 Hz, 3 dB pasmo przepustowe zaczynało się na częstotliwość f0=48.0f0​=48.0 Hz. Do zaprojektowania filtru wykorzystaj metodę okien czasowych i okno Hamminga. Jako rozwiązanie podaj sumę wartości wektora hh o indeksach z=[15,32,10]z=[15,32,10], to znaczy,
# ∑i=13hzi.
# i=1∑3​hzi​​.
          
function rozwiazanie(;
    order::Int = 50,
    fp::Float64 = 192.0,
    f0::Float64 = 48.0,
    z::Vector{Int} = [15, 32, 10],
)
    M = order
    F0 = f0/fp
    δ(n) = n==1 ? 1 : 0
    w = [0.54 + 0.46*cospi(2n/2(2M+1)) for n in -M÷2:M÷2]
    h = [(n!=0) ? δ(n)-2F0*sinc(2*F0*n) : 1-2F0 for n in -M÷2:M÷2]

    return sum((h .* w)[i] for i in z)
end
# rozwiazanie()