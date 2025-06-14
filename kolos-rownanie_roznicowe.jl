# rr, moc
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
# rozwiazanie()

# rr, rms
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

# rr, moc
function rozwiazanie(;
    b::Vector{Float64} = [0.6627634385563629, -2.386969976902741, 4.853880730576752, -5.920666906870027, 4.853880730576753, -2.3869699769027424, 0.6627634385563634],
    a::Vector{Float64} = [1.0, -3.1143057871887554, 5.4495020741233535, -5.784175869137851, 4.144565536196449, -1.7961252043489033, 0.4392207279464298],
    x::Vector{Float64} = [-0.89, 0.94, -0.23, -0.53, -0.07, 0.02, 0.09, -0.58, 0.73, 0.67, -0.05, -0.64, -0.98, -0.57, 0.85, -0.29, 0.37, -0.97, 0.25, 0.02, -0.64, -0.93],
    L::Int = 46,
)
    M = length(b)
    K = length(a)
    N = length(x)

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
    return sum(y .* y)/length(y)
end
# rozwiazanie()

# rr, moc
function rozwiazanie(;
    b::Vector{Float64} = [0.3086830209856637, -2.1783661412235524, 8.257373012318105, -20.936679297054166, 39.112116500166, -56.08688586466961, 63.144026107497716, -56.08688586466961, 39.112116500166, -20.936679297054166, 8.257373012318103, -2.1783661412235515, 0.30868302098566364],
    a::Vector{Float64} = [1.0, -5.704304078570587, 17.360262034122854, -35.503548986311266, 53.76506187599037, -62.88191324895661, 58.065950717583576, -42.544104993519596, 24.595822264801132, -10.968678204031226, 3.6179890060178987, -0.8013130945053464, 0.09528527592140434],
    x::Vector{Float64} = [-0.71, 0.16, 0.42, -0.66, 0.03, 0.29, 0.34, -0.13, -0.69, -0.54, 0.31, 0.18, 0.74, 0.12, 0.05, -0.57, 0.71, 0.17, -0.67, 0.27, -0.96, -0.52, 0.65, -0.38, -0.46, -0.43, -0.96, 0.03, -0.31, -0.56, 0.96, -0.16, 0.52, 0.19, -0.91, 0.12, -0.35],
    L::Int = 48,
)       
    M = length(b)
    K = length(a)
    N = length(x)
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
    return sum(y .* y)/length(y)
end
# rozwiazanie()

# rr, rms
function rozwiazanie(;
    b::Vector{Float64} = [0.17327242644714178, -0.8417443295265945, 1.8658198947340179, -2.3919295938041114, 1.8658198947340177, -0.8417443295265944, 0.17327242644714175],
    a::Vector{Float64} = [1.0, -1.9234073946110977, 2.3846685845801554, -1.7457173609326966, 0.836213821540958, -0.23354278972658468, 0.03005294382812647],
    x::Vector{Float64} = [-0.82, 0.05, -0.99, -0.94, 0.95, -0.05, -0.19, 0.54, 0.92, -0.18],
    L::Int = 33,
)
    M = length(b)
    K = length(a)
    N = length(x)
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

function rozwiazanie(;
    b::Vector{Float64} = [0.15220776723622292, -0.6306612471456423, 1.1561685097078895, -1.1561685097078895, 0.6306612471456425, -0.15220776723622292],
    a::Vector{Float64} = [1.0, -1.1397725742429585, 1.51371415095299, -0.3171939570028648, 0.16876931161974204, 0.2613749456390448],
    x::Vector{Float64} = [-0.12, -0.57, -0.18, -0.44, -0.71, -0.46, -0.38, 0.9, 0.43, 0.28, -0.41, 0.73, -0.7, -0.24, -0.63, -0.06, -0.85, 0.07, 0.76, -0.37, 0.99, 0.21, 0.26, -0.16, -0.01, -0.66, 0.48],
    L::Int = 54,
)
    M=length(b)
    K=length(a)
    N=length(x)
    y=zeros(Float64, L)

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
    return sum(y .* y)
end
# rozwiazanie()

function rozwiazanie(;
    b::Vector{Float64} = [0.3757899352158341, -1.9527551364270286, 5.285845756657121, -9.094927720306911, 10.852128285073531, -9.094927720306911, 5.28584575665712, -1.9527551364270275, 0.3757899352158339],
    a::Vector{Float64} = [1.0, -4.281959903460196, 9.568903988559057, -13.930574496897968, 14.43659413591283, -10.845331345268585, 5.887617239568579, -2.1526678868169586, 0.43046923518197955],
    x::Vector{Float64} = [-0.46, 0.52, 0.37, 0.28, -0.18, 0.64, 0.96, -0.58, -0.91, 0.49, 0.96, -0.24, 0.18, 0.28, 0.89, 0.7, 0.92, -0.85, 0.6, -0.94, 0.36, 0.48, -0.41, -0.87, -0.2, 0.32, 0.48, 0.12, 0.86, -0.42, -0.29, -0.52, -0.26, -0.39, -0.7, 0.61, -0.85, -0.84, 0.16, -0.88, 0.04, -0.59, -0.03, -0.42, -0.35, -0.11, -0.86, 0.44],
    L::Int = 91,
)
    M=length(b)
    K=length(a)
    N=length(x)
    y=zeros(Float64, L)

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
    return sqrt(sum(y.*y)/L)
end
# rozwiazanie()

function rozwiazanie(;
    b::Vector{Float64} = [0.0001824576217927976, 0.0, -0.0007298304871711904, 0.0, 0.0010947457307567856, 0.0, -0.0007298304871711904, 0.0, 0.0001824576217927976],
    a::Vector{Float64} = [1.0, -4.3462771138252085, 10.723355668678993, -17.10139259989802, 19.577251856850896, -15.985988924574716, 9.369620185372915, -3.5484690027943397, 0.7635949896497858],
    x::Vector{Float64} = [-0.85, -0.9, 0.15, 0.24, 0.54, -0.55, 0.35, 0.23, 0.58, 0.26, 0.72, 0.54, 0.25, 0.18, 0.65, -0.27, 0.96, -0.05, 0.26, 0.08, 0.75, -0.66, 0.68, 0.6, 0.82, -0.78, 0.85, -0.05, 0.25, -0.98, -0.82, 0.97, 0.08, -0.56, 0.2, -0.62, 0.13, 0.96, 0.28, 0.79, -0.92, -0.28, 0.13, -0.96, -0.11, -0.79, 0.24, 0.59, 0.29],
    L::Int = 68,
)
    M=length(b)
    K=length(a)
    N=length(x)
    y=zeros(Float64, L)

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
    return sum(y.*y)
end
# rozwiazanie()

function rozwiazanie(;
    b::Vector{Float64} = [0.005886216155083775, 0.0, -0.017658648465251326, 0.0, 0.017658648465251326, 0.0, -0.005886216155083775],
    a::Vector{Float64} = [1.0, -2.8390104753828656, 4.8984221995853146, -5.0962271628207905, 3.7271905661816733, -1.637353111316897, 0.43922072794642975],
    x::Vector{Float64} = [0.65, 0.72, 0.8, 0.58, -0.31, -0.92, -0.9, 0.21, 0.87, -0.52, -0.51, -0.37, 0.66, -0.29, -0.09, 0.44, -0.94, -0.25, -0.95, 0.15, -0.74, 0.28, -0.95, -0.99, 0.33, 0.42, -0.89, 0.89, -0.03],
    L::Int = 50,
)
    M=length(b)
    K=length(a)
    N=length(x)
    y=zeros(Float64, L)

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
# rozwiazanie()

function rozwiazanie(;
    b::Vector{Float64} = [0.5697113391879268, -3.250657215695663, 10.267592693373901, -21.46893144083625, 32.78492171829188, -37.53906825761523, 32.78492171829189, -21.46893144083625, 10.267592693373901, -3.2506572156956626, 0.5697113391879266],
    a::Vector{Float64} = [1.0, -5.068608469554383, 14.19202148337911, -26.37389188361491, 35.86216926073992, -36.65858175312293, 28.641101749714185, -16.81841424031455, 7.2245880063854955, -2.058749224072284, 0.32457100148871776],
    x::Vector{Float64} = [0.61, -0.87, -0.47, -0.42, -0.09, 0.87, 0.95, 0.75, 0.63, -0.14, -0.61, 0.24, 0.85, 0.8, 0.32, 0.14, 0.35, 0.99, 0.63],
    L::Int = 55,
)
    M=length(b)
    K=length(a)
    N=length(x)
    y=zeros(Float64, L)

    for n in 1:L
        for m in 1:M
            if 0 < n-m+1 <=N
                y[n]+=b[m]*x[n-m+1]
            end
        end
        for k in 2:K
            if 0 < n-k+1 <=L
                y[n]-=a[k]*y[n-k+1]
            end
        end
    end
    return sqrt(sum(y.*y)/L)
end
# rozwiazanie()

function rozwiazanie(;
    b::Vector{Float64} = [0.4502459803463013, -3.3159134256014844, 12.800325198551652, -32.85570786053933, 61.78676481661911, -88.951336129986, 100.24244837863881, -88.95133612998602, 61.78676481661912, -32.85570786053934, 12.800325198551654, -3.3159134256014844, 0.4502459803463013],
    a::Vector{Float64} = [1.0, -6.493802318349136, 22.079413867010267, -50.10480710110084, 83.67169483588289, -107.51680505930148, 108.76881149517737, -87.2330360058761, 55.210830496153505, -27.023446447346085, 9.814678874415318, -2.408637623492223, 0.3149989138594487],
    x::Vector{Float64} = [0.44, -0.98, -0.13, 0.4, 0.92, 0.29, -0.05, 0.47, 1.0, 0.55, -0.9, -0.12, 0.55, -0.46, -0.3, -0.95, 0.18, 0.29, -0.64, 0.4, 0.87, 0.18, 0.21, -0.31, 0.64, 0.12, -0.78, -0.36, -0.87, 0.19, 0.94, -0.42, -0.03, 0.84, 0.12, -0.95, 0.58, 0.22, -0.01, 0.48, 0.88, -0.71, 0.46, -0.72, -0.78, -0.21, -0.04],
    L::Int = 72,
)
    M=length(b)
    K=length(a)
    N=length(x)
    y=zeros(Float64, L)

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
# rozwiazanie()