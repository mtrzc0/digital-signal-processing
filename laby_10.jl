include("CPS.jl")

using Plots
using Random
using LinearAlgebra

#problem 10.1
function firwin_lp_I_test()
    M=20
    fir = CPS.firwin_lp_I(M, 0.4)
    plot(fir,
     seriestype = :stem,
     marker     = (:circle, 5, :blue),
     line       = (:solid, 1, :blue))
end
# firwin_lp_I_test()

#problem 10.2
function firwin_hp_I_test()
    m=20
    fir = cps.firwin_hp_I(m, 0.4)
    plot(fir,
     seriestype = :stem,
     marker     = (:circle, 5, :black),
     line       = (:solid, 1, :black))
end
# firwin_hp_I_test()

function firwin_bp_I_test()
    M=20
    fir = CPS.firwin_bp_I(M, 0.4, 0.2)
    plot(fir,
     seriestype = :stem,
     marker     = (:circle, 5, :black),
     line       = (:solid, 1, :black))
end
# firwin_bp_I_test()

function firwin_bs_I_test()
    M=20
    fir = CPS.firwin_bs_I(M, 0.4, 0.2)
    plot(fir,
     seriestype = :stem,
     marker     = (:circle, 5, :black),
     line       = (:solid, 1, :black))
end
firwin_bs_I_test()