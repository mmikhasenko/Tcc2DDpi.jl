md"""
# Effective range of the OPE model

Two close models are compared: the `model2``has the OPE (Dˣ⁺ + Dˣ⁺), while the `model3` does not.
The effective range flips the sign. The calculation show that the EFE does not work very well for the OPE model.
"""

using Pkg
cd(joinpath(@__DIR__, ".."))
Pkg.activate(".")
# 
using X2DDpi
using Parameters
using Measurements
using Interpolations
using Statistics
using LaTeXStrings

using Plots
theme(:wong2, size=(500, 350), minorticks=true, grid=false, frame=:box,
    guidefontvalign=:top, guidefonthalign=:right,
    foreground_color_legend=nothing,
    legendfontsize=9, legend=:topright,
    xlim=(:auto, :auto), ylim=(:auto, :auto))


#        _|              _|                
#    _|_|_|    _|_|_|  _|_|_|_|    _|_|_|  
#  _|    _|  _|    _|    _|      _|    _|  
#  _|    _|  _|    _|    _|      _|    _|  
#    _|_|_|    _|_|_|      _|_|    _|_|_|  

settings = transformdictrecursively!(readjson("settings.json"), ifstringgivemeasurement)
#
@unpack δm0 = settings["fitresults"]
const δm0_val = δm0.val
@unpack cutoff, estep = settings["phspmatching"]


#                            _|            
#    _|_|_|    _|_|      _|_|_|    _|_|    
#  _|        _|    _|  _|    _|  _|_|_|_|  
#  _|        _|    _|  _|    _|  _|        
#    _|_|_|    _|_|      _|_|_|    _|_|_|  


# πDD: Dˣ⁺ + Dˣ⁺
"""
The model with the π⁺D⁰D⁰ system with two Dˣ⁺ resonances.
"""
const model2 = let
    ch1 = πDD((m1=mπ⁺, m2=mD⁰, m3=mD⁰), BW(m=mDˣ⁺, Γ=ΓDˣ⁺), BW(m=mDˣ⁺, Γ=ΓDˣ⁺))
    iπDD2 = interpolated(
        ChannelWithIntegrationMethod(ch1, HookSqrtDalitzMapping{2}()),
        cutoff; estep=estep)
    iπDD3 = interpolated(
        ChannelWithIntegrationMethod(ch1, HookSqrtDalitzMapping{3}()),
        cutoff; estep=estep)
    Amplitude((iπDD2, iπDD3))
end

# πDD: Dˣ⁺
"""
The model with the π⁺D⁰D⁰ system with a single, unsymmetrized Dˣ⁺ resonance.
"""
const model3 = let
    ch = πDD((m1=mπ⁺, m2=mD⁰, m3=mD⁰), BW(m=mDˣ⁺, Γ=ΓDˣ⁺), ZeroBW(m=mDˣ⁺, Γ=ΓDˣ⁺))
    # 
    set3 = ChannelWithIntegrationMethod(ch, HookSqrtDalitzMapping{3}())
    ich3 = interpolated(set3, cutoff; estep=estep)
    # 
    Amplitude((ich3,))
end

function effrange(model; method)
    effrangepars = # 12s
        effectiverangeexpansion(
            Δe -> denominator_II(model, Eᵦˣ⁺ + Δe, δm0_val),
            Δe -> k3b(Eᵦˣ⁺ + Δe),
            method)
    # 
    (; tophysicsunits(effrangepars)..., effrangepars...)
end

md"""
## Expansion in the branch point

The expansion coeffiecients are obtained using the cauchy integral around the branch point.
"""

cauchysum = ComplexBranchPointExpansion(CircularSum(abs(imag(Eᵦˣ⁺)) / 2, 50))

efe2 = effrange(model2; method=cauchysum)
efe3 = effrange(model3; method=cauchysum)

efe2.r_fm
efe3.r_fm

Δr_ope = efe2.r_fm - efe3.r_fm

ere(k, p) = p.N * (p.a⁻¹ + p.r * k^2 / 2 - 1im * k)

real(efe2.r) + 1im * imag(efe2.r)

let (model, efe) = (model2, efe2)

    invD(Δe) = denominator_II(model, Δe, δm0_val)
    expansion(Δe) = ere(k3b(Δe), efe)
    minus_iNk(Δe) = ere(k3b(Δe); a⁻¹=0, r=0, N=efe.N)
    # 
    plot(layout=grid(1, 2), size=(800, 350), title=["Real" "Imag"])
    # 
    elim = (-1, 1)
    plot!(sp=1, Δe -> real(invD(Δe) - minus_iNk(Δe)), elim..., lab="1/D")
    plot!(sp=1, Δe -> real(expansion(Δe) - minus_iNk(Δe)), elim..., c=1, ls=:dash, lab="Expansion")
    # 
    plot!(sp=2, Δe -> imag(invD(Δe) - minus_iNk(Δe)), elim..., lab="1/D")
    plot!(sp=2, Δe -> imag(expansion(Δe) - minus_iNk(Δe)), elim..., c=1, ls=:dash, lab="Expansion")
end
savefig(joinpath("plots", "ope_effect_effrange_expansion.pdf"))


md"""
## Fit on the real axis

The expansion coeffiecients are obtained by fitting the inverse amplitude on the real axis around the expansion point.
"""

fitmethod = EffectiveRangeFit(range(-0.1, 0.1, 30),
    (a⁻¹=-0.027 + 0.0014, r=-3.3 + 1.6, N=0.007 + 0.0im))
efe2_fit = effrange(model2; method=fitmethod)

let (model, efe) = (model2, efe2_fit)

    invD(Δe) = denominator_II(model, Δe, δm0_val)
    expansion(Δe) = ere(k3b(Δe), efe)
    minus_iNk(Δe) = ere(k3b(Δe); a⁻¹=0, r=0, N=efe.N)
    # 
    plot(layout=grid(1, 2), size=(800, 350), title=["Real" "Imag"])
    # 
    elim = (-1, 1)
    plot!(sp=1, Δe -> real(invD(Δe) - minus_iNk(Δe)), elim..., lab="1/D")
    plot!(sp=1, Δe -> real(expansion(Δe) - minus_iNk(Δe)), elim..., c=1, ls=:dash, lab="Fit")
    # 
    plot!(sp=2, Δe -> imag(invD(Δe) - minus_iNk(Δe)), elim..., lab="1/D")
    plot!(sp=2, Δe -> imag(expansion(Δe) - minus_iNk(Δe)), elim..., c=1, ls=:dash, lab="Fit")
end
savefig(joinpath("plots", "ope_effect_effrange_fit.pdf"))

