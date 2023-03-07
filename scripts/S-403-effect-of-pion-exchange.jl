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
    ch1 = πDD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), BW(m=mDˣ⁺,Γ=ΓDˣ⁺), BW(m=mDˣ⁺,Γ=ΓDˣ⁺))
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


function effrange(model)
	effrangepars = # 12s
        effectiverangeexpansion(
	        Δe -> denominator_II(model, Eᵦˣ⁺ + Δe, δm0_val),
	        Δe -> k3b(Eᵦˣ⁺ + Δe),
	        ComplexBranchPointExpansion(CircularSum(abs(imag(Eᵦˣ⁺))/2, 50)))
    # 
	(; tophysicsunits(effrangepars)..., effrangepars...)
end

efe2 = effrange(model2)
efe3 = effrange(model3)

efe2.r_fm
efe3.r_fm

efe2.r_fm
efe3.r_fm

Δr_ope = efe2.r_fm - efe3.r_fm
Δr_ope = efe2.r_fm - efe3.r_fm

# nominal
# -0.9396341248075353 + 0.28446733393925405im
# 
# 0.95% width
# -0.8999020429222457 + 0.2583926613878467im
# 
# mass + 10e-6
# -0.9374367692088035 + 0.28323707879081084im

ere(k, p) = p.N*(p.a⁻¹ + p.r * k^2 / 2 - 1im* k)

real(efe2.r)+1im*imag(efe2.r)

let
    plot()
    plot!(Δe->real(denominator_II(model2, Δe, δm0_val)), -1, 1, lab="real")
    plot!(Δe->imag(denominator_II(model2, Δe, δm0_val)), -1, 1, lab="imag")
    #
    plot!(Δe->real(ere(k3b(Δe),efe2)), -1, 1, c=1, ls=:dash, lab="Eff.range")
    plot!(Δe->imag(ere(k3b(Δe),efe2)), -1, 1, c=2, ls=:dash, lab="")
    # 
    # efe2′ = (; efe2..., r = real(efe3.r)+1im*imag(efe3.r))
    # plot!(Δe->real(ere(k3b(Δe),efe2′)), -1, 1, c=1, ls=:dot, lab="r → -0.77fm")
    # plot!(Δe->imag(ere(k3b(Δe),efe2′)), -1, 1, c=2, ls=:dot, lab="")
end


fitmethod = EffectiveRangeFit(range(-0.1,0.1,10), efe2)


efe3_fit = let
    effrangepars = effectiverangeexpansion(
	    Δe -> denominator_II(model3, Eᵦˣ⁺ + Δe, δm0_val),
	    Δe -> k3b(Eᵦˣ⁺ + Δe),
        fitmethod)
    # 
	(; tophysicsunits(effrangepars)..., effrangepars...)
end

efe3
efe3_fit


efe3_fit.χ0
efe3_fit.minimum