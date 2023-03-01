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
using DataFrames
using Markdown


using Plots
theme(:vibrant, size=(500, 350), minorticks=true, grid=false, frame=:box,
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


# πDD: Dˣ⁺
"""
The model with the π⁺D⁰D⁰ system with a single, unsymmetrized Dˣ⁺ resonance.
"""
const model = let
    ch = πDD((m1=mπ⁺, m2=mD⁰, m3=mD⁰), BW(m=mDˣ⁺, Γ=ΓDˣ⁺), ZeroBW(m=mDˣ⁺, Γ=ΓDˣ⁺))
    # 
    set3 = ChannelWithIntegrationMethod(ch, HookSqrtDalitzMapping{3}())
    ich3 = interpolated(set3, cutoff; estep=estep)
    # 
    Amplitude((ich3,))
end


efe_cauchy = let
	effrangepars = # 12s
        effectiverangeexpansion(
	        Δe -> denominator_II(model, Eᵦˣ⁺ + Δe, δm0_val),
	        Δe -> k3b(Eᵦˣ⁺ + Δe),
	        ComplexBranchPointExpansion(CircularSum(abs(imag(Eᵦˣ⁺))/2, 50)))
    # 
	(; tophysicsunits(effrangepars)..., effrangepars...)
end

efe_fit = let
    fitmethod = EffectiveRangeFit(range(-1,1,10), efe2)
    # 
    effrangepars = effectiverangeexpansion(
	    Δe -> denominator_II(model, Δe, δm0_val),
	    Δe -> k3b(Δe),
        fitmethod)
    # 
	(; tophysicsunits(effrangepars)..., effrangepars...)
end


# summarize the resuls
df = DataFrame()
push!(df, (; method = :cauchy, efe_cauchy...))
push!(df, (; method = :fit, NamedTuple{keys(efe_cauchy)}(efe_fit)...))
select(df, :method, :r_fm, :a⁻¹)


# plot
md"""
The approximation using the fit is closer to the function at the range of the fit,
however, once you go to the complex plane, the functions start mismatching.
"""

let
    D(Δe) = denominator_II(model, Δe, δm0_val)
    minusiNk(Δe) = ere(k3b(Δe); a⁻¹=0, r=0, N=efe_cauchy.N)
    R(Δe) = D(Δe)-minusiNk(Δe)
    shift(Δe) = Δe + Eᵦˣ⁺
    # 
    R1(Δe) = ere(k3b(Δe), efe_cauchy)-minusiNk(Δe)
    R2(Δe) = ere(k3b(Δe), efe_fit)-minusiNk(Δe)
    # 
    plot(layout=grid(1,2), size=(700,300), yticks=false, leg=:topleft)
    plot!(sp=1, real ∘ R ∘ shift, -1, 1, title="Real", lab="Model")
    plot!(sp=2, imag ∘ R ∘ shift, -1, 1, title="Imag", lab="")
    #
    plot!(sp=1, real ∘ R1 ∘ shift, -1, 1, c=2, lab="Eff.range(cauchy)")
    plot!(sp=2, imag ∘ R1 ∘ shift, -1, 1, c=2, lab="")
    #
    plot!(sp=1, real ∘ R2 ∘ shift, -1, 1, c=2, ls=:dash, lab="Eff.range(fit)")
    plot!(sp=2, imag ∘ R2 ∘ shift, -1, 1, c=2, ls=:dash, lab="")
end
