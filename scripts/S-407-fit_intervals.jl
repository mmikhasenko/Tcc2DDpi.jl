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
import Plots.PlotMeasures: mm
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


# πDD: Dˣ⁺
"""
The model with the π⁺D⁰D⁰ system with a single, unsymmetrized Dˣ⁺ resonance.
"""
const model = let
    ch1 = πDD((m1=mπ⁺, m2=mD⁰, m3=mD⁰), BW(m=mDˣ⁺, Γ=ΓDˣ⁺), BW(m=mDˣ⁺, Γ=ΓDˣ⁺))
    iπDD2 = interpolated(
        ChannelWithIntegrationMethod(ch1, HookSqrtDalitzMapping{2}()),
        cutoff; estep=estep)
    iπDD3 = interpolated(
        ChannelWithIntegrationMethod(ch1, HookSqrtDalitzMapping{3}()),
        cutoff; estep=estep)
    Amplitude((iπDD2, iπDD3))
end


samples = [
    range(-0.1, 0.1, 10),
    range(-1, 1, 10),
    range(0.5, 1, 10),
    range(-0.5, -0.2, 10)]
# 

# efe_fit = let
function fitefe(sample)
    # 
    efe_ini = (a⁻¹=-0.027 + 0im, r=-3.3 + 1.7, N=0.007 + 0im)
    fitmethod = EffectiveRangeFit(sample, efe_ini)
    # 
    effrangepars = effectiverangeexpansion(
        Δe -> denominator_II(model, Δe, δm0_val),
        Δe -> k3b(Δe),
        fitmethod)
    # 
    (; tophysicsunits(effrangepars)..., effrangepars...)
end


efes = fitefe.(samples)


# summarize the resuls
df = DataFrame(efes)
df.minE = getindex.(samples, 1)
df.maxE = last.(samples)
select(df, :minE, :maxE,
    [:r_fm, :a_fm] .=> ByRow(x -> round(x; digits=2)) .=> [:r_fm, :a_fm])


# plot
md"""
The approximation using the fit is closer to the function at the range of the fit,
however, once you go to the complex plane, the functions start mismatching.
"""


let
    plot(layout=grid(1, 2), size=(700, 320),
        leg=:topright, link=:x, ylim=(-5e-6, 9e-6), yticks=false)

    D(Δe) = denominator_II(model, Δe, δm0_val)
    R(Δe) = D(Δe)
    # plot!(sp=1, real ∘ R, -1.1, 1.1, lw=2, title="Real", lab="Model")
    # plot!(sp=2, imag ∘ R, -1.1, 1.1, lw=2, title="Imag", lab="")
    plot!(sp=1, [-1.1, 1.1], [0, 0], lw=2, title="Re[Ere-M] ", lab="", legtitle="r [fm]")
    plot!(sp=2, [-1.1, 1.1], [0, 0], lw=2, title="Im[Ere-M] ", lab="")

    map(enumerate(eachrow(df))) do (c, fit)
        R2(Δe) = ere(k3b(Δe), fit) - D(Δe)
        # 
        plot!(sp=1, real ∘ R2, -1.1, 1.1, lab="", c=c + 1, lw=2)
        plot!(sp=2, imag ∘ R2, -1.1, 1.1, lab="", c=c + 1, lw=2)
        # 
        shift = c != 2 ? 0 : -1e-6
        plot!(sp=1, [fit.minE, fit.maxE], [0, 0] .+ shift, lw=12, c=c + 1, lab="$(round(fit.r_fm; digits=2))", alpha=0.7)
        plot!(sp=2, [fit.minE, fit.maxE], [0, 0] .+ shift, lw=12, c=c + 1, lab="", alpha=0.7)
        # 
    end
    vline!(sp=1, [δm0_val], lc=:red, lab="", ann=(δm0_val, 3e-6, text("Tcc", rotation=90, :bottom, :red, 10)))
    vline!(sp=2, [δm0_val], lc=:red, lab="", ann=(δm0_val, 3e-6, text("Tcc", rotation=90, :bottom, :red, 10)))
    plot!()
end
savefig(joinpath("plots", "effective-range-fitintervals.pdf"))

