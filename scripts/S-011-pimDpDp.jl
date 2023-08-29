using Pkg
cd(joinpath(@__DIR__,".."))
Pkg.activate(".")
# 
using X2DDpi
using Parameters
using Measurements
using Interpolations
using Statistics
using Plots

theme(:wong2, frame=:box, grid=false, minorticks=true, 
    guidefontvalign=:top, guidefonthalign=:right,
    xlim=(:auto,:auto))


settings = transformdictrecursively!(readjson("settings.json"), ifstringgivemeasurement)
@unpack cutoff, estep = settings["phspmatching"]
#
channels = [
    πDD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), (m=mDˣ⁺,Γ=ΓDˣ⁺), (m=mDˣ⁺,Γ=ΓDˣ⁺)),
    πDD((m1=mπ⁰,m2=mD⁺,m3=mD⁰), (m=mDˣ⁺,Γ=ΓDˣ⁺), (m=mDˣ⁰,Γ=ΓDˣ⁰)),
    γDD((m1=mγ, m2=mD⁺,m3=mD⁰), (m=mDˣ⁺,Γ=ΓDˣ⁺), (m=mDˣ⁰,Γ=ΓDˣ⁰)),
    πDD((m1=mπ⁺,m2=mD⁺,m3=mD⁺), (m=mDˣ⁰,Γ=ΓDˣ⁰), (m=mDˣ⁰,Γ=ΓDˣ⁰))]
#
ichannels13 = interpolated.(channels[1:3], cutoff; estep=estep)
ichannels4  = interpolated( channels[4], 50.0; estep=0.1)
# 
ichannels = [ichannels13..., ichannels4]
ampX0 = Amplitude(Tuple(ichannels), zero)
# ampX1 = Amplitude(Tuple(ichannels), zero)

plot( e->ρ_thr( channels[4],e), m2e(sum(channels[4].ms)):1.0:55, yscale=:log, ylim=(1e-8,:auto))
plot!(e->ρ_thr(ichannels[4],e), m2e(sum(channels[4].ms)):1.0:55, yscale=:log, leg=:left)


@unpack δm0 = settings["fitresults"]
δm0_val = δm0.val
#

caldata = let
    ev = -1:0.003:5.5
    calv = map(e->1/abs2(denominator_I(ampX0, e, δm0_val)), ev)
    (; xv = ev, yv = calv)
end
# 
phsps = let 
    ev = caldata.xv
    yvs = [map(e->ρ_thr(ich,e), ev) for ich in ichannels]
    (;xv = ev, yvs, labs=["π⁺D⁰D⁰","π⁰D⁺D⁰","γD⁺D⁰","10³ × π⁻D⁺D⁺"])
end

let
    plot(xlab="e [MeV]", ylab="|A|²ρᵢ")
    plot!(phsps.xv, caldata.yv .* phsps.yvs[1], lab=phsps.labs[1], yscale=:log10)
    plot!(phsps.xv, caldata.yv .* phsps.yvs[2], lab=phsps.labs[2], yscale=:log10)
    plot!(phsps.xv, caldata.yv .* phsps.yvs[3], lab=phsps.labs[3], yscale=:log10)
    # 
    m = phsps.xv .> m2e(sum(channels[4].ms))
    plot!(phsps.xv[m], 1e3 .* caldata.yv[m] .* phsps.yvs[4][m], lab=phsps.labs[4], yscale=:log10)
    # 
    vline!([0 m2e(mDˣ⁰+mD⁺) m2e(mπ⁺+2mD⁺)], lc=[:lightgreen :magenta], ls=:dash, lw=1.5, lab="")
    annotate!([
        (0, 2, text("Dˣ⁺D⁰",10, :left, :bottom, rotation=90)),
        (m2e(mDˣ⁰+mD⁺), 2, text("Dˣ⁰D⁺",10, :left, :bottom, rotation=90)),
        (m2e(mπ⁺+2mD⁺), 2, text("π⁺D⁻D⁺",10, :left, :bottom, rotation=90))
    ])
    plot!(ylim=(1e0,:auto))
end
savefig(joinpath("plots","nominal","intensityspectra.pdf"))

