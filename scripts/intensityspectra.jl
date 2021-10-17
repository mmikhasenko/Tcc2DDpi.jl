using Pkg
cd(joinpath(@__DIR__,".."))
Pkg.activate(".")
# 
using X2DDpi
using Parameters
using Measurements
using Interpolations
using Statistics

settings = transformdictrecursively!(readjson("settings.json"), ifstringgivemeasurement)
@unpack cutoff, estep = settings["phspmatching"]
#
channels = [
    πDD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), (m=mDˣ⁺,Γ=ΓDˣ⁺), (m=mDˣ⁺,Γ=ΓDˣ⁺)),
    πDD((m1=mπ⁰,m2=mD⁺,m3=mD⁰), (m=mDˣ⁺,Γ=ΓDˣ⁺), (m=mDˣ⁰,Γ=ΓDˣ⁰)),
    γDD((m1=mγ, m2=mD⁺,m3=mD⁰), (m=mDˣ⁺,Γ=ΓDˣ⁺), (m=mDˣ⁰,Γ=ΓDˣ⁰))]
#
ichannels = interpolated.(channels, cutoff; estep=estep) # cutoff

ampX0 = Amplitude(Tuple(ichannels), zero)

@unpack δm0 = settings["fitresults"]
δm0_val = δm0.val
#

caldata = let
    ev = -1:0.003:3
    calv = map(e->1/abs2(denominator_I(ampX0, e, δm0_val)), ev)
    (; xv = ev, yv = calv)
end
# 
phsps = let 
    ev = caldata.xv
    yvs = [map(e->ρ_thr(ich,e), ev) for ich in ichannels]
    (;xv = ev, yvs, labs=["π⁺D⁰D⁰","π⁰D⁺D⁰","γD⁺D⁰"])
end

let
    plot(xlab="e [MeV]", ylab="|A|²ρᵢ")
    plot!(phsps.xv, caldata.yv .* phsps.yvs[1], lab=phsps.labs[1], yscale=:log10)
    plot!(phsps.xv, caldata.yv .* phsps.yvs[2], lab=phsps.labs[2], yscale=:log10)
    plot!(phsps.xv, caldata.yv .* phsps.yvs[3], lab=phsps.labs[3], yscale=:log10)
    vline!([0 m2e(mDˣ⁰+mD⁺)], lc=[:lightgreen :magenta], ls=:dash, lw=1.5, lab="")
end
savefig(joinpath("plots","nominal","intensityspectra.pdf"))

