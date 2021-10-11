using Pkg
cd(joinpath(@__DIR__,".."))
Pkg.activate(".")
# 
using X2DDpi
using Parameters
using QuadGK

using Plots
theme(:wong2, frame=:box, grid=false, minorticks=true, 
    guidefontvalign=:top, guidefonthalign=:right,
    xlim=(:auto,:auto), ylim=(:auto,:auto))


settings = transformdictrecursively!(readjson("settings.json"), ifstringgivemeasurement)
@unpack cutoff, estep = settings["phspmatching"]
@unpack δm0 = settings["fitresults"]
const δm0_val = δm0.val
#
# retrieve the model
modelDict = readjson(joinpath("results","nominal","model.json"))
ichs1 = interpolated.(d2nt.(modelDict["ichannels"]))
# full = [
#     πDD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), (m=mDˣ⁺,Γ=ΓDˣ⁺), (m=mDˣ⁺,Γ=ΓDˣ⁺)),
#     πDD((m1=mπ⁰,m2=mD⁺,m3=mD⁰), (m=mDˣ⁺,Γ=ΓDˣ⁺), (m=mDˣ⁰,Γ=ΓDˣ⁰)),
#     γDD((m1=mγ, m2=mD⁺,m3=mD⁰), (m=mDˣ⁺,Γ=ΓDˣ⁺), (m=mDˣ⁰,Γ=ΓDˣ⁰))]
full = getproperty.(ichs1, :channel)

# 
sngl = [DˣD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), (m=mDˣ⁺,Γ=ΓDˣ⁺))]

# 
models = [sngl, full[1:1], full]
labels = ["(D⁰π⁺)D⁰", "(D⁰π⁺)D⁰+D⁰(π⁺D⁰)", "all channels"]

plt = let
    plot()
    for (f,l) in zip(models,labels)
        n = sum(ρ_thr.(f, 0.0))
        plot!(e->sum(ρ_thr.(f,e)) / n * (ΓDˣ⁺*1e6), -2:0.02:1, lab=l)
    end
    vline!([δm0_val], lab="", lw=2, leg=:topright)
    scatter!([0], [ΓDˣ⁺*1e6], lab="", m=(2, :red))
    scatter!([δm0_val δm0_val δm0_val],
        [sum(ρ_thr.(f, δm0_val)) / sum(ρ_thr.(f, 0.0)) * (ΓDˣ⁺*1e6) for f in models]',
        c=[1 2 3], lab="", xlab="δm [MeV]", ylab="ρ_thr / norm")
    lens!(-0.360 .+ 0.05 .* [-1,1], [0,30], inset = (1, bbox(0.1, 0.1, 0.5, 0.5)))
end
savefig(joinpath("plots","nominal","rhothrscaled.pdf"))

# plt = let
#     plot()
#     for (f,l) in zip([sngl, full[1:1], full],
#                      ["(D⁰π⁺)D⁰", "(D⁰π⁺)D⁰+D⁰(π⁺D⁰)", "all channels"])
#         n = 1
#         plot!(e->sum(ρ_thr.(f,e)) / n * (ΓDˣ⁺*1e6), -2:0.02:1, lab=l)
#     end
#     plot!()
# end



ichvs = [
    ichs1,
    ichs1[[1]],
    interpolated.(sngl, cutoff; estep=estep)]
# 
amps = [Amplitude(Tuple(ichv), zero) for ichv in ichvs]
# 

poles = pole_position.(amps, δm0_val)


p1 = let xmax = 0.5
    plot(xlim=(-0.5,xmax), ylim=(-0.06, 0.0),
        xlab="Re δ√s [MeV]", ylab="Im δ√s [MeV]")
    vline!([δm0_val], lab="δm₀=-369keV", l=(0.6,:gray))
    scatter!(
        getproperty.(poles, :m_pole)',
        getproperty.(poles, :half_Γ_pole)', mc=[2 3 4], ms=5,
        lab=permutedims(labels))
    plot!([0, xmax], -ΓDˣ⁺*1e3/2 .* [1, 1], lab="", lc=:red, lw=2)
    scatter!([0], [-ΓDˣ⁺*1e3/2], lab="", m=(:red, 6),
        top_margin=-1mm)
end

p2 = let
    plot()
    xv = range(-0.5, 0.5, length=300)
    yvs = [map(e->1/abs2(denominator_I(a, e, δm0_val)), xv) for a in amps]
    norm_yvs = yvs ./ maximum.(yvs)
    plot!(xv, [norm_yvs[1] norm_yvs[2] norm_yvs[3]],
        lab=permutedims(labels), c=[2 3 4],
        xtickfont=(:white,), bottom_margin=-1mm)
end

plot(
    p2,
    plot(p1, bottom_margin=8mm),
    layout=grid(2,1,heights=(0.3,0.7)), size=(500,600), link=:x)
savefig(joinpath("plots","nominal","polethreemodel.pdf"))

import Plots.PlotMeasures: mm
plot(rand(10), xtickfont=(:white,), bottom_margin=0mm)