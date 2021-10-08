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

sngl = [DˣD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), (m=mDˣ⁺,Γ=ΓDˣ⁺))]
# 
full = [
    πDD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), (m=mDˣ⁺,Γ=ΓDˣ⁺), (m=mDˣ⁺,Γ=ΓDˣ⁺)),
    πDD((m1=mπ⁰,m2=mD⁺,m3=mD⁰), (m=mDˣ⁺,Γ=ΓDˣ⁺), (m=mDˣ⁰,Γ=ΓDˣ⁰)),
    γDD((m1=mγ, m2=mD⁺,m3=mD⁰), (m=mDˣ⁺,Γ=ΓDˣ⁺), (m=mDˣ⁰,Γ=ΓDˣ⁰))]
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

