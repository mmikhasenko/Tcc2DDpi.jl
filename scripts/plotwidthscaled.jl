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



ichs1 = interpolated.(full, cutoff; estep=estep)
ichs2 = ichs1[[1]]
ichs3 = interpolated.(sngl, cutoff; estep=estep) # cutoff
# 
amp1 = Amplitude(Tuple(ichs1), zero)
amp2 = Amplitude(Tuple(ichs2), zero)
amp3 = Amplitude(Tuple(ichs3), zero)
# 

pole1 = pole_position(amp1, δm0_val)
pole2 = pole_position(amp2, δm0_val)
pole3 = pole_position(amp3, δm0_val)

let xmax = 0.5
    plot(xlim=(-0.5,xmax), ylim=(-0.06, 0.0), title="pole position",
        xlab="Re δ√s [MeV]", ylab="Im δ√s [MeV]")
    vline!([δm0_val], lab="δm₀=-369keV", l=(0.6,:gray))
    scatter!(
        [pole1.m_pole pole2.m_pole pole3.m_pole],
        [pole1.half_Γ_pole pole2.half_Γ_pole pole3.half_Γ_pole], mc=[2 3 4], ms=5,
        lab=["(D⁰π⁺)D⁰" "(D⁰π⁺)D⁰+D⁰(π⁺D⁰)" "all channels"])
    plot!([0, xmax], -ΓDˣ⁺*1e3/2 .* [1, 1], lab="", lc=:red, lw=2)
    scatter!([0], [-ΓDˣ⁺*1e3/2], lab="", m=(:red, 6))
end
savefig(joinpath("plots","nominal","polethreemodel.pdf"))
