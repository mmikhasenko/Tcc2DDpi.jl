using Pkg
cd(joinpath(@__DIR__,".."))
Pkg.activate(".")
# 
using X2DDpi
using Parameters
using QuadGK

using Plots
import Plots.PlotMeasures: mm
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



p1 = let xmax = 0.5
    plot(xlim=(-0.5,xmax), ylim=(-0.06, 0.0),
        xlab="Re δ√s [MeV]", ylab="Im δ√s [MeV]")
    vline!([δm0_val], lab="δm₀=-369keV", l=(0.6,:gray))

    for (i,δm) in enumerate([-0.42,δm0_val,-0.25,-0.15,-0.1, -0.05, -0.03, -0.02, -0.01, -0.001])
        poles = pole_position.(amps, δm)
        scatter!(
            getproperty.(poles, :m_pole)',
            getproperty.(poles, :half_Γ_pole)', mc=[2 3 4], ms=5,
            lab=(i==1 ? permutedims(labels) : ""))
    end
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
    plot(p1),
    layout=grid(2,1,heights=(0.3,0.7)), size=(500,600), link=:x)
savefig(joinpath("plots","nominal","polethreemodel.pdf"))

#

using FiniteDiff

ρ_thr(full[1],δm0_val)


derivative_data = let
    xv = range(-0.5, 2.5, length=300)
    yvs = [map(
        e->FiniteDiff.finite_difference_derivative(ex->sum(ρ_thr.(m,ex)), e),
        xv) for m in models]
    norm_yvs = yvs ./ [sum(ρ_thr.(m,δm0_val)) for m in models]
    (; xv, yvs=norm_yvs)
end

function relacemax(yv)
    _,i = findmax(yv)
    yv′ = copy(yv)
    yv′[i] = (yv′[i-1]+yv′[i+1])/2
    return yv′
end
function relacemin(yv)
    _,i = findmin(yv)
    yv′ = copy(yv)
    yv′[i] = (yv′[i-1]+yv′[i+1])/2
    return yv′
end

corr_derivative_data = let 
    @unpack xv, yvs = derivative_data
    yvs = [yvs[1],relacemin(relacemax(derivative_data.yvs[2])),relacemin(relacemax(derivative_data.yvs[3]))]
    (; xv, yvs)
end

let
    @unpack xv, yvs = corr_derivative_data
    plot(xlab="e [MeV]", ylab="ρ'(e) / ρ(δm₀)")
    plot!(xv, [yvs[1] yvs[2] yvs[3]],
        lab=permutedims(labels), c=[2 3 4])
end
savefig(joinpath("plots","nominal","rhoprimeoverrho.pdf"))


integrals_m1_p2h = 
    [sum(y ./
        (corr_derivative_data.xv .- δm0_val) .*
        (corr_derivative_data.xv .> 0.0)) for y in corr_derivative_data.yvs]
#

writejson(joinpath("results","nominal","widththreemodels.json"), Dict(
        :integral_m1_p2h => integrals_m1_p2h
        )
    )

derivative_data_middle = let
    xv = range(2.5, 10.5, length=30)
    yvs = [map(
        e->FiniteDiff.finite_difference_derivative(ex->sum(ρ_thr.(m,ex)), e),
        xv) for m in models]
    norm_yvs = yvs ./ [sum(ρ_thr.(m,δm0_val)) for m in models]
    (; xv, yvs=norm_yvs)
end
