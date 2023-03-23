using Pkg
cd(joinpath(@__DIR__, ".."))
Pkg.activate(".")
Pkg.instantiate()
# 
using X2DDpi
using Parameters
using Plots

theme(:vibrant, size=(500, 350), minorticks=true, grid=false, frame=:box,
    guidefontvalign=:top, guidefonthalign=:right,
    foreground_color_legend=nothing,
    legendfontsize=9, legend=:topright, lab="",
    xlim=(:auto, :auto), ylim=(:auto, :auto))
# 

settings = transformdictrecursively!(readjson("settings.json"),
    ifstringgivemeasurement)
@unpack δm0 = settings["fitresults"]
const δm0_val = δm0.val
@unpack cutoff, estep = settings["phspmatching"]
#

const model = let
    channels = [
        πDD((m1=mπ⁺, m2=mD⁰, m3=mD⁰), BW(m=mDˣ⁺, Γ=ΓDˣ⁺), ZeroBW(m=mDˣ⁺, Γ=ΓDˣ⁺)),
        πDD((m1=mπ⁰, m2=mD⁺, m3=mD⁰), ZeroBW(m=mDˣ⁺, Γ=ΓDˣ⁺), BW(m=mDˣ⁰, Γ=ΓDˣ⁰)),
        γDD((m1=mγ, m2=mD⁺, m3=mD⁰), ZeroBW(m=mDˣ⁺, Γ=ΓDˣ⁺), BW(m=mDˣ⁰, Γ=ΓDˣ⁰))]
    # 
    @time ichannels = interpolated.(channels, cutoff; estep=estep) # cutoff
    # 
    Amplitude(Tuple(ichannels), zero)
end

let xlim = m2e(mDˣ⁰ + mD⁺) .+ (-0.5, 0.5)
    plot(Δe -> ρ_thr(model.ik[3], Δe), xlim...)
end


function k3b_high(e)
    m = e2m(e)
    M = sqrt(mDˣ⁰^2 - 1im * mDˣ⁰ * ΓDˣ⁰) # taken at the Dˣ⁰ pole
    p = cis(π / 4) * sqrt((m - (M + mD⁺)) * cis(-π / 2)) *
        sqrt(m + (M + mD⁺)) * sqrt(m - (M - mD⁺)) * sqrt(m + (M - mD⁺)) / (2 * m)  # branch cut down
    return p
end


# denominator_II(model, Δe, δm0_val)

ere_guess = (N=3e-2 * cis(π / 2 - 0.9), r=-22.1, a⁻¹=-0.03)

efe_high_fit = let
    fitmethod = EffectiveRangeFit(m2e(mDˣ⁰ + mD⁺) .+ range(-0.1, 0.1, 50), ere_guess)
    # 
    effrangepars = effectiverangeexpansion(
        Δe -> denominator_II(model, Δe, δm0_val),
        Δe -> k3b_high(Δe),
        fitmethod)
    # 
    (; tophysicsunits(effrangepars)..., effrangepars...)
end

let xlim = m2e(mDˣ⁰ + mD⁺) .+ (-0.5, 0.5)
    # 
    D(Δe) = denominator_II(model, Δe, δm0_val)
    minusiNk(Δe) = ere(k3b_high(Δe); a⁻¹=0, r=0, N=efe_high_fit.N)
    # 
    R(Δe) = D(Δe) - minusiNk(Δe)
    shift(Δe) = Δe# + Eᵦˣ⁺
    # 
    R1(Δe) = ere(k3b_high(Δe), efe_high_fit) - minusiNk(Δe)

    plot(layout=grid(1, 2), size=(700, 300), yticks=false, leg=:topleft)
    # 
    # plot(R.(range(xlim..., 100)), lab="Model")
    # plot!(R1.(range(xlim..., 100)), lab="Eff.range(fit)")

    # plot!(sp=1, real ∘ minusiNk, xlim..., lab="iNk")
    # plot!(sp=2, imag ∘ minusiNk, xlim..., lab="iNk")
    plot!(sp=1, real ∘ R ∘ shift, xlim..., title="Real", lab="Model: D + ikN")
    plot!(sp=2, imag ∘ R ∘ shift, xlim..., title="Imag", lab="")
    #
    plot!(sp=1, real ∘ R1 ∘ shift, xlim..., c=2, lab="Eff.range(fit): N(1/a+rk^2)")
    plot!(sp=2, imag ∘ R1 ∘ shift, xlim..., c=2, lab="")
    #
    vline!(sp=1, [m2e(mDˣ⁰ + mD⁺)], lab="", l=(:dash, :black))
    vline!(sp=2, [m2e(mDˣ⁰ + mD⁺)], lab="", l=(:dash, :black))
end

