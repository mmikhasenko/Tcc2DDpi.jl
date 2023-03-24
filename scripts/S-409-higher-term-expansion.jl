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
The model with the π⁺D⁰D⁰ system with two Dˣ⁺ resonances.
"""
const model = let
    ch1 = πDD((m1=mπ⁺, m2=mD⁰, m3=mD⁰), BW(m=mDˣ⁺, Γ=ΓDˣ⁺), BW(m=mDˣ⁺, Γ=ΓDˣ⁺))
    # 
    @time iπDD2 = interpolated(
        ChannelWithIntegrationMethod(ch1, HookSqrtDalitzMapping{2}()),
        cutoff; estep=estep)
    @time iπDD3 = interpolated(
        ChannelWithIntegrationMethod(ch1, HookSqrtDalitzMapping{3}()),
        cutoff; estep=estep)
    # 
    Amplitude((iπDD2, iπDD3))
end

# ht_cauchy = let
function ht(factor, Npoints)
    effrangepars = # 12s
        highertermexpansion(
            Δe -> denominator_II(model, Eᵦˣ⁺ + Δe, δm0_val),
            Δe -> k3b(Eᵦˣ⁺ + Δe),
            ComplexBranchPointExpansion(CircularSum(abs(imag(Eᵦˣ⁺)) * factor, Npoints)))
    # 
    (; tophysicsunits(effrangepars)..., effrangepars...)
end

factors = [0.01, 0.1, 0.5]
Npoints = [50, 100]
fN = Tuple.(Pair.(factors', Npoints))
htes = map((f, N) -> ht(f, N), fN)

df = DataFrame(vcat(htes...))
df.factor = vcat(getindex.(fN, 1)...)
df.Npoints = vcat(getindex.(fN, 2)...)
df_summary = select(df, :factor, :Npoints,
    :N => ByRow(x -> round(x * 1e3; digits=5)) => :Nx1e3,
    [:a_fm, :r_fm] .=> ByRow(x -> round(x; digits=2)) .=> [:a_fm, :r_fm],
    :ξ, :ζ)

writejson(joinpath("results", "nominal", "higher-term-table.json"), df_summary)


let xlim = (-0.5, 0.5)
    plot(layout=grid(1, 2), size=(700, 300), yticks=false, leg=:topright,
        xlab="E (MeV)", bottom_margin=3mm)
    plot!(sp=1, legtitle="r [fm]")
    D(Δe) = denominator_II(model, Δe, δm0_val)
    shift(Δe) = Δe + Eᵦˣ⁺

    df_sel = df[df.Npoints.==100, :]
    for i in 1:size(df_sel, 1)
        dfi = df_sel[i, :]
        ht_cauchy = (; dfi...)
        # 
        minus_iNk(Δe) = hte(k3b(Δe); N=ht_cauchy.N)
        R1(Δe) = hte(k3b(Δe); NamedTuple{(:a⁻¹, :r, :N)}(ht_cauchy)...) - minus_iNk(Δe)
        R2(Δe) = hte(k3b(Δe); NamedTuple{(:a⁻¹, :r, :N, :ξ)}(ht_cauchy)...) - minus_iNk(Δe)
        R3(Δe) = hte(k3b(Δe), ht_cauchy) - minus_iNk(Δe)
        # 
        R(Δe) = D(Δe) - minus_iNk(Δe)
        # 
        plot!(sp=1, real ∘ R ∘ shift, xlim..., lw=2, c=1, title="Re[Rˣ]", lab="")
        plot!(sp=2, imag ∘ R ∘ shift, xlim..., lw=2, c=1, title="Im[Rˣ]", lab="")
        # 
        plot!(sp=1, real ∘ R1 ∘ shift, xlim..., lw=1, c=2, lab="$(round(dfi.r_fm; digits=2))")
        plot!(sp=2, imag ∘ R1 ∘ shift, xlim..., lw=1, c=2, lab="")
        #
        if i == 1
            p = plot!()
            plot_render = repr(MIME"image/svg+xml"(), p)
            yl1 = ylims(p[1])
            yl2 = ylims(p[2])
            plot!(sp=1, ylims=yl1)
            plot!(sp=2, ylims=yl2)
        end
        # 
        # plot!(sp=1, real ∘ R2 ∘ shift, xlim..., lw=1, c=3, lab="")
        # plot!(sp=2, imag ∘ R2 ∘ shift, xlim..., lw=1, c=3, lab="")
        #
        # plot!(sp=1, real ∘ R3 ∘ shift, xlim..., lw=1, c=4, lab="")
        # plot!(sp=2, imag ∘ R3 ∘ shift, xlim..., lw=1, c=4, lab="")
    end
    plot!()
end
savefig(joinpath("plots", "ope-R-circlesize.pdf"))
