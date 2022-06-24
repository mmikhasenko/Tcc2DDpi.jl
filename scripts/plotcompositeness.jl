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
using LaTeXStrings

theme(:wong2, frame=:box, grid=false, minorticks=true, 
    guidefontvalign=:top, guidefonthalign=:right)


effectiverangeresults = readjson(joinpath("results","nominal","effectiverangecauchy.json"))


@unpack Z_polosa_90, Z_polosa_95,
    Z_hanhart_90, Z_hanhart_95 = 
        effectiverangeresults["compositeness"]
# 
import Plots.PlotMeasures:mm
let
    plot(layout=grid(2,1, heights=(0.9,0.1)), link=:x, size=(500,150),
        yaxis=false, yticks=false, frame=:origin, xlim=(0,1)
        )
    plot!(sp=2, ylim=(-0.05,0.05), xticks=(collect(0:0.2:1.0), string.(1.0:-0.2:0.0)),
        ann=(1,0.05, text(L"1-Z",12,:bottom,:right)))
    plot!(sp=1, ylim=(-0.6,0.6),
        ann=(1,0.1, text(L"X",12,:bottom,:right)))
    plot!(sp=1, [Z_polosa_90...], 0.5 .* [1, 1], fill=0, c=2, α=0.2, lab="")
    plot!(sp=1, [Z_polosa_95...], 0.5 .* [0.8, 0.8], fill=0, c=2, α=0.2, lab="")
    annotate!(sp=1, [
        (Z_polosa_90[1], 0.5, text(string(round(Z_polosa_90[1],digits=2)),:bottom, 8)),
        (Z_polosa_90[2], 0.5, text(string(round(Z_polosa_90[2],digits=2)),:bottom, 8)),
        (Z_polosa_95[1], 0.4, text(string(round(Z_polosa_95[1],digits=2)),:bottom, 8)),
        ])
    plot!(sp=1, [Z_hanhart_90...], 0.5 .* [-1, -1], fill=0, c=3, α=0.2, lab="")
    plot!(sp=1, [Z_hanhart_95...], 0.5 .* [-0.8, -0.8], fill=0, c=3, α=0.2, lab="")
    annotate!(sp=1, [
        (Z_hanhart_90[1], -0.5, text(string(round(Z_hanhart_90[1],digits=2)),:top, 8)),
        (Z_hanhart_90[2], -0.5, text(string(round(Z_hanhart_90[2],digits=2)),:top, 8)),
        (Z_hanhart_95[1], -0.4, text(string(round(Z_hanhart_95[1],digits=2)),:top, 8)),
        ])
    plot!(sp=1, [0,0], [-0.06,0.06], c=:black, lab="", ann=(0.01,-0.06,text("0.0",:top,8)))
    plot!(sp=2, [0,0], [-0.06,0.06], c=:black, lab="", ann=(0.01,-0.06,text("1.0",:top,8)))
    annotate!(sp=1, [
        (0, -0.15, text(latexstring("[D^{*0}D^+\\mathrm{removed}]"), :top,    :left, 11))])#,
        # (0,  0.15, text(latexstring("[\\mathrm{experimental}\\,\\,r,a^{-1}]") ,:bottom, :left, 11))])
end
savefig(joinpath("plots", "compositeness.pdf"))

