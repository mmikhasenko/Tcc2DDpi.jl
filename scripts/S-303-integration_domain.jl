using Markdown

md"""
# Singularities in the integration domain

The path of the phase space integration has to be chosen such that no pole singularites
enter the integration domian.

The integration domain is σ2 x σ3.

We assume here that the inregrand has a pole singularity in the lower plane of σ3,
and in the upper plane of σ2.

Hence the `HookSqrtDalitzMapping{3}()` mapping is used.
It has the square hoop path in the σ3 variable.
A straight path in σ2 is used to connect end endpoints of the integral, σ3±

Two plotting recipies are defined to ease the plotting.
"""

using X2DDpi
# 
using Plots
import Plots.PlotMeasures: mm
using LaTeXStrings

theme(:vibrant, size=(500, 350), minorticks=true, grid=false, frame=:box,
    guidefontvalign=:top, guidefonthalign=:right,
    foreground_color_legend=nothing,
    legendfontsize=9, legend=:topright)
# 
# 

const ch = πDD((m1=mπ⁺, m2=mD⁰, m3=mD⁰), BW(m=mDˣ⁺, Γ=ΓDˣ⁺), BW(m=mDˣ⁺, Γ=ΓDˣ⁺))

@recipe function f(iσ::Int, e::Number, mapmethod::AbstractDalitzMapping)
    s = e2m(e)^2
    xv = range(1e-6, 1, length=703)
    yv = [0, 1]
    # 
    calv = integrandpath.(xv', yv, s, Ref(mapmethod))
    calv = getindex.(calv, iσ)
    v = vcat(permutedims(hcat(permutedims(calv), ones(size(calv, 2)) .* NaN))...)
    #
    (v,)
end

@recipe function f(e::Number, mapmethod::AbstractDalitzMapping)

    s = e2m(e)^2
    xv = range(1e-6, 1, length=303)
    calv = integrand.(xv', xv, s, Ref(mapmethod))

    seriestype --> :heatmap
    (xv, xv, (imag.(calv)))
    # (xv,xv, map(x->sign(imag(x))*log10(abs2(x)), calv))
end

# 
integrandpath(x, y, s, mapmethod=HookSqrtDalitzMapping{3}()) =
    mapdalitz(mapmethod, (x, y), masses(ch), s)[1]

function integrand(x, y, s, mapmethod=HookSqrtDalitzMapping{3}()) # 
    (σ3, σ2) = integrandpath(x, y, s, mapmethod)
    𝔐² = X2DDpi.covertapply(X2DDpi.πDD_𝔐²_nonana3, ch, s, σ3, σ2)
    return 𝔐²
end

let e = 0.98 - 0.1im
    plot(layout=grid(1, 2), size=(1000, 400))
    # 
    plot!(sp=1, 2, e, HookSqrtDalitzMapping{3}(), lab="", linealpha=0.1, lc=2)
    # 
    plot!(sp=1, 1, e, HookSqrtDalitzMapping{3}(), lab="", lc=4)
    map([(0, 0), (1, 1)]) do (x, y)
        integrandpath(x, y, e2m(e)^2, HookSqrtDalitzMapping{3}())[1]
    end .|> Complex{Float64} |> v -> scatter!(v; sp=1, c=4, ms=4, lab="")
    # 
    plot!(sp=2, e, HookSqrtDalitzMapping{3}(), colorbar=false,
        ann=(0.5, 0.5, text("conjucated pole", 10, :left)))
    annotate!(sp=2, [(0.4, 0.7, text("lower-plane pole", 10, :right))])
    # 
    plot!(sp=1, [pole_position(ch.R13)], m=(7, :star5), mc=4, lab="")
    plot!(sp=1, [pole_position(ch.R13)'], m=(7, :star5), mc=:black, lab="")
    plot!(sp=1, xlab="Re m(D⁰π⁻)", ylab="Im m(D⁰π⁻)", margin=4mm)
end
savefig(joinpath("plots", "complex-integration-path.pdf"))
