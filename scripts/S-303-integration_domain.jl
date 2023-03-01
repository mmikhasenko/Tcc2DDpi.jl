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
import Plots.PlotMeasures:mm
using LaTeXStrings

theme(:vibrant, size=(500,350), minorticks=true, grid=false, frame=:box,
    guidefontvalign=:top, guidefonthalign=:right,
    foreground_color_legend = nothing,
    legendfontsize=9, legend =:topright,
    xlim=(:auto,:auto), ylim=(:auto,:auto))
# 
# 

const ch = πDD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), BW(m=mDˣ⁺,Γ=ΓDˣ⁺), BW(m=mDˣ⁺,Γ=ΓDˣ⁺))

@recipe function f(iσ::Int, e::Number, mapmethod::AbstractDalitzMapping)
 
    @series begin
        seriestype := :scatter
        markersize --> 8
        markerstype --> :d
        markercolor --> :red
        label := "D* poles"
        # 
        ([pole_position(ch.R13)'],) #pole_position(ch.R13), 
    end

    s = e2m(e)^2
    xv = range(1e-6,1,length=703)
    yv = [0, 1]
    calv = integrandpath.(xv', yv, s, Ref(mapmethod))
    calv = getindex.(calv, iσ)
    v = vcat(permutedims(hcat(permutedims(calv), ones(size(calv,2)) .* NaN))...)
    # 
    (v,)
end

@recipe function f(e::Number, mapmethod::AbstractDalitzMapping)
 
    s = e2m(e)^2
    xv = range(1e-6,1,length=303)
    calv = integrand.(xv', xv, s, Ref(mapmethod))

    seriestype --> :heatmap
    (xv,xv, (imag.(calv)))
    # (xv,xv, map(x->sign(imag(x))*log10(abs2(x)), calv))
end

# 
integrandpath(x,y, s, mapmethod = HookSqrtDalitzMapping{3}()) =
    mapdalitz(mapmethod, (x,y), masses(ch), s)[1]

function integrand(x,y, s, mapmethod = HookSqrtDalitzMapping{3}()) # 
    (σ3,σ2) = integrandpath(x,y,s,mapmethod)
    𝔐² = X2DDpi.covertapply(X2DDpi.πDD_𝔐²_nonana3, ch, s,σ3,σ2)
    return 𝔐²
end

let e = 0.9-0.1im
    plot(layout=grid(1,2), size=(1000,400),
        plot(2, e, HookSqrtDalitzMapping{3}(), lab=""),
        plot(e, HookSqrtDalitzMapping{3}(), colorbar=false))
end

