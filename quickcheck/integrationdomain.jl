using X2DDpi
# 
using Plots
import Plots.PlotMeasures:mm
using LaTeXStrings
theme(:wong2, size=(500,350), minorticks=true, grid=false, frame=:box,
    guidefontvalign=:top, guidefonthalign=:right,
    foreground_color_legend = nothing,
    legendfontsize=9, legend =:topright,
    xlim=(:auto,:auto), ylim=(:auto,:auto))
# 
# 


@recipe function f(iσ::Int, e::Number, mapmethod::AbstractDalitzMapping)
 
    @series begin
        seriestype := :scatter
        markersize --> 5
        markercolor --> :red
        label --> "D* poles"
        # 
        ([pole_position(ch.R12), pole_position(ch.R12)'],)
    end

    s = e2m(e)^2
    xv = range(1e-6,1,length=103)
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



const ch = πDD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), BW(m=mDˣ⁺,Γ=ΓDˣ⁺), BW(m=mDˣ⁺,Γ=ΓDˣ⁺))

# 
integrandpath(x,y, s, mapmethod = HookSqrtDalitzMapping{3}()) =
    mapdalitz(mapmethod, (x,y), masses(ch), s)[1]

function integrand(x,y, s, mapmethod = HookSqrtDalitzMapping{3}()) # 
    (σ3,σ2) = integrandpath(x,y,s,mapmethod)
    𝔐² = X2DDpi.covertapply(X2DDpi.πDD_𝔐²_nonana3, ch, s,σ3,σ2)
    return 𝔐²
end

let e = 0.9-0.1im
    plot(layout=grid(1,2), size=(1000,600),
        plot(2, e, HookSqrtDalitzMapping{3}()),
        plot(e, HookSqrtDalitzMapping{3}(), colorbar=false))
end

