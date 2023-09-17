using Markdown

md"""
# ERE: $ik N(s)$ extension

The most general form of the ERE has to include a polynomial dependence of the term $N$.

```math
\begin{align}
m^2-s- i\Sigma = R(s) - ik N(s)\,,
\end{align}
```
where both $R(s)$ and N(s) are polynomials.
"""


using Pkg
cd(joinpath(@__DIR__, ".."))
Pkg.activate(".")
Pkg.instantiate()
# 
using X2DDpi
using Parameters
using Measurements
using Interpolations
using Statistics
using LaTeXStrings
using DataFrames
using CSV

using Plots
theme(:wong2, size=(500, 350), minorticks=true, grid=false, frame=:box,
    guidefontvalign=:top, guidefonthalign=:right,
    foreground_color_legend=nothing,
    legendfontsize=9, legend=:topright, lab="",
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


function effrange_denominator_II_k3b(model; method)
    effrangepars =
        effectiverangeexpansion(
            Δe -> denominator_II(model, Eᵦˣ⁺ + Δe, δm0_val),
            Δe -> k3b(Eᵦˣ⁺ + Δe),
            method)
    # 
    (; tophysicsunits(effrangepars)..., effrangepars...)
end


steps = [0.001, 0.003, 0.007]
dfv = map(steps) do _estep
    @show _estep

    filename = joinpath("results","nominal","df_radius_cutoff=$(cutoff)_estep=$(_estep).csv")
    if isfile(filename)
        println("Using the file")
        return CSV.File(filename) |> DataFrame
    end
    
    model = let
        ch1 = πDD((m1=mπ⁺, m2=mD⁰, m3=mD⁰), BW(m=mDˣ⁺, Γ=ΓDˣ⁺), BW(m=mDˣ⁺, Γ=ΓDˣ⁺))
        iπDD2 = interpolated(
            ChannelWithIntegrationMethod(ch1, HookSqrtDalitzMapping{2}()),
            cutoff; estep=_estep)
        iπDD3 = interpolated(
            ChannelWithIntegrationMethod(ch1, HookSqrtDalitzMapping{3}()),
            cutoff; estep=_estep)
        Amplitude((iπDD2, iπDD3))
    end

    ϵ0 = abs(imag(Eᵦˣ⁺))
    efv = [0.001, 0.002, 0.003, 0.005, 0.007, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2]
    df_radius = map(efv) do ϵf
        efe = effrange_denominator_II_k3b(model;
            method=ComplexBranchPointExpansion(CircularSum(ϵf * ϵ0, 50)))
        # 
        (; ϵf, efe...)
    end |> DataFrame
    
    df_radius = DataFrame()
    CSV.write(filename, df_radius)

    df_radius
end


dfvv = map(dfv) do df
    transform(df,
        :r => ByRow(x->eval(Meta.parse(x))),
        :a⁻¹ => ByRow(x->eval(Meta.parse(x))), renamecols=false)
end

let
    plot()
    map(enumerate(dfvv)) do (c,df_radius)
        xv = sqrt.(df_radius.ϵf)
        yv = df_radius.r .|> real
        plot!(xv, yv; c, lab="")
        scatter!(xv, yv; c)
    end
    plot!()
end

df_radius

for (i,df) in enumerate(dfvv)
    df.estep .= steps[i]
end

let
    plot()
    dfall = vcat(dfvv...)
    transform!(dfall, :ϵf => x->x.>0.04)
    combine(groupby(dfall, :ϵf),
        [:r,:ϵf,:estep] => (r,ϵf,estep)-> let
        plot!(estep, (r .- r[end]) .|> real, lab="$(ϵf[1])")  
        plot!(estep, (r .- r[end]) .|> imag, ls=:dash)  
    end)
    plot!(legendtitle="ϵf")
end
