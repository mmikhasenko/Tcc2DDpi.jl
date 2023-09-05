using Markdown

md"""
# ERE for a toy log ∘ sqrt function
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

using Plots
theme(:wong2, size=(500, 350), minorticks=true, grid=false, frame=:box,
    guidefontvalign=:top, guidefonthalign=:right,
    foreground_color_legend=nothing,
    legendfontsize=9, legend=:topright, lab="",
    xlim=(:auto, :auto), ylim=(:auto, :auto))


#                            _|            
#    _|_|_|    _|_|      _|_|_|    _|_|    
#  _|        _|    _|  _|    _|  _|_|_|_|  
#  _|        _|    _|  _|    _|  _|        
#    _|_|_|    _|_|      _|_|_|    _|_|_|  


logP(x, ϕ) = log(x * cis(-(ϕ + π))) + 1im * (ϕ + π)
sqrtP(x, ϕ) = sqrt(x * cis(-(ϕ + π))) * cis((ϕ + π) / 2)

"""
The toy model has only sqrt branch point on the first sheet.
The log branch point is sitting on the other sheet at `x_log`
"""
const toy_model = let
    x_log = 1 - 1im
    x0 = 0
    # 
    f_I(x, ϕ0, ϕ_log) = logP(-sqrtP(x - x0, ϕ0) - sqrt(x_log - x0), ϕ_log)
    f_II(x, ϕ0, ϕ_log) = logP(sqrtP(x - x0, ϕ0) - sqrt(x_log - x0), ϕ_log)
    # 
    x -> f_I(x, -π / 2, 0)
end


heatmap(range(-3, 3, 200), range(-3, 3, 200), (x, y) -> imag(toy_model(x + 1im * y)))
heatmap(range(-3, 3, 200), range(-3, 3, 200), (x, y) -> imag(sqrtP(x + 1im * y, -π / 2)))

md"""
## Standard threshold expansion

Assuming the $N$ to be a constant,
the circular Cauchy integrals are related to the expansion parameters of $R(s)$.

We find that the expansion of the function is not great in the complex plane in the range ~0.1 MeV.
"""

function effrange(model; method)
    effrangepars =
        effectiverangeexpansion(
            x -> model(x),
            x -> sqrtP(x, -π / 2),
            method)
    # 
    (; tophysicsunits(effrangepars)..., effrangepars...)
end


ϵ0 = 1.0
cauchysum = ComplexBranchPointExpansion(CircularSum(ϵ0 / 20, 50))
efe2 = effrange(toy_model; method=cauchysum)

"""
## The function shows the same effect
"""

df_radius = map([0.0001, 0.001, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2]) do ϵf
    efe = effrange(toy_model;
        method=ComplexBranchPointExpansion(CircularSum(ϵf * ϵ0, 50)))
    # 
    (; ϵf, efe...)
end |> DataFrame

df_selected = select(df_radius, :ϵf,
    :a⁻¹ => ByRow(x -> real(x)) => :Re_a⁻¹_x1e3,
    :a⁻¹ => ByRow(x -> imag(x)) => :Im_a⁻¹_x1e3,
    :r => ByRow(real) => :Re_r,
    :r => ByRow(imag) => :Im_r,
    :N => ByRow(real) => :Re_N,
    :N => ByRow(imag) => :Im_N
)


ana_ere = let
    b = sqrt(1 - 1im)
    N = 1im / b
    a⁻¹ = log(-b) / N
    r = -1 / b^2 / N
    (; N, a⁻¹, r)
end


let
    k(z) = sqrtP(z, -π / 2)
    path(t) = t * cis(-π)
    ere_series(k, pars) = pars.N * (pars.a⁻¹ + pars.r / 2 * k^2 - 1im * k)
    # 
    plot(legtitle="r")
    plot!(x -> x |> path |> toy_model |> real, 0, 1, lab="", lw=2)
    Nrow = size(df_radius, 1)
    for (i, ere) in enumerate(eachrow(df_radius[1:Nrow-1, :]))
        # ere = df_radius[5, :]
        ere_model = Base.Fix2(ere_series, ere)
        rRe = round(real(ere.r), digits=2)
        plot!(x -> x |> path |> k |> ere_model |> real, 0, 1, lab="$(rRe)", c=cgrad(:matter)[(i-1)/(Nrow-1)])
    end
    ere_model_ana = Base.Fix2(ere_series, ana_ere)
    plot!(x -> x |> path |> k |> ere_model_ana |> real, 0, 1, c=:red, lw=2)
    # 
    plot!(xlab="path variable, t", ylab="f")
end

md"""
The tailor expansion gives an approximation of the function in small vicinity of
the branch point. In presence of additional singularities in the target function,
the difference between the function and its series blows up quickly.
"""