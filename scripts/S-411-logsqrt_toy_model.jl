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


"""
The toy model has only sqrt branch point on the first sheet.
The log branch point is sitting on the other sheet at `x_log`
"""
const toy_model = let
    sqrtP(x, ϕ) = sqrt(x * cis(-(ϕ + π))) * cis((ϕ + π) / 2)
    logP(x, ϕ) = log(x * cis(-(ϕ + π))) + 1im * (ϕ + π)
    # 
    x_log = 0.032 - 0.04im
    x0 = Eᵦˣ⁺
    # 
    f_I(x, ϕ0, ϕ_log) = logP(-sqrtP(x - x0, ϕ0) - sqrt(x_log - x0), ϕ_log)
    f_II(x, ϕ0, ϕ_log) = logP(sqrtP(x - x0, ϕ0) - sqrt(x_log - x0), ϕ_log)
    # 
    x -> f_I(x, -π / 2, 0)
end


md"""
## Standard threshold expansion

Assuming the $N$ to be a constant,
the circular Cauchy integrals are related to the expansion parameters of $R(s)$.

We find that the expansion of the function is not great in the complex plane in the range ~0.1 MeV.
"""

function effrange(model; method)
    effrangepars =
        effectiverangeexpansion(
            Δe -> model(Eᵦˣ⁺ + Δe),
            Δe -> k3b(Eᵦˣ⁺ + Δe),
            method)
    # 
    (; tophysicsunits(effrangepars)..., effrangepars...)
end


ϵ0 = abs(imag(Eᵦˣ⁺))
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
    :a⁻¹ => ByRow(x -> real(1e3x)) => :Re_a⁻¹_x1e3,
    :a⁻¹ => ByRow(x -> imag(1e3x)) => :Im_a⁻¹_x1e3,
    :r => ByRow(real) => :Re_r,
    :r => ByRow(imag) => :Im_r,
    :N => ByRow(real) => :Re_N,
    :N => ByRow(imag) => :Im_N
)

show(df_selected, allrows=true)