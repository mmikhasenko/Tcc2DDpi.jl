md"""
# Visual Validate Eff. Expansion

The expansion is plotted over the amplitude. 
The notebook requires the models to be defined: should be run after with S-400.
"""

using Pkg
cd(joinpath(@__DIR__, ".."))
Pkg.activate(".")
# 
using X2DDpi
using Parameters
using Measurements
using Interpolations
using Statistics

using Plots
using LaTeXStrings
import Plots.PlotMeasures: mm

theme(:wong2, size=(500, 350), minorticks=true, grid=false, frame=:box,
    guidefontvalign=:top, guidefonthalign=:right,
    foreground_color_legend=nothing,
    legendfontsize=9, legend=:topright,
    xlim=(:auto, :auto), ylim=(:auto, :auto))

# 
dfin = readjson(joinpath("results", "nominal", "effective-range-table.json")) |> DataFrame
function d2c(x)
    !(x isa Dict{String,Any}) && return x
    return x["re"] + 1im * x["im"]
end
select!(dfin, names(dfin) .=> ByRow(d2c) .=> names(dfin))


"""
The line comes from S-400 notebook where the model are computed.
"""
const models = df.model

struct ERA
    a⁻¹::Complex{Float64}
    r::Complex{Float64}
    N::Complex{Float64}
end
# 
struct EFRMatching
    𝒜::Amplitude
    era::ERA
end
EFRMatching(𝒜, nt) = EFRMatching(𝒜, ERA(nt.a⁻¹, nt.r, nt.N))

function X2DDpi.denominator_II(era::ERA, e)
    @unpack a⁻¹, r, N = era
    return N * (a⁻¹ + r / 2 * k3b(e)^2 - 1im * k3b(e))
end

using RecipesBase
@recipe function f(m::EFRMatching)
    Δev = range(-0.01, 0.0251, length=30)
    @series begin
        calv = [denominator_II(m.era, Eᵦˣ⁺ + Δe) for Δe in Δev]
        label := "effective range"
        (Δev, real(calv))
    end
    calv = [denominator_II(m.𝒜, Eᵦˣ⁺ + Δe, δm0_val) for Δe in Δev]
    (Δev, real(calv))
end

let
    plot(size=(200 * 5, 150), layout=grid(1, 5), yaxis=false, yticks=false)
    plot!(sp=1, EFRMatching(models[1], df[1, :]), lab=string(dfin.modelnames[1]))
    plot!(sp=2, EFRMatching(models[2], df[2, :]), lab=string(dfin.modelnames[2]))
    plot!(sp=3, EFRMatching(models[3], df[3, :]), lab=string(dfin.modelnames[3]))
    plot!(sp=4, EFRMatching(models[4], df[4, :]), lab=string(dfin.modelnames[4]))
    plot!(sp=5, EFRMatching(models[5], df[5, :]), lab=string(dfin.modelnames[5]))
end
savefig(joinpath("plots", "testmatchefr.pdf"))

# 
@unpack w_matching, rho_inf =
    readjson(joinpath("results", "nominal", "effective_range.json"))["effective_range_parameters"]["technical"]
# 
N′ = 1 / (w_matching * rho_inf)
let i = 1
    plot(xlab=L"\delta' m\,\,(\mathrm{MeV})", ylab=L"\mathrm{Re}\,\mathcal{A}^{-1}(black),\,\,\mathrm{Im}\,\mathcal{A}^{-1}(red)")
    plot!(Δe -> N′ * real(denominator_II(models[i], Δe, δm0_val)), -1, 3, lab="LHCb Model", lc=:black)
    plot!(Δe -> N′ * imag(denominator_II(models[i], Δe, δm0_val)), -1, 3, lab="", lc=:red)
    # 
    @unpack a⁻¹, r, N = df[i, :]
    era = ERA(a⁻¹, r, N)
    plot!(Δe -> N′ * real(denominator_II(era, Δe)), -1, 3, lab="Eff.-range exp.", ls=:dash, lc=:black)
    plot!(Δe -> N′ * imag(denominator_II(era, Δe)), -1, 3, lab="", ls=:dash, lc=:red)
    # plot!(Δe->imag(2*denominator_II(dfin.model[4], Eᵦˣ⁺+Δe, δm0_val)),-0.01, 0.1, lab=string(dfin.modelnames[4]))
    vline!([0], lab="", c=:green, lw=1)
    # scatter!([0.0], [N′*real(denominator_II(era, 0))], ms=5, mc=:black, lab="exp. point")
    # scatter!([0.0], [N′*imag(denominator_II(era, 0))], ms=5, mc=:red, lab="")
end
savefig(joinpath("plots", "nominal", "effectiverangecauchy.pdf"))
