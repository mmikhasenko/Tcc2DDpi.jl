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

const πi = π * 1im

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


@with_kw struct LogSqrtModel
    x_log::ComplexF64
    x_sqrt::ComplexF64
    ϕ_log::Float64
    ϕ_sqrt::Float64
end

f_I(x, x0, x_log, ϕ0, ϕ_log) = logP(-sqrtP(x - x0, ϕ0) - sqrt(x_log - x0), ϕ_log)
f_II(x, x0, x_log, ϕ0, ϕ_log) = logP(sqrtP(x - x0, ϕ0) - sqrt(x_log - x0), ϕ_log)
(m::LogSqrtModel)(x) = f_I(x, m.x_sqrt, m.x_log, m.ϕ_sqrt, m.ϕ_log)

toy_model = LogSqrtModel(
    x_log=1 - 1e-6im,
    x_sqrt=0,
    ϕ_log=0,
    ϕ_sqrt=-π / 2
)


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
        method=ComplexBranchPointExpansion(CircularSum(ϵf * ϵ0, 1050)))
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
    b = toy_model.x_log
    N = 1im / b
    a⁻¹ = log(-b) / N
    r = -1 / b^2 / N
    (; N, a⁻¹, r)
end


let
    k(z) = sqrtP(z, -π / 2)
    path(t) = t * cis(3π / 2)
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

md"""
## Correction terms

for the function R = log(√x+1)-√x,
one can compute the integral of R/x^2 analytically and plot corrections vs circle radius
"""


indefint0(x) = -x / 2 + sqrt(x) + (x - 1) * log(sqrt(x) + 1)
disc_indefint0(x) = (2 * sqrt(x) + (x - 1) * log(sqrt(x) + 1) - (x - 1) * log(-sqrt(x) + 1))




indefint2(x) = -1 / sqrt(x) + (x - 1) * log(sqrt(x) + 1) / x - log(x) / 2

"""
indefint2_no_sqrt = Int [ log(sqrt(x) + 1) - sqrt(x) ] dx

given that,

Int sqrt(x)/x^2 dx = -2/sqrt(x)

we removed 

"""
indefint2_no_sqrt(x) = indefint2(x) + 2 / sqrt(x)
# 
"""
    disc_indefint2_no_sqrt

Note: -2πi / 2 is given by disc of log(x) / 2
"""
disc_indefint2_no_sqrt(x) =
    2 / sqrt(x) +
    (x - 1) * log(sqrt(x) + 1) / x -
    (x - 1) * log(-sqrt(x) + 1) / x -
    (-2πi) / 2 # clockwise


let cim = CircularSum(0.4ϵ0, 50)
    N(ϵ) = disc_indefint0(ϵ) / (ϵ * sqrt(ϵ) * 4 / 3) / (-1im)
    N_ = N(-1im * cim.r)
    cauchyintegral′(x -> toy_model(x) / N_ + 1im * sqrtP(x, toy_model.ϕ_sqrt), 0.0, cim),
    -disc_indefint2_no_sqrt(-1im * cim.r) / 2πi / N_
end


effrange(toy_model; method=ComplexBranchPointExpansion(CircularSum(0.4 * ϵ0, 50)))

plot(c, 0, 1)
plot!(x -> x * sqrt(x) * 4 / 3, 0, 1)

plot!(x -> c(x) / (x * sqrt(x) * 4 / 3), 0, 2)

let ϵ = -1im * cauchysum.cim.r
    disc_indefint0(ϵ) / (ϵ * sqrt(ϵ) * 4 / 3) / (-1im)
end

let df = df_radius[1:end-1, :]
    xv = df.ϵf .* ϵ0
    plot()
    # plot!(xv, df.N .|> imag, lab="num")
    plot!(xv, imag.(df.r), lab="num")
    # 
    N(ϵ) = disc_indefint0(ϵ) / (ϵ * sqrt(ϵ) * 4 / 3) / (-1im)
    R′(ϵ) = -disc_indefint2_no_sqrt(ϵ) / (2πi) / N(ϵ)
    r(ϵ) = 2 * R′(ϵ)
    # 
    # plot!(ϵ -> -1im * ϵ |> N |> imag, xv, lab="ana")
    @show N(-0.4im)
    @show R′(-0.4im)
    @show r(-0.4im)

    plot!(x -> imag(r(-1im * x)), xv, lab="ana")
    plot!()
end


md"""
### Test with sqrt(x): F = -2/sqrt(x)
"""
cauchyintegral′(x -> sqrtP(x, -π / 2), 0.0, cauchysum.cim)

let ϵ = cauchysum.cim.r
    F(x) = -2 / sqrtP(x, -π / 2)
    ΔF = (F(-1im * ϵ + 1e-7) - F(-1im * ϵ - 1e-7)) / 2πi
    disc(x) = -4 / sqrt(x) / 2πi
    ΔF, disc(-1im * ϵ)
end

md"""
everything matches given the direction anticlockwise here
"""

