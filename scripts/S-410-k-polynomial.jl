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
using Statistics


using Plots
theme(:wong2, size=(500, 350), minorticks=true, grid=false, frame=:box,
    guidefontvalign=:top, guidefonthalign=:right,
    foreground_color_legend=nothing,
    legendfontsize=9, legend=:topright, lab="",
    xlim=(:auto, :auto), ylim=(:auto, :auto))

# 
#      _|_|
#    _|      _|    _|  _|_|_|      _|_|_|
#  _|_|_|_|  _|    _|  _|    _|  _|
#    _|      _|    _|  _|    _|  _|
#    _|        _|_|_|  _|    _|    _|_|_|


const RegressionResult = NamedTuple{(:α, :β, :SE_α, :SE_β)}
import Base: *
*(nt::RegressionResult, f::Number) = RegressionResult(collect(nt) .* f)

function Base.show(io::IO, nt::RegressionResult)
    print(io,
        "Interpolation to zero:",
        "\n\tintercept: ", round(nt.β, digits=2),
        "\n\tinterpolation error: ", round(nt.SE_β, digits=2),
        "\n\tslope: ", round(nt.α, digits=2),"\n")
end


"""
    linear_regression(xv, yv)

Perform linear regression on the given data points `xv` and `yv`.

# Arguments
- `xv::Vector`: A vector of x-values.
- `yv::Complex`: A vector of corresponding y-values.

# Returns
- `α::Complex`: The slope of the linear regression.
- `β::Complex`: The intercept of the linear regression.
- `SE_α_over_α::Complex`: The standard error of the slope as δα_re/|α| + j*δα_im/|α|.
- `SE_β_over_β::Complex`: The standard error of the intercept as δβ_re/|α| + j*δβ_im/|α|.

# Example
```julia
x = [1.0, 2.0, 3.0, 4.0, 5.0]
y = [2.1, 4.0, 6.1, 8.0, 10.2]
result = linear_regression(x, y)

# Note
The function computes the intercept and slope using the formulae derived from the method of least squares.
    It also computes the standard errors for the intercept and slope to provide a measure of the uncertainty associated with these estimates. 
"""
function linear_regression(xv, yv)
    Ex = mean(xv)
    Ey = mean(yv)
    Exy = mean(xv .* yv)
    Ex² = mean(xv .^ 2)
    Dx = Ex² - Ex^2

    # α is a slope, β is the intercept 
    α, β = [Exy - Ex * Ey, Ey * Ex² - Ex * Exy] ./ Dx

    # Compute residuals
    residuals = yv .- (β .+ α .* xv)

    # Compute residual sum of squares
    RSS_real = sum(real.(residuals) .^ 2)
    RSS_imag = sum(imag.(residuals) .^ 2)

    # Compute standard error of the regression
    n = length(xv)
    dof = n - 2
    SER_real = sqrt(RSS_real / dof)
    SER_imag = sqrt(RSS_imag / dof)

    # Compute standard error of the slope
    SE_β_real = SER_real / sqrt(sum((xv .- Ex) .^ 2))
    SE_β_imag = SER_imag / sqrt(sum((xv .- Ex) .^ 2))

    # Compute standard error of the intercept
    SE_α_real = SE_β_real * sqrt(1/n + Ex^2 / sum((xv .- Ex) .^ 2))
    SE_α_imag = SE_β_imag * sqrt(1/n + Ex^2 / sum((xv .- Ex) .^ 2))

    # SE complex
    SE_α = SE_α_real + 1im*SE_α_imag
    SE_β = SE_β_real + 1im*SE_β_imag

    return RegressionResult((; α, β, SE_α, SE_β))
end



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



#                            _|            
#    _|_|_|    _|_|      _|_|_|    _|_|    
#  _|        _|    _|  _|    _|  _|_|_|_|  
#  _|        _|    _|  _|    _|  _|        
#    _|_|_|    _|_|      _|_|_|    _|_|_|  


# πDD: Dˣ⁺ + Dˣ⁺
md"""
The model with the π⁺D⁰D⁰ system with two Dˣ⁺ resonances.
"""
const model2 = let
    ch1 = πDD((m1=mπ⁺, m2=mD⁰, m3=mD⁰), BW(m=mDˣ⁺, Γ=ΓDˣ⁺), BW(m=mDˣ⁺, Γ=ΓDˣ⁺))
    iπDD2 = interpolated(
        ChannelWithIntegrationMethod(ch1, HookSqrtDalitzMapping{2}()),
        cutoff; estep=estep)
    iπDD3 = interpolated(
        ChannelWithIntegrationMethod(ch1, HookSqrtDalitzMapping{3}()),
        cutoff; estep=estep)
    Amplitude((iπDD2, iπDD3))
end



md"""
## Standard threshold expansion

Assuming the $N$ to be a constant,
the circular Cauchy integrals are related to the expansion parameters of $R(s)$.

We find that the expansion of the function is not great in the complex plane in the range ~0.1 MeV.
"""

function effrange_denominator_II_k3b(model, δm0; method)
    effrangepars =
        effectiverangeexpansion(
            Δe -> denominator_II(model, Eᵦˣ⁺ + Δe, δm0),
            Δe -> k3b(Eᵦˣ⁺ + Δe),
            method)
    # 
    (; tophysicsunits(effrangepars)..., effrangepars...)
end


ϵ0 = abs(imag(Eᵦˣ⁺))

md"""
## Size of the effect 

Before further studies, let check how strong the dependce on r, that we are fighting!
"""

df_radius = map([0.0001, 0.001, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2]) do ϵf
    efe = effrange_denominator_II_k3b(model2, δm0_val;
        method=ComplexBranchPointExpansion(CircularSum(ϵf * ϵ0, 50)))
    # 
    (; ϵf, efe...)
end |> DataFrame

select(df_radius, :ϵf,
    :a⁻¹ => ByRow(x -> real(1e3x)) => :Re_a⁻¹_x1e3,
    :a⁻¹ => ByRow(x -> imag(1e3x)) => :Im_a⁻¹_x1e3,
    :r => ByRow(real) => :Re_r,
    :r => ByRow(imag) => :Im_r,
    :N => ByRow(x -> real(1e3x)) => :Re_N_x1e3,
    :N => ByRow(x -> imag(1e3x)) => :Im_N_x1e3
) |> DataFrames.PrettyTables.pretty_table



md"""
### Consistency check

We are going to expand
"""

cauchysum = ComplexBranchPointExpansion(CircularSum(ϵ0 / 20, 50))
efe2 = effrange_denominator_II_k3b(model2, δm0_val; method=cauchysum)

C0(f, ϵ, N=50) = circleintegral(f, 0.0, CircularSum(ϵ, N))

C0_k = C0(cauchysum.cim.r) do Δe
    k3b(Eᵦˣ⁺ + Δe)
end
C0_D = C0(cauchysum.cim.r) do Δe
    denominator_II(model2, Eᵦˣ⁺ + Δe, δm0_val)
end
@assert efe2.N == C0_D / C0_k / (-1im)


import Plots.PlotMeasures: mm
let
    selected_indices = 4:size(df_radius,1)-2
    # 
    ϵfv = df_radius.ϵf[selected_indices]
    xv = sqrt.(ϵfv)
    yv = df_radius.r[selected_indices]
    regression = linear_regression(xv, yv)
    @unpack α, β = regression
    # 
    plot(layout=grid(1,2), size=(900,400),
        title = "Regression",
        ylab=["Re r" "Im r"] .* " [fm]", xlab="scan scale, √ϵf",
        bottom_margin = 3mm, left_margin = 4mm)
    # 
    plot!(sp=1, x-> real(β + α * x), range(xv[[1,end]]..., 30))
    scatter!(sp=1, xv, yv .|> real, m=(7, :d))
    #
    plot!(sp=2, x-> imag(β + α * x), range(xv[[1,end]]..., 30))
    scatter!(sp=2, xv, yv .|> imag, m=(7, :d))
    #
    plot!()
end



let
    model = model2
#     efe = df_radius[4, :]
    efe = efe2
#     dx = efe.ϵf * ϵ0
    dx = ϵ0
    elim = (-1, 1) .* (2dx)

    invD(Δe) = denominator_II(model, Δe, δm0_val)
    expansion(Δe) = ere(k3b(Δe); a⁻¹=efe.a⁻¹, r=efe.r * (1 - 1e-10 * cis(0.1π)), N=efe.N)
    minus_iNk(Δe) = ere(k3b(Δe); a⁻¹=0, r=0, N=efe.N)
    # 
    plot(layout=grid(1, 2), size=(800, 350), title=["Real" "Imag"], grid=true)
    # 
    xv = range(elim..., 55)
    invD_yv = map(xv .+ Eᵦˣ⁺) do Δe
        invD(Δe) - minus_iNk(Δe)
    end
    efe_yv = map(xv .+ Eᵦˣ⁺) do Δe
        expansion(Δe) - minus_iNk(Δe)
    end
    plot!(sp=1, xv, real(invD_yv), lab="1/D")
    plot!(sp=1, xv, real(efe_yv), c=1, ls=:dash, lab="Expansion")
    # 
    plot!(sp=2, xv, imag(invD_yv), lab="1/D")
    plot!(sp=2, xv, imag(efe_yv), c=1, ls=:dash, lab="Expansion")
end

md"""
The plot shows a striking difference between the expansion series and the target function.
 - The ERE is a ~linear on energy
 - The target function has a significant non-linear behavior.
 - The scale of the deviation is 
"""


md"""
## Explore how well branch cut is approximated with the $N=$const assumption.
"""

md"""
### Exploration of the integral convergence
"""

C0_N_data = let ϵ = abs(imag(Eᵦˣ⁺)) / 2
    Nv = 10:10:70
    C0_kv = map(Nv) do N
        C0(ϵ, N) do Δe
            k3b(Eᵦˣ⁺ + Δe)
        end
    end
    (; Nv, C0_kv)
end

let
    plot(xlab="integral discretization, N",
        title="Circle integral discretization")
    plot!(C0_N_data.Nv, 1e6 * C0_N_data.C0_kv .|> real, lab="|Re A|")
    plot!(C0_N_data.Nv, 1e6 * C0_N_data.C0_kv .|> imag .|> abs, lab="|Im A|")
    vline!([50], lab="N=50", ylab="C0 × 10⁶")
end

md"""
Convergence of the circular integral is checked by varying the number of sampled points.
The used N=50 seems to be a good compromise.
"""


md"""
Now, the circle radius is changed modifying the integral path over the discontinuity.
The comparison with the model tells how good the N=cost approximation is.
"""

C0_data = let
    ϵv = [0.1, 0.2, 0.5, 1, 2]
    C0_kv = map(ϵv) do ϵ
        C0(ϵ) do Δe
            k3b(Eᵦˣ⁺ + Δe)
        end
    end .* 1e3
    C0_Dv = map(ϵv) do ϵ
        C0(ϵ) do Δe
            denominator_II(model2, Eᵦˣ⁺ + Δe, δm0_val)
        end
    end ./ (-1im) .* 1e3
    (; ϵv, C0_kv, C0_Dv)
end


let
    N0 = C0_data.C0_Dv[end] / C0_data.C0_kv[end]
    plot(layout=grid(3, 1, heights=(0.6, 0.2, 0.2)), size=(700, 700),
        xlab="Cauchy integral radius, ϵ")
    plot!(sp=1, C0_data.ϵv, C0_data.C0_kv .|> real, lw=1.5, lab="Re C₀(k)")
    plot!(sp=1, C0_data.ϵv, (C0_data.C0_Dv ./ N0) .|> real, lw=1.5, lab="Re C₀(A)")
    plot!(sp=1, title="Comparison of C₀ integral of A and k")
    # 
    plot!(sp=1, C0_data.ϵv, C0_data.C0_kv .|> imag, ls=:dash, lw=1.5, lab="Im C₀(k)")
    plot!(sp=1, C0_data.ϵv, (C0_data.C0_Dv ./ N0) .|> imag, ls=:dash, lw=1.5, lab="Im C₀(A)")
    # 
    plot!(sp=2, C0_data.ϵv, (C0_data.C0_Dv ./ N0 .- C0_data.C0_kv) .|> real, lc=2, lw=1.5)
    plot!(sp=2, C0_data.ϵv, (C0_data.C0_Dv ./ N0 .- C0_data.C0_kv) .|> imag, lc=4, lw=1.5, ls=:dash)
    # 
    plot!(sp=3, C0_data.ϵv, (C0_data.C0_Dv ./ C0_data.C0_kv) .|> imag, lc=2, lw=1.5)
    # N0
end

md"""
The discontinuity across the cut matches the ikN with N being constant well.
The deviation is of an order

|δA/A| < 1e-7 / 1e-4

is seen. It is much smaller the deviation effect in Re Im pars.
"""


md"""
# Summary

The dependence on the integral radius is natural and understood.
The correction is linear in √circle_r. The functions below do the intepolation

The extraction algorithm will be:
1. build model,
2. scan in ϵf with `effective_range_scan`,
3. extrapolate to zero with `effective_range_extrapolation`
"""


function effective_range_scan(model, settings)
    @unpack ϵ0 = settings
    @unpack ϵf_grid = settings
    @unpack circular_sum_n = settings
    @unpack model_δm0 = settings
    # 
    method(r) = ComplexBranchPointExpansion(CircularSum(r, circular_sum_n))

    scan = map(ϵf_grid) do ϵf
        circle_r = ϵf * ϵ0
        effrangepars =
            effectiverangeexpansion(
                Δe -> denominator_II(model, Eᵦˣ⁺ + Δe, model_δm0),
                Δe -> k3b(Eᵦˣ⁺ + Δe),
                method(circle_r))
        efe = (; tophysicsunits(effrangepars)..., effrangepars..., ϵf, circle_r)
        efe
    end

    return scan
end

function linear_regression_in_sqrt_ϵf(scan, quantity)
    ϵfv = getproperty.(scan, :ϵf)
    sqrt_ϵfv = sqrt.(ϵfv)
    # 
    calv = getproperty.(scan, quantity)
    regression = linear_regression(sqrt_ϵfv, calv)
    # 
    return regression
end

function effective_range_extrapolation(scan)
    r_regression = linear_regression_in_sqrt_ϵf(scan, :r)
    inva_regression = linear_regression_in_sqrt_ϵf(scan, :a⁻¹)
    # 
    (; r_regression, inva_regression)
end


# apply to the model as a test

effective_range_intepolation_settings = ( 
    ϵ0 = abs(imag(Eᵦˣ⁺)),
    ϵf_grid = [0.02, 0.05, 0.1, 0.2, 0.5],
    circular_sum_n = 50,
    model_δm0 = δm0_val)

ere_scan = effective_range_scan(model2,
    effective_range_intepolation_settings)
# 
ssc = effective_range_extrapolation(ere_scan)

@show ssc.r_regression;
@show ssc.inva_regression * 1000;
