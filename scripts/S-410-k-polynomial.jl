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

function effrange(model; method)
    effrangepars =
        effectiverangeexpansion(
            Δe -> denominator_II(model, Eᵦˣ⁺ + Δe, δm0_val),
            Δe -> k3b(Eᵦˣ⁺ + Δe),
            method)
    # 
    (; tophysicsunits(effrangepars)..., effrangepars...)
end


ϵ0 = abs(imag(Eᵦˣ⁺))
cauchysum = ComplexBranchPointExpansion(CircularSum(ϵ0 / 20, 50))
efe2 = effrange(model2; method=cauchysum)

md"""
## Size of the effect that we are fighting
"""

df_radius = map([0.0001, 0.001, 0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1, 2]) do ϵf
    efe = effrange(model2;
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
C0(f, ϵ, N=50) = circleintegral(f, 0.0, CircularSum(ϵ, N))

C0_k = C0(cauchysum.cim.r) do Δe
    k3b(Eᵦˣ⁺ + Δe)
end
C0_D = C0(cauchysum.cim.r) do Δe
    denominator_II(model2, Eᵦˣ⁺ + Δe, δm0_val)
end
@assert efe2.N == C0_D / C0_k / (-1im)

let (model, efe) = (model2, efe2)

    invD(Δe) = denominator_II(model, Δe, δm0_val)
    expansion(Δe) = ere(k3b(Δe); a⁻¹=efe.a⁻¹, r=efe.r * (1 - 1e-10 * cis(0.1π)), N=efe.N)
    minus_iNk(Δe) = ere(k3b(Δe); a⁻¹=0, r=0, N=efe.N)
    # 
    plot(layout=grid(1, 2), size=(800, 350), title=["Real" "Imag"], grid=true)
    # 
    elim = (-1, 1) .* (ϵ0)
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
    plot(xlab="integral discretization, N")
    plot!(C0_N_data.Nv, C0_N_data.C0_kv .|> real, lab="|Re[A]|")
    plot!(C0_N_data.Nv, C0_N_data.C0_kv .|> imag .|> abs, lab="|Im[A]|")
    vline!([50], lab="N=50")
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
    plot!(C0_data.ϵv, C0_data.C0_kv .|> real, lw=1.5)
    plot!(C0_data.ϵv, (C0_data.C0_Dv ./ N0) .|> real, lw=1.5)
    # 
    plot!(C0_data.ϵv, C0_data.C0_kv .|> imag, ls=:dash, lw=1.5)
    plot!(C0_data.ϵv, (C0_data.C0_Dv ./ N0) .|> imag, ls=:dash, lw=1.5)
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