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
"""
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
    effrangepars = # 12s
        effectiverangeexpansion(
            Δe -> denominator_II(model, Eᵦˣ⁺ + Δe, δm0_val),
            Δe -> k3b(Eᵦˣ⁺ + Δe),
            method)
    # 
    (; tophysicsunits(effrangepars)..., effrangepars...)
end


cauchysum = ComplexBranchPointExpansion(CircularSum(abs(imag(Eᵦˣ⁺)) / 2, 50))
efe2 = effrange(model2; method=cauchysum)
@assert efe2.N == C0_D / C0_k / (-1im)

let (model, efe) = (model2, efe2)

    invD(Δe) = denominator_II(model, Δe, δm0_val)
    expansion(Δe) = ere(k3b(Δe), efe)
    minus_iNk(Δe) = ere(k3b(Δe); a⁻¹=0, r=0, N=efe.N)
    # 
    plot(layout=grid(1, 2), size=(800, 350), title=["Real" "Imag"])
    # 
    elim = (-1, 1)
    plot!(sp=1, Δe -> real(invD(Δe) - minus_iNk(Δe)), elim..., lab="1/D")
    plot!(sp=1, Δe -> real(expansion(Δe) - minus_iNk(Δe)), elim..., c=1, ls=:dash, lab="Expansion")
    # 
    plot!(sp=2, Δe -> imag(invD(Δe) - minus_iNk(Δe)), elim..., lab="1/D")
    plot!(sp=2, Δe -> imag(expansion(Δe) - minus_iNk(Δe)), elim..., c=1, ls=:dash, lab="Expansion")
end


md"""
## Explore how well branch cut is approximated with the $N=$const assumption.
"""

C0(f, ϵ=abs(imag(Eᵦˣ⁺)) / 2, N=50) = circleintegral(f, 0.0, CircularSum(ϵ, N))

md"""
### Consistency check
"""

C0_k = C0() do Δe
    k3b(Eᵦˣ⁺ + Δe)
end

C0_D = C0() do Δe
    denominator_II(model2, Eᵦˣ⁺ + Δe, δm0_val)
end

md"""
### Exploration of dependence
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

plot(layout=grid(1, 2), size=(900, 350), ylim=(0, :auto),
    plot(C0_N_data.Nv, C0_N_data.C0_kv .|> real, lc=2),
    plot(C0_N_data.Nv, C0_N_data.C0_kv .|> imag .|> abs, lc=2))

C0_data = let
    ϵv = range(1e-4, abs(imag(Eᵦˣ⁺)) * 2, 10)
    C0_kv = map(ϵv) do ϵ
        C0(ϵ) do Δe
            k3b(Eᵦˣ⁺ + Δe)
        end
    end
    C0_Dv = map(ϵv) do ϵ
        C0(ϵ) do Δe
            denominator_II(model2, Eᵦˣ⁺ + Δe, δm0_val)
        end
    end ./ (-1im)
    (; ϵv, C0_kv, C0_Dv)
end


let
    N0 = C0_data.C0_Dv[end] / C0_data.C0_kv[end]
    plot(layout=grid(3, 1, heights=(0.6, 0.2, 0.2)), size=(700, 700))
    plot!(C0_data.ϵv, C0_data.C0_kv .|> real)
    plot!(C0_data.ϵv, (C0_data.C0_Dv ./ N0) .|> real)
    # 
    plot!(C0_data.ϵv, C0_data.C0_kv .|> imag, ls=:dash)
    plot!(C0_data.ϵv, (C0_data.C0_Dv ./ N0) .|> imag, ls=:dash)
    # 
    plot!(sp=2, C0_data.ϵv, (C0_data.C0_Dv ./ N0 .- C0_data.C0_kv) .|> real, lc=2)
    plot!(sp=3, C0_data.ϵv, (C0_data.C0_Dv ./ N0 .- C0_data.C0_kv) .|> imag, lc=4)
end

