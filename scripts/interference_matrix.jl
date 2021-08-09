using Pkg
cd(joinpath(@__DIR__,".."))
Pkg.activate(".")
# 
using X2DDpi
using Parameters
using Measurements

# 
using Plots
using LaTeXStrings
theme(:wong2, size=(500,350), minorticks=true, grid=false, frame=:box,
    guidefontvalign=:top, guidefonthalign=:right,
    foreground_color_legend = nothing,
    legendfontsize=9, legend =:topright,
    xlim=(:auto,:auto), ylim=(:auto,:auto))


using AlgebraPDF
# 


settings = transformdictrecursively!(readjson("settings.json"), ifstringgivemeasurement)
#
@unpack δm0 = settings["fitresults"]

channels = [
    πDD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), (m=mDˣ⁺,Γ=ΓDˣ⁺), (m=mDˣ⁺,Γ=ΓDˣ⁺)),
    πDD((m1=mπ⁰,m2=mD⁺,m3=mD⁰), (m=mDˣ⁺,Γ=ΓDˣ⁺), (m=mDˣ⁰,Γ=ΓDˣ⁰)),
    γDD((m1=mγ, m2=mD⁺,m3=mD⁰), (m=mDˣ⁺,Γ=ΓDˣ⁺), (m=mDˣ⁰,Γ=ΓDˣ⁰))]


import X2DDpi: J_I, J_II
J_I(σ::Float64,pars::NamedTuple{(:zero,)}) = 0.0
J_II(σ::Float64,pars::NamedTuple{(:zero,)}) = 0.0
# 
function J_I(σ,pars::NamedTuple{(:im,:iΓ)})
	m, Γ = pars
	FF = 1.0 # (σ-(mD⁰+mπ⁺)^2) / (m^2-(mD⁰+mπ⁺)^2)
	1im/(m^2 - σ + 1im*m*Γ*FF)
end
function J_II(σ,pars::NamedTuple{(:im,:iΓ)})
	m, Γ = pars
	FF = 1.0 # (σ-(mD⁰+mπ⁺)^2) / (m^2-(mD⁰+mπ⁺)^2)
	-1im/(m^2 - σ - 1im*m*Γ*FF)
end

const mΓDˣ⁺ = (m=mDˣ⁺, Γ=ΓDˣ⁺)
const mΓDˣ⁰ = (m=mDˣ⁰, Γ=ΓDˣ⁰)
const noR = (zero=0,)
# 
const imΓDˣ⁺ = (im=mDˣ⁺, iΓ=ΓDˣ⁺)
const imΓDˣ⁰ = (im=mDˣ⁰, iΓ=ΓDˣ⁰)
# 
function scattmatrix(all_12_13)
    F, a, b, f12, f13 = all_12_13
    Δ = (F-a-b) - 1im*(f12-a-b)
    [a Δ/2; Δ'/2 b]
end
# 
function normmatrix(all_12_13)
    F, a, b, f12, f13 = all_12_13
    Δ = (F-a-b) - 1im*(f12-a-b)
    [1 Δ/2/sqrt(a*b); Δ'/2/sqrt(a*b) 1]
end
# 
function add_error(m0, m1, m2)
    Δ = (m1-m2)/2
    .±(real.(m0),real.(Δ)) + 
      .±(imag.(m0),imag.(Δ)).*1im
end

function matrix_and_interference(expfunction)
    C = add_error(scattmatrix(expfunction(δm0.val)),
                  scattmatrix(expfunction(δm0.val-δm0.err)),
                  scattmatrix(expfunction(δm0.val+δm0.err)))

    P = add_error(normmatrix(expfunction(δm0.val)),
                  normmatrix(expfunction(δm0.val-δm0.err)),
                  normmatrix(expfunction(δm0.val+δm0.err)))
    #
    return (; C, I = P[1,2])
end

total_interference = x->2real(x)/(2+2real(x))

#
##########################################################
# 
# 
MI_π⁺D⁰D⁰ = maexprho_π⁺D⁰D⁰(δm) = [1e6*ρ_thr(πDD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), ξs...), δm)
for ξs in [(mΓDˣ⁺, mΓDˣ⁺), (mΓDˣ⁺,noR), (noR, mΓDˣ⁺), (imΓDˣ⁺, mΓDˣ⁺), (mΓDˣ⁺, imΓDˣ⁺)]]
# 
trix_and_interference(exprho_π⁺D⁰D⁰)
MI_π⁺D⁰D⁰.I |> total_interference
#
##########################################################
# 
# 
exprho_π⁰D⁺D⁰(δm) = [1e6*ρ_thr(πDD((m1=mπ⁰,m2=mD⁺,m3=mD⁰), ξs...), δm)
    for ξs in [(mΓDˣ⁺, mΓDˣ⁰), (mΓDˣ⁺, noR), (noR, mΓDˣ⁰), (imΓDˣ⁺, mΓDˣ⁰), (mΓDˣ⁺, imΓDˣ⁰)]]
# 
MI_π⁰D⁺D⁰ = matrix_and_interference(exprho_π⁰D⁺D⁰)
MI_π⁰D⁺D⁰.I |> total_interference
#
##########################################################
# 
# 
exprho_γD⁺D⁰(δm) = [1e6*ρ_thr(γDD((m1=mγ, m2=mD⁺,m3=mD⁰), ξs...), δm)
    for ξs in [(mΓDˣ⁺, mΓDˣ⁰), (mΓDˣ⁺, noR), (noR, mΓDˣ⁰), (imΓDˣ⁺, mΓDˣ⁰), (mΓDˣ⁺, imΓDˣ⁰)]]
#
MI_γD⁺D⁰ = matrix_and_interference(exprho_γD⁺D⁰)

MI_γD⁺D⁰.I |> total_interference