using Pkg
cd(joinpath(@__DIR__,".."))
Pkg.activate(".")

using X2DDpi

using Plots
using LaTeXStrings
theme(:wong2, size=(450,400), minorticks=true, grid=false, frame=:box,
    guidefontvalign=:top, guidefonthalign=:right,
    foreground_color_legend = nothing,
    legendfontsize=9, legend =:topright,
    xlim=(:auto,:auto), ylim=(:auto,:auto),
    xlab=L"m_{D^0\pi^+}^2\,\,[\mathrm{GeV}]", ylab=L"m_{D^0\pi^+}^2\,\,[\mathrm{GeV}]")


# 
using Parameters
using Measurements
# 
using ThreeBodyDecay
using PartialWaveFunctions


settings = transformdictrecursively!(readjson("settings.json"), ifstringgivemeasurement)
#
@unpack δm0 = settings["fitresults"]

# retrieve the model
modelDict = readjson(joinpath("results","nominal","model.json"))
ichannels = interpolated.(d2nt.(modelDict["ichannels"]))
channels = getproperty.(ichannels, :channel)

const ch1 = πDD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), (m=mDˣ⁺,Γ=ΓDˣ⁺), (m=mDˣ⁺,Γ=ΓDˣ⁺))

const m0 = e2m(δm0.val) 
const s0 = e2m(δm0.val)^2 
const ms = ThreeBodyMasses(ch1.ms...;m0)

###################################################################

import X2DDpi: decay_matrix_element_squared
import X2DDpi: AbstractxDD

@with_kw struct πDD_heli <: AbstractxDD
    j0::Int = 1
    j::Int = 1
    L::Int = 0
    plus::Bool = true
end

pDπ(σ) = σ > (mD⁰+mπ⁺)^2 ? sqrt(λ(σ,mD⁰^2,mπ⁺^2)/(4*σ)) : 0.0
pDˣD(σ) = σ < (m0-mD⁰)^2 ? sqrt(λ(s0,σ,mD⁰^2)/(4*s0)) : 0.0

function decay_matrix_element_squared(d::πDD_heli,s,σ3,σ2)
    s != s0 && println("\nHm, s != s0 !\n")
    # 
    @unpack j, j0, L, plus = d
    # 
	J12_I, J12_II = J_I(σ3,ch1.R12), J_II(σ3,ch1.R12)
	J13_I, J13_II = J_I(σ2,ch1.R13), J_II(σ2,ch1.R13)
    p12_j = pDπ(σ3)^j; p3_L = (L==0 ? 1.0 : pDˣD(σ3)^L)
    p13_j = pDπ(σ2)^j; p2_L = (L==0 ? 1.0 : pDˣD(σ2)^L)
    # 
    msq4 = ms^2
    σs = Invariants(ms; σ3,σ2)

    zhat23 = cosθhat23(σs,msq4)
    z12 = cosθ12(σs,msq4)
    z31 = cosθ31(σs,msq4)
    H(ν) = clebschgordan(L,0,j,ν,j0,ν)


    𝔐² =
        sum(abs2, wignerd(j,λ,0,z12)*H(λ)*p12_j*p3_L for λ in -j:j) * J12_I*J12_II +
        sum(abs2, wignerd(j,λ,0,z31)*H(λ)*p13_j*p2_L for λ in -j:j) * J13_I*J13_II +
            (plus ? 1 : -1) * (J13_I * J12_II +  J12_I * J13_II) *
            sum(wignerd(j,λ,0,z12)*H(λ) * p12_j * p3_L *
                wignerd(j,ν,0,z31)*H(ν) * p13_j * p2_L *
                wignerd(j0,λ,ν,zhat23)
                    for (ν,λ) in Iterators.product(-j:j,-j:j))
    return 𝔐²
end


hel = πDD_heli(;plus=false)
let σs = randomPoint(ms)
    decay_matrix_element_squared(hel,s0, σs.σ3,σs.σ2)
end

dalitzplot(hel; what2apply=real)
dalitzplot(hel; what2apply=log ∘ real)
dalitzplot(ch1; what2apply=log ∘ real)

# interesting_combinations
interesting_combinations = [
    (;j0=1,j=1,L=0,plus=false)
    (;j0=1,j=1,L=0,plus=true )
    (;j0=1,j=1,L=1,plus=false)
    (;j0=1,j=1,L=1,plus=true )
    (;j0=0,j=1,L=1,plus=false)
    (;j0=0,j=1,L=1,plus=true )
    (;j0=2,j=1,L=1,plus=false)
    (;j0=2,j=1,L=1,plus=true )
]

for ic in interesting_combinations
    @unpack j, j0, L, plus = ic
    sign = plus ? '+' : '-'
    dalitzplot(πDD_heli(;j0,j,L,plus);
        title=latexstring("X\\to (D^0\\pi^+)_\\mathrm{P}D^0\\textrm{$(sign)}D^0(\\pi^+D^0)_\\mathrm{P}"))
    annotate!(4.018, 4.038, text(latexstring("j_X=$(j0), L_{D^*D}=$(L)"), :left))
    savefig(joinpath("plots","nominal","dalitz","dalitz_j0=$(j0)_j=$(j)_L=$(L)_signplus=$(plus).pdf"))
end


ps = [let
    @unpack j, j0, L, plus = ic
    sign = plus ? '+' : '-'
    dalitzplot(πDD_heli(;j0,j,L,plus); ticks=false, xlab="", ylab="",
        title=latexstring("X\\to (D^0\\pi^+)_\\mathrm{P}D^0\\textrm{$(sign)}D^0(\\pi^+D^0)_\\mathrm{P}"))
    annotate!(4.018, 4.038, text(latexstring("j_X=$(j0), L_{D^*D}=$(L)"), :left))
end for ic in interesting_combinations]

plot(ps..., layout=grid(4,2), size=(450*2, 400*4))
savefig(joinpath("plots","nominal","dalitz","dalitz_summary.pdf"))
