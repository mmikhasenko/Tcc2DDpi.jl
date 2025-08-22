using Pkg
cd(joinpath(@__DIR__, ".."))
Pkg.activate(".")

using X2DDpi

using Plots
using LaTeXStrings
theme(:wong2, size = (450, 400), minorticks = true, grid = false, frame = :box,
    guidefontvalign = :top, guidefonthalign = :right,
    foreground_color_legend = nothing,
    legendfontsize = 9, legend = :topright,
    xlim = (:auto, :auto), ylim = (:auto, :auto),
    xlab = L"m_{D^0\pi^+}^2\,\,[\mathrm{GeV}]", ylab = L"m_{D^0\pi^+}^2\,\,[\mathrm{GeV}]")


# 
using Parameters
using Measurements
# 
using ThreeBodyDecay
using PartialWaveFunctions


###################################################################

@userplot Dalitz
@recipe function f(hp::Dalitz; what2apply = real)
    model, = hp.args
    iσx := 2
    iσy := 3
    density = σs -> what2apply(X2DDpi.decay_matrix_element_squared(model, s0, σs.σ3, σs.σ2))
    (ms, density)
end

###################################################################

settings = transformdictrecursively!(readjson("settings.json"), ifstringgivemeasurement)
#
@unpack δm0 = settings["fitresults"]

# const m0 = e2m(δm0.val) 
# const s0 = e2m(δm0.val)^2 
# const ms = ThreeBodyMasses(ch1.ms...;m0)

###################################################################

abstract type AbstractLineshape end

@with_kw struct BW <: AbstractLineshape
    m::Float64
    Γ::Float64
end

import X2DDpi: J_I, J_II

function J_I(pars::BW, σ)
    @unpack m, Γ = pars
    FF = 1.0 # (σ-(mD⁰+mπ⁺)^2) / (m^2-(mD⁰+mπ⁺)^2)
    1 / (m^2 - σ + 1im * m * Γ * FF)
end

function J_II(pars::BW, σ)
    @unpack m, Γ = pars
    FF = 1.0 # (σ-(mD⁰+mπ⁺)^2) / (m^2-(mD⁰+mπ⁺)^2)
    1 / (m^2 - σ - 1im * m * Γ * FF)
end


abstract type AbstractDecayChain end

@with_kw struct decaychain{T} <: AbstractDecayChain
    k::Int
    me::T
    j::Int
    L::Int
    relativeplus::Bool = true
end

index(dc::decaychain) = dc.k

# 
decaychain_3(v...; kw...) = decaychain(v...; k = 3, kw...)
decaychain_2(v...; kw...) = decaychain(v...; k = 2, kw...)
# 
decaychain_12(v...; kw...) = decaychain(v...; k = 3, kw...)
decaychain_31(v...; kw...) = decaychain(v...; k = 2, kw...)
# 
decaychain_3(me = BW(m = mDˣ⁺, Γ = ΓDˣ⁺), j = 1, L = 0)


import X2DDpi: AbstractxDD


@with_kw struct decaymodel <: AbstractxDD
    ms::NamedTuple{(:m1, :m2, :m3), NTuple{3, Float64}}
    Rs::Vector{decaychain{T} where T <: AbstractLineshape}
    #
    j0::Int = 1
end

ch1 = decaymodel(
    j0 = 1,
    ms = (m1 = mπ⁺, m2 = mD⁰, m3 = mD⁰),
    Rs = [
        decaychain_12(me = BW(m = mDˣ⁺, Γ = ΓDˣ⁺), j = 1, L = 0),
        decaychain_31(me = BW(m = mDˣ⁺, Γ = ΓDˣ⁺), j = 1, L = 0)]) # reference frame for 0 helicity
# 
ch2 = decaymodel(
    j0 = 1,
    ms = (m1 = mπ⁰, m2 = mD⁺, m3 = mD⁰),
    Rs = [
        decaychain_12(me = BW(m = mDˣ⁺, Γ = ΓDˣ⁺), j = 1, L = 0),
        decaychain_31(me = BW(m = mDˣ⁰, Γ = ΓDˣ⁰), j = 1, L = 0)]) # reference frame for 0 helicity
# 
# ch3 = decaymodel(
#     j0 = 1,
#     ms = (m1=mγ, m2=mD⁺,m3=mD⁰),
#     Rs = [
#         decaychain_12(me = BW(m=mDˣ⁺,Γ=ΓDˣ⁺), j=1, L=0)
#         decaychain_31(me = BW(m=mDˣ⁰,Γ=ΓDˣ⁰), j=1, L=0)],
#     ref0 = 1) # reference frame for 0 helicity
# # 



# const ch1 = πDD(
#     (m1=mπ⁺,m2=mD⁰,m3=mD⁰),  # massed and order
#     (12, BW(m=mDˣ⁺,Γ=ΓDˣ⁺)), # channel 12 with D*
#     (31, BW(m=mDˣ⁺,Γ=ΓDˣ⁺))) # channel 31 with D*
# 

function pk(k, σs, msq)
    i, j, _ = ijk(k)
    sqrt(ThreeBodyDecays.λ(σs[k], msq[i], msq[j]) / (4σs[k]))
end
qk(k, σs, msq) = sqrt(ThreeBodyDecays.λ(msq[4], σs[k], msq[k]) / (4msq[4]))


using ThreeBodyDecays
import X2DDpi: decay_matrix_element_squared
#

function decay_matrix_element_squared(d::decaymodel, s, σ3, σ2)
    @unpack ms, j0, Rs = d
    #
    H(ν, j, L) = clebschgordan(L, 0, j, ν, j0, ν)
    # 
    ms4 = ThreeBodyMasses(m0 = √s; ms...)
    msq4 = ms4^2
    σs = Invariants(ms4; σ3, σ2)
    # 
    𝔐² = 0.0im
    for (a, Ra) in enumerate(Rs),
        (b, Rb) in enumerate(Rs)
        # 
        a < b && continue
        # 
        ka = index(Ra)
        σa = σs[ka]
        #
        Ja_I, Ja_II = J_I(Ra.me, σa), J_II(Ra.me, σa)
        #
        pʲa = (Ra.j == 0 ? 1.0 : pk(ka, σs, msq4)^Ra.j)
        qᴸa = (Ra.L == 0 ? 1.0 : qk(ka, σs, msq4)^Ra.L)
        #
        cosθa = ThreeBodyDecays.cosθ(ka, σs, msq4)
        # 
        if a == b
            pʲ, qᴸ = pʲa, qᴸa
            @unpack j, L = Ra
            #
            𝔐² += sum(abs2,
                      wignerd(j, λ, 0, za) * H(λ, j, L) * pʲ * qᴸ
                      for λ in -j:j) * Ja_I * Ja_II +
                  continue
        end
        #
        Jb_I, Jb_II = J_I(Rb.me, σk), J_II(Rb.me, σk)
        #
        pʲb = (Rb.j == 0 ? 1.0 : pk(kb, σs, msq4)^Rb.j)
        qᴸb = (Rb.L == 0 ? 1.0 : qk(kb, σs, msq4)^Rb.L)
        # 
        w0 = wr(Ra.k, Rb.k, 0)
        cosζ0 = cosζ(w0, σs, msq4)
        cosθb = cosθ(Rb.k, σs, msq4)
        #
        𝔐² +=
            Ra.relativeplus * Rb.relativeplus *
            # 
            (Jb_I * Ja_II + Ja_I * Jb_II) *
            # 
            sum(wignerd(j, λ, 0, cosθa) * H(λ, Ra.j, Ra.L) * pʲa * qᴸa *
                wignerd(j, ν, 0, cosθb) * H(ν, Rb.j, Rb.L) * pʲb * qᴸb *
                wignerd_sign(j0, λ, ν, cosζ0, ispositive(w0))
            # 
                for (ν, λ) in Iterators.product(-Rb.j:Rb.j, -Ra.j:Ra.j))
    end
    #
    return 𝔐²
end

J_I
# (ch1.Rs[1].me, 2.1)
J_II(ch1.Rs[1].me, 2.1)

# ms=
# hel = πDD_heli(;plus=false)
let
    s0 = e2m(-0.3)^2
    # 
    ms4 = ThreeBodyMasses(m0 = √s0; ch1.ms...)
    msq4 = ms4^2
    # 
    σs = randomPoint(ms4)
    # decay_matrix_element_squared(hel,s0, σs.σ3,σs.σ2)
    decay_matrix_element_squared(ch1, s0, σs.σ3, σs.σ2)
end

dalitzplot(hel; what2apply = real)
dalitzplot(hel; what2apply = log ∘ real)
dalitzplot(ch1; what2apply = log ∘ real)

# interesting_combinations
interesting_combinations = [
    (; j0 = 1, j = 1, L = 0, plus = false);
    (; j0 = 1, j = 1, L = 0, plus = true);
    (; j0 = 1, j = 1, L = 1, plus = false);
    (; j0 = 1, j = 1, L = 1, plus = true);
    (; j0 = 0, j = 1, L = 1, plus = false);
    (; j0 = 0, j = 1, L = 1, plus = true);
    (; j0 = 2, j = 1, L = 1, plus = false);
    (; j0 = 2, j = 1, L = 1, plus = true)
]

for ic in interesting_combinations
    @unpack j, j0, L, plus = ic
    sign = plus ? '+' : '-'
    dalitzplot(πDD_heli(; j0, j, L, plus);
        title = latexstring("X\\to (D^0\\pi^+)_\\mathrm{P}D^0\\textrm{$(sign)}D^0(\\pi^+D^0)_\\mathrm{P}"))
    annotate!(4.018, 4.038, text(latexstring("j_X=$(j0), L_{D^*D}=$(L)"), :left))
    savefig(joinpath("plots", "nominal", "dalitz", "dalitz_j0=$(j0)_j=$(j)_L=$(L)_signplus=$(plus).pdf"))
end


ps = [
    let
        @unpack j, j0, L, plus = ic
        sign = plus ? '+' : '-'
        dalitzplot(πDD_heli(; j0, j, L, plus); ticks = false, xlab = "", ylab = "",
            title = latexstring("X\\to (D^0\\pi^+)_\\mathrm{P}D^0\\textrm{$(sign)}D^0(\\pi^+D^0)_\\mathrm{P}"))
        annotate!(4.018, 4.038, text(latexstring("j_X=$(j0), L_{D^*D}=$(L)"), :left))
    end for ic in interesting_combinations
]

plot(ps..., layout = grid(4, 2), size = (450 * 2, 400 * 4))
savefig(joinpath("plots", "nominal", "dalitz", "dalitz_summary.pdf"))
