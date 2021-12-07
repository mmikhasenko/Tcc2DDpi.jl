using X2DDpi
using QuadGK

using Plots
using LaTeXStrings
theme(:wong2, frame=:box, grid=false, minorticks=true, 
    guidefontvalign=:top, guidefonthalign=:right,
    xlim=(:auto,:auto), ylim=(:auto,:auto), lw=1.2)


channel(Γ) = DˣD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), BW(m=mDˣ⁺; Γ))
# test
channel(ΓDˣ⁺)


function branch_position(Γ)
    ch = channel(Γ)
    σpole = pole_position(ch.R)
    mbranch = sqrt(σpole)+ch.ms[3]
    m2e(mbranch)
end
# test
@assert branch_position(ΓDˣ⁺) == m2e(sqrt(mDˣ⁺^2 - 1im * mDˣ⁺ * ΓDˣ⁺)+mD⁰)

# 
underthrconst(Γ) = Γ*ρ_thr(channel(Γ), branch_position(Γ)-1e-7)

let
    plot(e->ρ_thr(channel(ΓDˣ⁺), e), -1, 1)
    plot!(e->0.01*ρ_thr(channel(0.01ΓDˣ⁺), e), -1, 1)
end

let
    pv = -1:2
    factors = 1e3 .^ pv
    calv = underthrconst.(factors .* ΓDˣ⁺)
    plot(pv, real.(abs.(calv)), yscale=:log10, ylab=L"\mathrm{under\,\,threshold\,\,const}",
        xticks=(pv,[latexstring("10^{$(3i)}\\Gamma_{D^*}") for i in pv]))
end

using FiniteDiff

const Eᵦ =branch_position(ΓDˣ⁺)
const approx_Eᵦ = Eᵦ-1e-7

#

function k(e; ϕ=-π)
    μ = mD⁰*mDˣ⁰ / (mD⁰+mDˣ⁰)
    Eᵦ = branch_position(ΓDˣ⁺)
    cis(-ϕ/2)*sqrt(2μ*(e-Eᵦ)/1e3*cis(ϕ))
end
# 
function k3b(e; ϕ=-π)
    m = e2m(e)
    M = sqrt(mDˣ⁺^2-1im*mDˣ⁺*ΓDˣ⁺)
    cis(-ϕ/2)*sqrt((m-(M+mD⁰))*cis(ϕ))*sqrt(m+(M+mD⁰))*sqrt(m-(M-mD⁰))*sqrt(m+(M-mD⁰))/(2*m)
end

heatmap(-1:0.1:1, -1:0.1:1, (x,y)->imag(k3b(Eᵦ+(1e3*ΓDˣ⁰)*(x+1im*y))))
heatmap(-1:0.1:1, -1:0.1:1, (x,y)->imag(k(Eᵦ+(1e3*ΓDˣ⁰)*(x+1im*y))))


k3b(-1e-5im)
k(-1e-5im)


Σ′ = FiniteDiff.finite_difference_derivative(
    e->ρ_thr(channel(ΓDˣ⁺), e),
    approx_Eᵦ)
# 
k3b′ = FiniteDiff.finite_difference_derivative(
    k3b,
    approx_Eᵦ)
# 
k′ = FiniteDiff.finite_difference_derivative(
    k,
    approx_Eᵦ)
#
N = Σ′/k3b′
# 
Σsum(e) = -1im*ρ_thr(channel(ΓDˣ⁺), e) + 1im*N*k3b(e;ϕ=-0.999π)
#

# 
ϕv = range(-π/2, 3π/2, length=600)
ev = Eᵦ .+ 0.001*(ΓDˣ⁺ * 1e3).*cis.(ϕv)

# let f(e) = k3b(e)*N
#     calv = f.(ev)
#     plot(ϕv/π*180, [real.(calv) imag.(calv)])
# end
# let f(e) = ρ_thr(channel(ΓDˣ⁺), e)
#     calv = f.(ev)
#     plot!(ϕv/π*180, [real.(calv) imag.(calv)])
# end

let f = Σsum
    calv = f.(ev)
    i = div(length(ev),2)
    Ni, Nr = real.(calv)[i], imag(calv)[i]
    @show Ni, Nr
    plot(ϕv/π*180, [real.(calv)/Ni imag.(calv)/Nr])
end

