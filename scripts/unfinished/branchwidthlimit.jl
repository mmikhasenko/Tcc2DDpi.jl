using X2DDpi
using QuadGK
using FiniteDiff

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

underthrconst(Γ) = Γ*ρ_thr(channel(Γ), branch_position(Γ)-1e-7)
# to justify the factor Γ above 
let fΓ = 0.01
    plot(e->ρ_thr(channel(ΓDˣ⁺), e), -1, 1)
    plot!(e->fΓ*ρ_thr(channel(fΓ*ΓDˣ⁺), e), -1, 1)
end

let
    pv = -1:2
    factors = 1e3 .^ pv
    calv = underthrconst.(factors .* ΓDˣ⁺)
    plot(pv, real.(abs.(calv)), yscale=:log10, ylab=L"\mathrm{under\,\,threshold\,\,const}",
        xticks=(pv,[latexstring("10^{$(3i)}\\Gamma_{D^*}") for i in pv]))
end

const Eᵦ =branch_position(ΓDˣ⁺)
const approx_Eᵦ = Eᵦ-1e-7

function k3b(e; ϕ=-π)
    m = e2m(e)
    M = sqrt(mDˣ⁺^2-1im*mDˣ⁺*ΓDˣ⁺)
    p = cis(-ϕ/2)*sqrt((m-(M+mD⁰))*cis(ϕ))*
        sqrt(m+(M+mD⁰))*sqrt(m-(M-mD⁰))*sqrt(m+(M-mD⁰))/(2*m)
    return p
end
# simplified
function k(e; ϕ=-π)
    μ = mD⁰*mDˣ⁰ / (mD⁰+mDˣ⁰)
    Eᵦ = branch_position(ΓDˣ⁺)
    p = cis(-ϕ/2)*sqrt(2μ*(e-Eᵦ)/1e3*cis(ϕ))
    return p
end
k3b(-1e-5im)
k(-1e-5im)

heatmap(-1:0.1:1, -1:0.1:1, (x,y)->imag(k3b(Eᵦ+(1e3*ΓDˣ⁰)*(x+1im*y))))
# heatmap(-1:0.1:1, -1:0.1:1, (x,y)->imag(k(Eᵦ+(1e3*ΓDˣ⁰)*(x+1im*y))))



let f(e) = ρ_thr(channel(ΓDˣ⁺), e)
    # 
    ϕv = range(-7.4e-3, -6.9e-3, length=50)
    ev = Eᵦ .+ 0.001*(ΓDˣ⁺ * 1e3).*cis.(ϕv)
    # 
    calv = f.(ev)
    plot(1e3ϕv, [real.(calv) imag.(calv)], layout=grid(2,1))
    plot!(sp=1, title="rho3")
end
let f(e) = k3b(e; ϕ=-0.99774π)
    # 
    ϕv = range(-7.4e-3, -6.9e-3, length=30)
    ev = Eᵦ .+ 0.001*(ΓDˣ⁺ * 1e3).*cis.(ϕv)
    # 
    calv = f.(ev)
    plot(1e3ϕv, [real.(calv) imag.(calv)], layout=grid(2,1))
    plot!(sp=1, title="k3")
end

const ϕjump = -0.99774π


function Δk3b(r=1e-3,ϕ₀=-π-ϕjump,ϵ=1e-3)
    ev = Eᵦ .+ r*(1e3ΓDˣ⁺) .* cis.(ϕ₀ .+ ϵ.*[-1, 1])
    f(e) = k3b(e;ϕ=ϕjump)
    return diff(f.(ev))[1]
end
function Δρ_thr(r=1e-3,ϕ₀=-π-ϕjump,ϵ=1e-3)
    ev = Eᵦ .+ r*(1e3ΓDˣ⁺) .* cis.(ϕ₀ .+ ϵ.*[-1, 1])
    f(e) = ρ_thr(channel(ΓDˣ⁺), e)
    return diff(f.(ev))[1]
end

N = Δρ_thr(1e-4)/Δk3b(1e-4)
Σsum(e) = -1im*ρ_thr(channel(ΓDˣ⁺), e) + 1im*N*k3b(e;ϕ=ϕjump)


# check if the path is continues
ϕv = range(-π/2, 3π/2, length=600)
ev = Eᵦ .+ 0.1*(ΓDˣ⁺ * 1e3).*cis.(ϕv)

let f = Σsum
    calv = f.(ev)
    i = div(length(ev),2)
    Ni, Nr = real.(calv)[i], imag(calv)[i]
    plot(ϕv/π*180, [real.(calv)/Ni imag.(calv)/Nr])
end




Σ′ = FiniteDiff.finite_difference_derivative(
    e->ρ_thr(channel(ΓDˣ⁺), e),
    approx_Eᵦ)
# 
k3b′ = FiniteDiff.finite_difference_derivative(k3b, approx_Eᵦ)
# 
k′ = FiniteDiff.finite_difference_derivative(k, approx_Eᵦ)
