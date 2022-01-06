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

const Eᵦ = branch_position(ΓDˣ⁺)
# const approx_Eᵦ = Eᵦ-1e-7

function k3b(e)
    m = e2m(e)
    M = sqrt(mDˣ⁺^2-1im*mDˣ⁺*ΓDˣ⁺)
    p = cis(π/4)*sqrt((m-(M+mD⁰))*cis(-π/2))*
        sqrt(m+(M+mD⁰))*sqrt(m-(M-mD⁰))*sqrt(m+(M-mD⁰))/(2*m)
    return p
end

# simplified
function k(e)
    μ = mD⁰*mDˣ⁰ / (mD⁰+mDˣ⁰)
    Eᵦ = branch_position(ΓDˣ⁺)
    p = cis(π/4)*sqrt(2μ*(e-Eᵦ)/1e3*cis(-π/2))
    return p
end

heatmap(-1:0.1:1, -1:0.1:1, (x,y)->imag(k3b(Eᵦ+(1e3*ΓDˣ⁰)*(x+1im*y))))
heatmap(-1:0.1:1, -1:0.1:1, (x,y)->imag(k(Eᵦ+(1e3*ΓDˣ⁰)*(x+1im*y))))


decay_matrix_element_squared(channel(ΓDˣ⁺),e2m(branch_position(ΓDˣ⁺))^2,σ3,0.0)

ρ_thr(channel(ΓDˣ⁺), 0.01-0.01im)
ρ_thr(channel(ΓDˣ⁺), 0.01-0.01im)
ρ_thr(channel(ΓDˣ⁺), 0.01-0.01im)

heatmap(-1:0.09:1, -1:0.09:1, (x,y)->real(ρ_thr(channel(ΓDˣ⁺), (0.0+(1e3*ΓDˣ⁰)*(x+1im*y)))))



function Δk3b(r=1e-3,ϵ=1e-3)
    ev = Eᵦ .+ r*(1e3ΓDˣ⁺) .* cis.(-π/2 .+ ϵ.*[-1, 1])
    return diff(k3b.(ev))[1]
end
function Δρ_thr(r=1e-3,ϵ=1e-3)
    ev = Eᵦ .+ r*(1e3ΓDˣ⁺) .* cis.(-π/2 .+ ϵ.*[-1, 1])
    f(e) = ρ_thr(channel(ΓDˣ⁺), e)
    return diff(f.(ev))[1]
end

N = Δρ_thr(1e-4)/Δk3b(1e-4)
Σsum(e) = -1im*ρ_thr(channel(ΓDˣ⁺), e) + 1im*N*k3b(e)



# check if the path is continues
ϕv = range(-π, π, length=600)
ev = Eᵦ .+ 0.1*(ΓDˣ⁺ * 1e3).*cis.(ϕv)

let f = Σsum
    calv = f.(ev)
    i = div(length(ev),2)
    Ni, Nr = real.(calv)[i], imag(calv)[i]
    plot(ϕv/π*180, [real.(calv)/Ni imag.(calv)/Nr])
end


ρ_thr(1e-4)

ρ_thr(channel(ΓDˣ⁺), Eᵦ)


ρ_thr.(Ref(channel(ΓDˣ⁺)), Eᵦ + abs(imag(Eᵦ))/2 * cis.(-π/2)-1e-12)

ρ_thr.(Ref(channel(ΓDˣ⁺)), Eᵦ .+ abs(imag(Eᵦ))/2 .* cis.(-π:0.1:π))

c1 = cauchy(e->ρ_thr(channel(ΓDˣ⁺), e), Eᵦ, abs(imag(Eᵦ))/2000)
c2 = cauchy(k3b, Eᵦ, abs(imag(Eᵦ))/2000)

c1/c2
N

ρ_thr(channel(ΓDˣ⁺), Eᵦ-0.1im)
ρ_thr(channel(ΓDˣ⁺), Eᵦ-0.01im)

let 
    d = channel(ΓDˣ⁺)
    s = e2m(Eᵦ-0.1im)^2
    # 
    intgr(x) = X2DDpi.integrand_mapped_thr(d,s,x)
    # quadgk(intgr, 0, 1)
    # intgr(0.6471418489543712)
    xF = range(0.6471418489543694, 0.647141848954373, length=10)[3]
    # intgr(xF)
    # 
    σ3 = 4.04112838143364 - 0.00016765533372im
    # X2DDpi.decay_matrix_element_squared(d,s,σ3,1.0)
    J_I(σ3,d.R)
end

1/(0.0)
1im/(0.0)
1/(0.0+0.0im)

Inf + Inf*im