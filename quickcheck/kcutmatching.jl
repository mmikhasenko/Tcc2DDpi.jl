using X2DDpi
using Parameters


const Eᵦ = X2DDpi.Eᵦˣ⁺

# matching discontinuity
function Δk3b(r=1e-3,ϵ=1e-3)
    ev = Eᵦ .+ r*(1e3ΓDˣ⁺) .* cis.(-π/2 .+ ϵ.*[-1, 1])
    return diff(k3b.(ev))[1]
end
function Δρ_thr(r=1e-3,ϵ=1e-3)
    ev = Eᵦ .+ r*(1e3ΓDˣ⁺) .* cis.(-π/2 .+ ϵ.*[-1, 1])
    f(e) = ρ_thr(channel(ΓDˣ⁺), e)
    return diff(f.(ev))[1]
end

N = Δρ_thr(1e-4)/Δk3b(1e-4) # 0.003475797316511961 - 3.7679086257784454e-5im
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

let
    c1 = circleintegral(e->ρ_thr(channel(ΓDˣ⁺), Eᵦ+e), abs(imag(Eᵦ))/20)
    c2 = circleintegral(e->k3b(Eᵦ+e), abs(imag(Eᵦ))/20)
    c1/c2
end # 0.003475797316511961 - 3.7679086257784454e-5im
N   # 0.0034757971014234288 - 3.766391144490617e-5im

