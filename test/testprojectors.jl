using X2DDpi


function bw_e(s,e0,Γ)
    m = e2m(e0)
    return sqrt(m*(Γ*1e-3) / π) / (m^2-s-1im*m*(Γ*1e-3))
end
# 
bw_e(s,e0,Γ) = J_I(s, BW(m=e2m(e0), Γ=Γ, c=sqrt(m*(Γ*1e-3) / π)))


import X2DDpi: projectto3
projectto3(d, σ3) = projectto3(d, σ3,
    s->bw_e(s,δm0_val,Γ0),
    e2m(δm0_val-2Γ0)^2, e2m(δm0_val+2Γ0)^2)

const δm0_val = -359.2e-3
π⁺D⁰D⁰ = πDD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), (m=mDˣ⁺,Γ=ΓDˣ⁺), (m=mDˣ⁺,Γ=ΓDˣ⁺))

@time projectto3(π⁺D⁰D⁰, e2m(δm0_val)^2, 2.007^2)
@time projectto3(π⁺D⁰D⁰, 2.007^2)

# should look similar
let d = π⁺D⁰D⁰, m = e2m(δm0_val)
    plot(e3->projectto3(d, m^2, e3^2),
        2.004, 2.011, lab="")
# 
    plot!(e3->projectto3(d, e3^2)*1.19,
        2.004, 2.011, lab="")
end
