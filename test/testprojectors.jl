using X2DDpi
using Plots
# 
theme(:wong2, minorticks=true, grid=false, frame=:box,
    guidefontvalign=:top, guidefonthalign=:right,
    foreground_color_legend = nothing,
    legendfontsize=9, legend =:topright,
    xlim=(:auto,:auto), ylim=(:auto,:auto), lab="")

const δm0_val = -359.2e-3 # MeV
const Γ0 = 48e-3 # MeV
# 
bw_e(s,e0,Γ) = J_I(s,
    BW_norm(m=e2m(e0), Γ=Γ*1e-3, n=sqrt(e2m(e0)*(Γ*1e-3) / π)))
# 
const srange = (e2m.(δm0_val .+ (-1, 1) .* 2Γ0)) .^ 2
const bwn = quadgk(s->abs2(bw_e(s,δm0_val,Γ0)), e2m(δm0_val-2Γ0)^2, e2m(δm0_val+2Γ0)^2)[1]


import X2DDpi: projectto3
projectto3(d, σ3) = projectto3(d, σ3,
    s->bw_e(s,δm0_val,Γ0),
    srange...)

π⁺D⁰D⁰ = πDD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), BW(m=mDˣ⁺,Γ=ΓDˣ⁺), BW(m=mDˣ⁺,Γ=ΓDˣ⁺))

@time projectto3(π⁺D⁰D⁰, e2m(δm0_val)^2, 2.007^2)
@time projectto3(π⁺D⁰D⁰, 2.007^2)

# should look similar
let d = π⁺D⁰D⁰, m = e2m(δm0_val)
    plot(e3->projectto3(d, m^2, e3^2),
        2.004, 2.011, lab="")
# 
    plot!(e3->projectto3(d, e3^2) / bwn,
        2.004, 2.011, lab="")
end



import X2DDpi: projectto1
projectto1(d, σ1) = projectto1(d, σ1,
    s->bw_e(s,δm0_val,Γ0),
    srange...)
#
# 
@time projectto1(π⁺D⁰D⁰, e2m(δm0_val)^2, 3.730^2)
@time projectto1(π⁺D⁰D⁰, 3.730^2)

# should look similar
let d = π⁺D⁰D⁰, m = e2m(δm0_val)
    plot(e1->projectto1(d, m^2, e1^2),
        3.729, 3.737, lab="")
# 
    plot!(e1->projectto1(d, e1^2) / bwn,
        3.729, 3.737, lab="")
end
