using X2DDpi
using QuadGK

using Plots
using LaTeXStrings
theme(:wong2, frame=:box, grid=false, minorticks=true, 
    guidefontvalign=:top, guidefonthalign=:right,
    xlim=(:auto,:auto), ylim=(:auto,:auto), lw=1.2)


channel(Γ) = DˣD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), BW(m=mDˣ⁺; Γ))

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
