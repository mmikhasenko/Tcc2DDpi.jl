using X2DDpi

using Plots
using LaTeXStrings
theme(:wong2, frame=:box, grid=false, minorticks=true, 
    guidefontvalign=:top, guidefonthalign=:right,
    xlim=(:auto,:auto), ylim=(:auto,:auto), lw=1.2)

const ch₀ = DˣD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), BW(m=mDˣ⁺, Γ=ΓDˣ⁺))
ρ₀(e) = ρ_thr(ch₀, e)

let
    plot(ylab=L"\mathrm{Real\,\,part}", xlab=L"e / \Gamma^{*+}")
    plot!(-1:0.09:1, x->real(k3b((1e3*ΓDˣ⁰)*x)/real(k3b(0.001))))
    plot!(-1:0.09:1, x->real(kNR((1e3*ΓDˣ⁰)*x)/real(kNR(0.001))))
    plot!(-1:0.09:1, x->real(ρ₀( (1e3*ΓDˣ⁰)*x)/real(ρ₀( 0.001))))
end

let
    plot(ylab=L"\mathrm{Imaginary\,\,part}", xlab=L"e\,\,[\mathrm{MeV}]")
    plot!(-0.1:0.009:0.1, x->imag(k3b(x)/real(k3b(0.001))))
    plot!(-0.1:0.009:0.1, x->imag(kNR(x)/real(kNR(0.001))))
    plot!(-0.1:0.009:0.1, x->imag(ρ₀( x)/real(ρ₀( 0.001))))
end

let
    plot(size=(1000,300), layout=grid(1,3),
        ylab=L"\mathrm{Im}\,e / \Gamma^{*+}",
        xlab=L"\mathrm{Re}\,e / \Gamma^{*+}",
        colorbar=false)
    heatmap!(sp=1, -1:0.1:1, -1:0.1:1,   (x,y)->real(k3b(Eᵦˣ⁺ +(1e3*ΓDˣ⁰)*(x+1im*y))))
    heatmap!(sp=2, -1:0.1:1, -1:0.1:1,   (x,y)->real(kNR(Eᵦˣ⁺ +(1e3*ΓDˣ⁰)*(x+1im*y))))
    heatmap!(sp=3, -1:0.09:1, -1:0.09:1, (x,y)->real(ρ₀( Eᵦˣ⁺ +(1e3*ΓDˣ⁰)*(x+1im*y))))
end
    

