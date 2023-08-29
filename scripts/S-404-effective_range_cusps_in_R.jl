using X2DDpi
using Parameters

using Plots
using LaTeXStrings
theme(:wong2, frame=:box, grid=false, minorticks=true, 
    guidefontvalign=:top, guidefonthalign=:right,
    xlim=(:auto,:auto), ylim=(:auto,:auto), lw=1.2)


settings = transformdictrecursively!(readjson("settings.json"), ifstringgivemeasurement)
@unpack δm0 = settings["fitresults"]
const δm0_val = δm0.val
@unpack cutoff, estep = settings["phspmatching"]

ere(k; a⁻¹, r, N) = N*(a⁻¹ + r * k^2 / 2 - 1im* k)
ere(k, p) = ere(k, a⁻¹ = p.a⁻¹, r = p.r, N = p.N)


# DˣD
const A₀_DˣD = Amplitude(
    interpolated(
        DˣD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), BW(m=mDˣ⁺, Γ=ΓDˣ⁺)),
        cutoff; estep=0.2*estep))
# 
ERE_DˣD = let
    @time effrangepars = # 12s
        effectiverangeexpansion(
            Δe->denominator_II(A₀_DˣD, Eᵦˣ⁺+Δe, δm0_val),
            Δe->k3b(Eᵦˣ⁺+Δe),
            ComplexBranchPointExpansion(CircularSum(abs(imag(Eᵦˣ⁺))/10, 150)))
    # 
    (; tophysicsunits(effrangepars)..., effrangepars...)
end

let
    plot(size=(700,300), title=["Real" "Imag"], yticks=false,
        plot(Δe->real(denominator_II(A₀_DˣD, Δe+Eᵦˣ⁺, δm0_val)-ere(k3b(Δe+Eᵦˣ⁺); a⁻¹=0, r=0, N=ERE_DˣD.N)), -1.0, 1.0, lab="Model"),
        plot(Δe->imag(denominator_II(A₀_DˣD, Δe+Eᵦˣ⁺, δm0_val)-ere(k3b(Δe+Eᵦˣ⁺); a⁻¹=0, r=0, N=ERE_DˣD.N)), -1.0, 1.0, lab="Model"))
    #
    plot!(sp=1, Δe->real(ere(k3b(Δe+Eᵦˣ⁺),ERE_DˣD)-ere(k3b(Δe+Eᵦˣ⁺); a⁻¹=0, r=0, N=ERE_DˣD.N)), -1.0, 1.0, lab="Eff.Range")
    plot!(sp=2, Δe->imag(ere(k3b(Δe+Eᵦˣ⁺),ERE_DˣD)-ere(k3b(Δe+Eᵦˣ⁺); a⁻¹=0, r=0, N=ERE_DˣD.N)), -1.0, 1.0, lab="Eff.Range")
    # 
    scatter!(sp=1, [0], [real(ere(k3b(Eᵦˣ⁺), ERE_DˣD))], lab="", m=(:red,5))
    scatter!(sp=2, [0], [imag(ere(k3b(Eᵦˣ⁺), ERE_DˣD))], lab="", m=(:red,5))
end



# DˣD
const A₀_DˣD_course = Amplitude(
    interpolated(
        DˣD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), BW(m=mDˣ⁺, Γ=ΓDˣ⁺)),
        cutoff; estep=estep))
# 
ERE_DˣD_course = let
    @time effrangepars = # 12s
        effectiverangeexpansion(
            Δe->denominator_II(A₀_DˣD_course, Eᵦˣ⁺+Δe, δm0_val),
            Δe->k3b(Eᵦˣ⁺+Δe),
            ComplexBranchPointExpansion(CircularSum(abs(imag(Eᵦˣ⁺))/10, 150)))
    # 
    (; tophysicsunits(effrangepars)..., effrangepars...)
end

let
    plot!(sp=1, Δe->real(denominator_II(A₀_DˣD_course, Δe+Eᵦˣ⁺, δm0_val)-ere(k3b(Δe+Eᵦˣ⁺); a⁻¹=0, r=0, N=ERE_DˣD_course.N)), -1.0, 1.0, lab="", c=1, ls=:dash)
    plot!(sp=2, Δe->imag(denominator_II(A₀_DˣD_course, Δe+Eᵦˣ⁺, δm0_val)-ere(k3b(Δe+Eᵦˣ⁺); a⁻¹=0, r=0, N=ERE_DˣD_course.N)), -1.0, 1.0, lab="", c=1, ls=:dash)
    #
    plot!(sp=1, Δe->real(ere(k3b(Δe+Eᵦˣ⁺),ERE_DˣD_course)-ere(k3b(Δe+Eᵦˣ⁺); a⁻¹=0, r=0, N=ERE_DˣD_course.N)), -1.0, 1.0, lab="", c=2, ls=:dash)
    plot!(sp=2, Δe->imag(ere(k3b(Δe+Eᵦˣ⁺),ERE_DˣD_course)-ere(k3b(Δe+Eᵦˣ⁺); a⁻¹=0, r=0, N=ERE_DˣD_course.N)), -1.0, 1.0, lab="", c=2, ls=:dash)
end