### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 5ac3c9c0-0ce8-11ed-102b-03e4d63d6c1b
# ╠═╡ show_logs = false
begin
	cd(joinpath(@__DIR__, ".."))
	import Pkg
	Pkg.activate(".")
	Pkg.instantiate()
	# 
	using X2DDpi
	using Parameters	
	using Plots
	using LaTeXStrings
end

# ╔═╡ d095fe63-50e0-4069-ab16-9bb66f63839e
md"""
# Effective range: cusps in $R^*(s)$

The effective range expansion reads

$A^{-1}(s) = N^*\,(R^*(s) - ik+O(a^3 k^4))$

where $R(s)$ is supposed to be regular at the expansion point, $k=0$.
When looking carefully, at the $R$, one finds irregular cusps near the expantion point. This notebook shows that the origin of the irregularities is the discrete grid for the dispersion integral. 

When decrising the step side, the function gets smooth.
"""

# ╔═╡ 13804d6e-fee0-47a2-a0f2-4883ce14bcb6
theme(:wong2, size=(500, 350), minorticks=true, grid=false, frame=:box,
    guidefontvalign=:top, guidefonthalign=:right,
    foreground_color_legend=nothing,
    legendfontsize=9, legend=:topright,
    xlim=(:auto, :auto), ylim=(:auto, :auto))

# ╔═╡ 04808f8e-25bf-4c6c-9794-5a9d59b197e8
begin
	settings = transformdictrecursively!(readjson("settings.json"), 
		ifstringgivemeasurement)
	#
	@unpack δm0 = settings["fitresults"]
	const δm0_val = δm0.val
	@unpack cutoff, estep = settings["phspmatching"]
end

# ╔═╡ 3fdd51a5-8cdf-4716-81ac-da3cbe37ed76
begin
	ere(k; a⁻¹, r, N) = N*(a⁻¹ + r * k^2 / 2 - 1im* k)
	ere(k, p) = ere(k, a⁻¹ = p.a⁻¹, r = p.r, N = p.N)
end

# ╔═╡ 8ff3f33b-eda9-413e-bef8-a259ee5e3345
md"""
## Construction of models
"""

# ╔═╡ 654f8c7c-c14c-4934-b876-8fab5a2575f3
const model = Amplitude(
    interpolated(
        DˣD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), BW(m=mDˣ⁺, Γ=ΓDˣ⁺)),
        cutoff; estep=0.2*estep));

# ╔═╡ 8f1cba70-78c6-4e36-96f8-490ea7cabac7
const model_course = Amplitude(
    interpolated(
        DˣD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), BW(m=mDˣ⁺, Γ=ΓDˣ⁺)),
        cutoff; estep=estep));

# ╔═╡ 549a4df0-75e5-4ab3-bad4-7ce0f173323a
md"""
## Effective range
"""

# ╔═╡ d888b34c-1e00-407f-85bf-a48a8ef354d0
ERE = let
    effrangepars = # 12s
        effectiverangeexpansion(
            Δe->denominator_II(model, Eᵦˣ⁺+Δe, δm0_val),
            Δe->k3b(Eᵦˣ⁺+Δe),
            ComplexBranchPointExpansion(CircularSum(abs(imag(Eᵦˣ⁺))/10, 150)))
    # 
    (; tophysicsunits(effrangepars)..., effrangepars...)
end

# ╔═╡ 680aff95-d592-4a38-996c-d94d5c007858
ERE_course = let
    effrangepars = # 12s
        effectiverangeexpansion(
            Δe->denominator_II(model_course, Eᵦˣ⁺+Δe, δm0_val),
            Δe->k3b(Eᵦˣ⁺+Δe),
            ComplexBranchPointExpansion(CircularSum(abs(imag(Eᵦˣ⁺))/10, 150)))
    # 
    (; tophysicsunits(effrangepars)..., effrangepars...)
end

# ╔═╡ 01b6fb58-8f80-4adc-b5df-85baa2e63911
let elim = (-0.5, 0.5)
	D1(Δe) = denominator_II(model, Δe+Eᵦˣ⁺, δm0_val)
	minus_Nik1(Δe) = ere(k3b(Δe+Eᵦˣ⁺); a⁻¹=0, r=0, N=ERE.N)
	R1(Δe) = D1(Δe) - minus_Nik1(Δe)
	# 
    plot(size=(700,300), title=["Real" "Imag"], yticks=false,
        plot(real ∘ R1, elim..., lab="Model"),
        plot(imag ∘ R1,elim..., lab="Model"))
    #
	ERE1(Δe) = ere(k3b(Δe+Eᵦˣ⁺), ERE)
	RERE1(Δe) = ERE1(Δe) - minus_Nik1(Δe)
	# 
    plot!(sp=1, real ∘ RERE1, elim..., lab="Eff.Range", leg=:topleft)
    plot!(sp=2, imag ∘ RERE1, elim..., lab="")
    #
	D2(Δe) = denominator_II(model_course, Δe+Eᵦˣ⁺, δm0_val)
	minus_Nik2(Δe) = ere(k3b(Δe+Eᵦˣ⁺); a⁻¹=0, r=0, N=ERE_course.N)
	R2(Δe) = D2(Δe) - minus_Nik2(Δe)	
	# 
	plot!(sp=1, real ∘ R2, elim..., lab="", c=1, ls=:dash)
    plot!(sp=2, imag ∘ R2, elim..., lab="Model(course grid)", c=1, ls=:dash)
    #
	ERE2(Δe) = ere(k3b(Δe+Eᵦˣ⁺), ERE_course)
	RERE2(Δe) = ERE2(Δe) - minus_Nik2(Δe)
	# 
    plot!(sp=1, real ∘ RERE2, elim..., lab="", c=2, ls=:dash)
    plot!(sp=2, imag ∘ RERE2, elim..., lab="", c=2, ls=:dash)
	# 
    scatter!(sp=1, [0], [real(ere(k3b(Eᵦˣ⁺), ERE))], lab="", m=(:red,3))
    scatter!(sp=2, [0], [imag(ere(k3b(Eᵦˣ⁺), ERE))], lab="", m=(:red,3))
end

# ╔═╡ Cell order:
# ╟─d095fe63-50e0-4069-ab16-9bb66f63839e
# ╠═5ac3c9c0-0ce8-11ed-102b-03e4d63d6c1b
# ╠═13804d6e-fee0-47a2-a0f2-4883ce14bcb6
# ╠═04808f8e-25bf-4c6c-9794-5a9d59b197e8
# ╠═3fdd51a5-8cdf-4716-81ac-da3cbe37ed76
# ╟─8ff3f33b-eda9-413e-bef8-a259ee5e3345
# ╠═654f8c7c-c14c-4934-b876-8fab5a2575f3
# ╠═8f1cba70-78c6-4e36-96f8-490ea7cabac7
# ╟─549a4df0-75e5-4ab3-bad4-7ce0f173323a
# ╠═d888b34c-1e00-407f-85bf-a48a8ef354d0
# ╠═680aff95-d592-4a38-996c-d94d5c007858
# ╠═01b6fb58-8f80-4adc-b5df-85baa2e63911
