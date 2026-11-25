### A Pluto.jl notebook ###
# v0.16.1

using Markdown
using InteractiveUtils

# ╔═╡ 61b01d0b-1f8f-4e9d-805a-b533cd39410c
begin
    import Pkg
    Pkg.activate(@__DIR__)
end

# ╔═╡ 8a8ef3e6-1277-407c-b906-7fb173230a7d
begin
	using Measurements
	using Unitful
end

# ╔═╡ 74db3011-7c88-40dc-ae78-bd9c3e6be856
md"""
# Effective range for $\chi_{c1}(3872)$
"""

# ╔═╡ 033a5372-bdd0-4065-8147-7a9a375ed3d4
begin
	mX = -7.18u"MeV"
	gLHCb = 0.108
end;

# ╔═╡ a32cb39c-5ad0-457e-8e78-6dc1afe795b8
begin
	const mD⁰ = 1.86483u"GeV"
	const mD⁺ = 1.86965u"GeV"
	# 
	const mDˣ⁺ = (1.86483+145.4258e-3)*1u"GeV" # m(D) + Δm(D*,D) from PDG
	const mDˣ⁰ = 2.00685u"GeV"
	# 
	const fm_times_mev = 197.3269804u"fm*MeV"
end

# ╔═╡ 1c80959c-ed0a-44e9-8da1-b16e1d0e3bc0
begin
	μ = mD⁰*mDˣ⁰ / (mD⁰+mDˣ⁰)
	μ⁺ = mD⁺*mDˣ⁺ / (mD⁺+mDˣ⁺)
	δ = (mD⁺ + mDˣ⁺)-(mDˣ⁰ + mD⁰)
	(;μ, μ⁺, δ)
end

# ╔═╡ 1eb7633e-4c75-4d18-96d5-bf2fadf8e5cd
md"""
## Scattering length
"""

# ╔═╡ cc0af199-3b96-4f8b-a31a-706bdf716365
inva0_1, inva0_2 = 2mX/gLHCb, sqrt(2μ⁺*δ)

# ╔═╡ 5f9294d6-68a6-4cde-b59d-4a9790d63668
inva0 = uconvert(u"MeV", inva0_1 + inva0_2)

# ╔═╡ e73dbf75-53d8-4ec4-b23e-3a15abb3a28b
md"""
## Effective range
"""

# ╔═╡ 7f7af133-c60c-4317-b357-9840eebb4785
md"""
nominal values
"""

# ╔═╡ 0e275bfc-9d93-4d02-957b-de9fc77c1964
r0_1, r0_2 =
	-2/(μ*gLHCb),
	-sqrt(μ⁺/(2*μ^2*δ))

# ╔═╡ 782412e1-b150-459a-b973-3e7e27c20408
r0 = uconvert(u"fm", (r0_1 + r0_2) * fm_times_mev)

# ╔═╡ 19740af1-19cc-482f-83fe-613ada9a55f8
md"""
as a function of g
"""

# ╔═╡ 43bcd184-ffe7-40da-832c-025028649272
r(g) = uconvert(u"fm", 
	(-2/(μ*g) - sqrt(μ⁺/(2*μ^2*δ))) *
		fm_times_mev)

# ╔═╡ 89d036a1-1c15-40ee-890a-cecfe8857297
invg_66CL = (0.27, 33.3)

# ╔═╡ 42c40192-3ac9-4259-be11-8fc508f9611e
r_66CL = r.(1 ./ invg_66CL)

# ╔═╡ 562878be-01c2-4aa9-9547-7b63dd9edd64
md"""
#### Result: $(round(u"fm", r_66CL[2], digits=2)) < r < $(round(u"fm", r_66CL[1], digits=2)) with 66% CL
"""

# ╔═╡ a35aabc0-1e04-4b8e-9883-610209029d1f
md"""
## Compositeness
$1-Z = \sqrt{\frac{1}{1+2r/a}}$
"""

# ╔═╡ 1186fdbf-46ef-42bf-9ecf-90ae649ed899
Z(g) = 1-sqrt(1/(1+2r(g)*inva0/fm_times_mev))

# ╔═╡ e6a09750-cce8-4728-8b8f-4fccfa2d258d
Z0 = Z(gLHCb)

# ╔═╡ fcd8338f-b68e-401f-8832-8800110db9c1
Z_66CL = Z.(1 ./ invg_66CL)

# ╔═╡ f026f9cd-8b91-480e-9b09-bb8beb1cea84
md"""
#### Result: $(round(Z_66CL[1], digits=2)) < Z < $(round(Z_66CL[2], digits=2)) with 66% CL
"""

# ╔═╡ Cell order:
# ╠═61b01d0b-1f8f-4e9d-805a-b533cd39410c
# ╟─74db3011-7c88-40dc-ae78-bd9c3e6be856
# ╠═8a8ef3e6-1277-407c-b906-7fb173230a7d
# ╠═033a5372-bdd0-4065-8147-7a9a375ed3d4
# ╠═a32cb39c-5ad0-457e-8e78-6dc1afe795b8
# ╠═1c80959c-ed0a-44e9-8da1-b16e1d0e3bc0
# ╟─1eb7633e-4c75-4d18-96d5-bf2fadf8e5cd
# ╠═cc0af199-3b96-4f8b-a31a-706bdf716365
# ╠═5f9294d6-68a6-4cde-b59d-4a9790d63668
# ╟─e73dbf75-53d8-4ec4-b23e-3a15abb3a28b
# ╟─7f7af133-c60c-4317-b357-9840eebb4785
# ╠═0e275bfc-9d93-4d02-957b-de9fc77c1964
# ╠═782412e1-b150-459a-b973-3e7e27c20408
# ╟─19740af1-19cc-482f-83fe-613ada9a55f8
# ╠═43bcd184-ffe7-40da-832c-025028649272
# ╠═89d036a1-1c15-40ee-890a-cecfe8857297
# ╠═42c40192-3ac9-4259-be11-8fc508f9611e
# ╟─562878be-01c2-4aa9-9547-7b63dd9edd64
# ╟─a35aabc0-1e04-4b8e-9883-610209029d1f
# ╠═1186fdbf-46ef-42bf-9ecf-90ae649ed899
# ╠═e6a09750-cce8-4728-8b8f-4fccfa2d258d
# ╠═fcd8338f-b68e-401f-8832-8800110db9c1
# ╟─f026f9cd-8b91-480e-9b09-bb8beb1cea84
