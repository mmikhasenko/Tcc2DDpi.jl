### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 2fae0580-b7b9-11ed-1301-cd79a7c1ff1b
# ╠═╡ show_logs = false
begin
	cd(joinpath(@__DIR__, ".."))
	import Pkg
	Pkg.activate(".")
	Pkg.instantiate()
	# 
	using SymPy
	#
	using Plots
	using Parameters
	using LaTeXStrings
	using JSON
end

# ╔═╡ 0e34c3ac-2601-495f-a538-cae3fe21622e
md"""
# Triangle singularity in $\pi^+D^0D^0$

Accroding to [Schmid theorem](https://inspirehep.net/literature/1474944) the triangle singulary is not physics for the elastic scattering. However, it is still located in the complex plane on the unphysical sheet, hidden by the $D^{*+}D^0$ cut.

This notebook finds the exact location of the logarithmic branch points using the [Landau rules](https://inspirehep.net/literature/4561) and plots its location on the complex plane energy plane.
"""

# ╔═╡ 7a1b1e69-e221-478f-a6b0-4edd1cff7c0a
theme(:vibrant, size=(500, 350), minorticks=true, grid=false, frame=:box,
    guidefontvalign=:top, guidefonthalign=:right,
    foreground_color_legend=nothing,
    legendfontsize=9, legend=:topright, lab="")

# ╔═╡ 89b0dc9f-6f32-48b5-8b12-bb53e6d86197
begin
	const mD⁰ = 1.86483
	const mDˣ⁺ = 1.86483+145.4258e-3 # m(D) + Δm(D*,D) from PDG
	const ΓDˣ⁺ = 83.4e-6
	const mπ⁺ = 0.13957039
	const mth = mDˣ⁺ + mD⁰
	const mcth = sqrt(mDˣ⁺^2 - 1im*mDˣ⁺*ΓDˣ⁺) + mD⁰
	# 
	m2e(m) = (m-mth)*1e3
end;

# ╔═╡ 3e86486b-70d6-4318-bb29-1f25683bae99
begin
	@syms k1 k2 k3
	@syms s t
	@syms μ1 μ2 μ3
	@syms m0 m1 m2 m3
	@syms x # m0²
	@syms z # mR²
end;

# ╔═╡ aded2e35-52d6-4dd9-a48a-77ee507fb4ae
begin
	plot(xaxis=false, yaxis=false)
	plot!([-1,1+2im,1-2im,-1, -2], c=2, l=(3))
	plot!([2+2im,1+2im,1-2im,2-2im], l=(3))
	scatter!([-1,1+2im,1-2im,-1], mc=2, ms=10)
	annotate!([
		(-0.3, 1.2,"Dˣ⁺: k₁",:right),
		(-0.3,-1.2,"D⁰: k₃",:right),
		(1.2,0.0,"π⁺: k₂",:left)])
	plot!(xlab="", ylab="")
end

# ╔═╡ 112ee1a2-597d-4289-aab4-c923d188d638
M = [k1 k2 k3] .* [k1, k2, k3]

# ╔═╡ aa9e9fc0-83fe-4eea-a870-cdac1ee99f2c
repl = Dict(
	k1^2=>μ1^2,
	k2^2=>μ2^2,
	k3^2=>μ3^2,
	k1*k3=>(μ1^2+μ3^2-m0^2)/2,
	k1*k2=>(μ1^2+μ2^2-m1^2)/2,
	k2*k3=>(μ2^2+μ3^2-s)/2);

# ╔═╡ b6aa7110-ef01-4cc9-bb9b-c7a79811da97
vars = Dict(μ1^2=>t, μ2^2=>m2^2, μ3^2=>m3^2, m0^2=>x);

# ╔═╡ 3e55b255-66a0-482f-bf2d-72be39817b9f
md"""
Here is the Landau condition: the determinant of the matrix is equal to zero.
"""

# ╔═╡ 62138207-78b6-495e-aea4-b4e14dae46b9
Mrepl = map(M) do Mi
	Mi.xreplace(repl).xreplace(vars)
end

# ╔═╡ e898eecd-27ee-4ba3-932c-2a274c4af085
detM = sympy.det(Mrepl);

# ╔═╡ 62f6f333-3f08-4a17-8e34-07bfca80ee8a
x12 = sympy.solve(detM, x)

# ╔═╡ ba5fa039-fa1e-4139-9e0f-0eefb58fca8e
md"""
**A check**: the [a₁(1420) triangle singularity](https://inspirehep.net/literature/1341619) is located around 1.42GeV,
when the masses correspond to the KˣK->f₀π scattering
"""

# ╔═╡ c4f67c75-1ff1-4b1c-a0b8-5ac5f6d5bcf3
repl_a11420 = Dict(t=>0.892^2, m1=>0.15, m2=>0.49, m3=>0.49, s=>0.981^2)

# ╔═╡ 2032f939-5100-4e84-9753-f58fd928513b
map(x12) do x
	N(sympy.sqrt(x.xreplace(repl_a11420)))
end

# ╔═╡ 2bf70eee-c06e-44cb-9f41-e07bc288918a
md"""
Now, the masses for $D^{*+}D^{0}\to D^{*+}D^{0}$ are inserted with a complex value for the mass of the $D^{*+}$ meson.
"""

# ╔═╡ 9c040074-49b3-4419-bfe5-51cb390f49b3
repl_Tcc = Dict(t=>z, m1=>mD⁰, m2=>mπ⁺+0.005, m3=>mD⁰, s=>z); 

# ╔═╡ a2f1e1d6-1407-4ef4-9536-e6349f21bec0
repl_Dx = Dict(z=>mDˣ⁺^2 - 1im*mDˣ⁺*ΓDˣ⁺,);

# ╔═╡ 7e9f91b8-9fe1-4901-a9ad-8b5d36e727e7
TSpoints = map(x12) do x
	N(sympy.sqrt(x.xreplace(repl_Tcc).xreplace(repl_Dx)))
end .|> m2e

# ╔═╡ 391f21fe-dcc7-4f39-a104-2247a32df51c
let ymin=-0.1
	plot(title="Analytic structure", xlim=(-0.1,0.3), ylim=(ymin,0.01), frame=:origin)
	plot!([m2e(mcth), real(m2e(mcth))+1im*ymin], l=(3, :red), lab="cut")
	scatter!([m2e(mcth)], m=(6, :red), lab="threshold")
	plot!(TSpoints, l=(3,:dash), lc=3, lab="log cut")
	scatter!(TSpoints, ms=8, m=:d, mc=3, lab="triangle singularity")
	plot!(xlab="Re E", ylab="Im E", leg=:bottomright)
	#
	savefig(joinpath("plots", "trainglesingularity.pdf"))
	plot!()
end

# ╔═╡ Cell order:
# ╟─0e34c3ac-2601-495f-a538-cae3fe21622e
# ╠═2fae0580-b7b9-11ed-1301-cd79a7c1ff1b
# ╠═7a1b1e69-e221-478f-a6b0-4edd1cff7c0a
# ╠═89b0dc9f-6f32-48b5-8b12-bb53e6d86197
# ╠═3e86486b-70d6-4318-bb29-1f25683bae99
# ╠═aded2e35-52d6-4dd9-a48a-77ee507fb4ae
# ╠═112ee1a2-597d-4289-aab4-c923d188d638
# ╠═aa9e9fc0-83fe-4eea-a870-cdac1ee99f2c
# ╠═b6aa7110-ef01-4cc9-bb9b-c7a79811da97
# ╟─3e55b255-66a0-482f-bf2d-72be39817b9f
# ╠═62138207-78b6-495e-aea4-b4e14dae46b9
# ╠═e898eecd-27ee-4ba3-932c-2a274c4af085
# ╠═62f6f333-3f08-4a17-8e34-07bfca80ee8a
# ╟─ba5fa039-fa1e-4139-9e0f-0eefb58fca8e
# ╠═c4f67c75-1ff1-4b1c-a0b8-5ac5f6d5bcf3
# ╠═2032f939-5100-4e84-9753-f58fd928513b
# ╟─2bf70eee-c06e-44cb-9f41-e07bc288918a
# ╠═9c040074-49b3-4419-bfe5-51cb390f49b3
# ╠═a2f1e1d6-1407-4ef4-9536-e6349f21bec0
# ╠═7e9f91b8-9fe1-4901-a9ad-8b5d36e727e7
# ╠═391f21fe-dcc7-4f39-a104-2247a32df51c
