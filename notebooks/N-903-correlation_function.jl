### A Pluto.jl notebook ###
# v0.19.12

using Markdown
using InteractiveUtils

# ╔═╡ 61b01d0b-1f8f-4e9d-805a-b533cd39410c
begin
    import Pkg
    Pkg.activate(@__DIR__)
end

# ╔═╡ e0bd7550-0ebb-11ee-2f33-451a4d06e98a
begin
	using Plots
	using LaTeXStrings
end

# ╔═╡ 7ea80f9f-abd5-4217-a716-9fe0e1b54c85
md"""
# Mass spectrum of correlated $D^{*+}D^0$

The notebook computes how the mass spectrum of the $D^{*+}D^0$ would look like give the correlation function computed using the femtoscopy technique.
"""

# ╔═╡ be68d089-6ca5-41a0-ba2f-cd28e447fa08
theme(:wong2, xlim=(:auto,:auto), frame=:box)

# ╔═╡ 0e57bf54-f0a9-48ae-b420-f0178ab50257
md"""
## Correlation function
I make up the function look like the computatation in Fig.3 of [2303.06079](https://arxiv.org/pdf/2303.06079.pdf).
"""

# ╔═╡ 26a520f5-36da-4ce7-97aa-a4dc9d5c4321
begin
	g(x) = (15/(1+x) - 10/(1+x^2))
	f(x) = g(x/80-0.75) /  g(300/80-0.75)
end

# ╔═╡ 7d90ead9-616e-49f9-818f-ac307ccc6239
begin
	plot(f, 0, 300, xlab="kˣ (MeV)", ylab="C(kˣ)", lab="Dˣ⁺D⁰")
	hline!([f(300)], lab="=1")
end

# ╔═╡ 0a0a0711-de80-425b-aee9-e8a3de0034ea
begin
	const mD⁰ = 1.86483
	const mDˣ⁺ = 1.86483 + 145.4258e-3 # m(D) + Δm(D*,D) from PDG
end

# ╔═╡ 967ac230-2e05-480f-a126-b539cf6e3ba7
begin
	Kallen(x,y,z) = x^2+y^2+z^2-2x*y-2y*z-2z*x
	k(m) = sqrt(Kallen(m^2,mD⁰^2,mDˣ⁺^2))/2m * 1e3
	e2m(e) = e/1e3+mD⁰+mDˣ⁺
end

# ╔═╡ be619072-a714-4107-94e1-3acc78e800d0
begin
	plot(xlab="m(Dˣ⁺D⁰)-m(Dˣ⁺)-m(D⁰) (MeV)")
	plot!(f ∘ k ∘ e2m, 0, 12, ylab="C", lab="Dˣ⁺D⁰")
	hline!([f(300)], lab="=1")
end

# ╔═╡ f0d66366-a5c4-4a0a-b095-be74fe58eea6
begin
	plot(xlab="m(Dˣ⁺D⁰)-m(Dˣ⁺)-m(D⁰) (MeV)", ylab="#N(Dˣ⁺D⁰)")
	plot!(e->f(k(e2m(e)))*k(e2m(e)), 0, 12, lab="Dˣ⁺D⁰", l=(2,:red))
	plot!(e->k(e2m(e)), 0, 12, lab="ph.sp.", l=(1,:gray))
end

# ╔═╡ Cell order:
# ╠═61b01d0b-1f8f-4e9d-805a-b533cd39410c
# ╟─7ea80f9f-abd5-4217-a716-9fe0e1b54c85
# ╠═e0bd7550-0ebb-11ee-2f33-451a4d06e98a
# ╠═be68d089-6ca5-41a0-ba2f-cd28e447fa08
# ╟─0e57bf54-f0a9-48ae-b420-f0178ab50257
# ╠═26a520f5-36da-4ce7-97aa-a4dc9d5c4321
# ╠═7d90ead9-616e-49f9-818f-ac307ccc6239
# ╠═0a0a0711-de80-425b-aee9-e8a3de0034ea
# ╠═967ac230-2e05-480f-a126-b539cf6e3ba7
# ╠═be619072-a714-4107-94e1-3acc78e800d0
# ╠═f0d66366-a5c4-4a0a-b095-be74fe58eea6
