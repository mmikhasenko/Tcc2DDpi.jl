### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ ac869de2-bdbf-11ed-112d-15f08551dee2
# ╠═╡ show_logs = false
begin
	using Pkg
	Pkg.add([
		Pkg.PackageSpec("PyCall"),
		Pkg.PackageSpec("SymPy"),
		Pkg.PackageSpec("Parameters"),
		Pkg.PackageSpec("Plots"),
		Pkg.PackageSpec("DataFrames"),
		])
	# 
	using Parameters
	using SymPy
	using Plots
	import Plots.PlotMeasures: mm
	using DataFrames
	# 
	import PyCall
	PyCall.pyimport_conda("sympy.physics.wigner",       "sympy")
	PyCall.pyimport_conda("sympy.physics.quantum.spin", "sympy")
	# 
	import_from(sympy.physics.wigner)
	# import_from(sympy.physics.quantum.cg, )
    import_from(sympy.physics.quantum.spin, (:WignerD, :CG), typ=:Any)
end

# ╔═╡ 0ae0767b-6b8c-490f-8d5c-8dd89b662ade
md"""
# Effective range with stable $D^*$

Expression for the effective range is analytic for two-particle scattering.
In this notebook, we compute it using sympy.
"""

# ╔═╡ 1936cdb9-70e2-4716-9778-2a4eee9a0103
theme(:vibrant, size=(500, 350), minorticks=true, grid=false, frame=:box,
    guidefontvalign=:top, guidefonthalign=:right,
    foreground_color_legend=nothing,
    legendfontsize=9, legend=:topright, lab="",
    xlim=(:auto, :auto), ylim=(:auto, :auto))

# ╔═╡ 64d893db-8971-45f6-b01e-bf5d61879bba
Kallen(x,y,z) = x^2+y^2+z^2-2x*y-2y*z-2z*x

# ╔═╡ a4031742-03ef-4c96-ac56-60a344af9d17
q(s,M1,m1) = sqrt(Kallen(s,M1^2,m1^2))/(2sqrt(s))

# ╔═╡ 955d81ad-6d32-441b-ac73-9860be3801a5
Σ(s,m1,m2) =
	sqrt(Kallen(s,m1^2,m2^2)) / s *
	log((m1^2+m2^2-s+sqrt(Kallen(s,m1^2,m2^2)))/(2m1*m2)) -
	(m1^2-m2^2)/s*log(m1/m2)

# ╔═╡ ab4ad8bb-5c14-4165-b73c-0f7027acc850
imΣ(s,m1,m2) = sqrt(Kallen(s,m1^2,m2^2)) / s * π

# ╔═╡ 92f9527a-b8c5-453f-bd65-adafb07b7400
reΣ(s,m1,m2) = Σ(s,m1,m2) - 1im*imΣ(s,m1,m2) 

# ╔═╡ e1305876-a53e-4221-a6b2-8496ca8c41f4
let
	plot()
	plot!(s->real(Σ(s+1e-7im, 1, 1)), -1, 7, lab="real")
	plot!(s->imag(Σ(s+1e-7im, 1, 1)), -1, 7, lab="imag")
end

# ╔═╡ 5ee4de99-290d-4b7d-8276-d59fe548560c
md"""
The reaction amplitude is written as $1/D(s)$, where the denominator
incorporates the dynamics of two coupled channels $(M_1,m_1)$, and $(M_2,m_2)$.

The function reads as a two two-body loop function, the subtraction point is set by the binding energy. A pole below the lowest threshold is set by making the real part vanish at $s=m_0^2$.
"""

# ╔═╡ 075ae401-b8ac-4f3f-bc35-11c1f9fe2c01
D(s,M1,m1,M2,m2,m0) =
	real(Σ(m0^2,M1,m1)+Σ(m0^2,M2,m2)) -
	(Σ(s,M1,m1)+Σ(s,M2,m2))

# ╔═╡ 673a2696-8fa4-4fdb-b93d-2b72b0afc593
begin
	R1(s,M1,m1,M2,m2,m0) = D(s,M1,m1,M2,m2,m0) + 1im*(imΣ(s,M1,m1))
	R2(s,M1,m1,M2,m2,m0) = D(s,M1,m1,M2,m2,m0) + 1im*(imΣ(s,M2,m2))
end ;

# ╔═╡ f6182faf-5a3d-4b69-bd2c-590aa2050394
md"""
## Effective range expansion
"""

# ╔═╡ f25eec38-79c5-444b-accc-855b3aa07d1a
begin
	ere(k; a⁻¹, r, N) = N*(a⁻¹ + r * k^2 / 2 - 1im* k)
	ere(k, p) = ere(k, a⁻¹ = p.a⁻¹, r = p.r, N = p.N)
	ereR(k, p) = ere(k, p)-ere(k, (;p..., r=0, a⁻¹=0))
end ;

# ╔═╡ eceebe12-2394-4a57-8d03-d3a91702e8aa
md"""
### Symbolic computation

Symbolic expression for the effective range at the first threshold

The expansion for the regular part of the denominator reads:

$R = N \left( \frac{1}{a} + r \frac{k^2}{2}+O(a^3k^4) \right)$

Hence, the expansion coeffients are found as 

```math
\begin{align}
a^{-1} &= R(e_\text{th})/N\,,\\
r &=  \left.\frac{dR}{d(k^2)}\right|_{e_\text{th}} =
	\left.\frac{d R}{d e}\right|_{e_\text{th}}
	\bigg/
	\left.\frac{d k^2}{d e}\right|_{e_\text{th}}\,.
\end{align}
```
Here is a symbolic expression for the $r$:
"""

# ╔═╡ 74abf809-f627-4ea5-9f31-0b94c82faab2
@syms e ;

# ╔═╡ 3b55da67-e23c-4a8f-83c3-b5d702b8e205
masssymbols = let
	@syms M1::positive m1::positive M2::positive m2::positive m0::real
	(; M1,m1,M2,m2,m0)
end ;

# ╔═╡ 1eea9231-fcd5-4d82-86e0-d188de8810f1
r1_ana = let
	@unpack M1, m1, M2, m2, m0 = masssymbols
	# 
	N1 = 2PI/(M1+m1)
	# 
	dk1² = sympy.diff(q((M1+m1+e)^2,M1,m1)^2, e) |> simplify
	dR1 = sympy.diff(R1((M1+m1+e)^2,M1,m1,M2,m2,m0), e)
	# 
	lambdify(2 / N1 * dR1 / dk1², (e, M1,m1,M2,m2))
end

# ╔═╡ 46742146-1341-4083-9c4c-31e7c755b839
r2_ana = let
	@unpack M1, m1, M2, m2, m0 = masssymbols
	# 
	N2 = 2PI/(M2+m2)
	# 
	dk2² = sympy.diff(q((M1+m1+e)^2,M2,m2)^2, e) |> simplify
	dR2 = sympy.diff(R2((M1+m1+e)^2,M1,m1,M2,m2,m0), e)
	# 
	lambdify(2 / N2 * dR2 / dk2², (e, M1,m1,M2,m2))
end

# ╔═╡ fc3a7973-ed0e-4a8c-aa36-504ceb0917fa
md"""
## Example of $T_{cc}^+$
"""

# ╔═╡ 5466eb08-17c8-47ec-ac3e-94aea5d6dca5
begin
	const fm_times_mev = 197.3269804
	# masses
	const mπ⁰ = 0.1349768
	const mπ⁺ = 0.13957039
	const mD⁰ = 1.86483
	const mD⁺ = 1.86965
	# 
	const mDˣ⁺ = 1.86483+145.4258e-3 # m(D) + Δm(D*,D) from PDG
	const mDˣ⁰ = 2.00685
	# 
	const ΓDˣ⁺ = 83.4e-6
	const ΓDˣ⁰ = 55.2e-6
	# 
	const δm = -0.3592 # MeV
	# 
	const iϵ = 1e-8im
	# 
	e2m(e) = (mDˣ⁺+mD⁰)+e*1e-3
	# 
	eth1 = 0.0+iϵ
	eth2 = ((mDˣ⁰+mD⁺) - (mDˣ⁺+mD⁰))*1e3+iϵ
	# 
	
	massvalues_Tcc = let
		@unpack M1, m1, M2, m2, m0 = masssymbols
		Dict(M1=>mDˣ⁺, m1=>mD⁰, M2=>mDˣ⁰, m2=>mD⁺)
	end
end ;

# ╔═╡ 41c17061-0543-4d11-9199-a069dd6a82e6
D_Tcc(e) = D(e2m(e)^2,
	mDˣ⁺,mD⁰,
	mDˣ⁰,mD⁺,
	e2m(δm)+iϵ) # m0

# ╔═╡ 3a0736a0-635a-4a88-aaf4-f1f62d3c43c5
begin
	R1_Tcc(e) = R1(e2m(e)^2, mDˣ⁺,mD⁰, mDˣ⁰,mD⁺, e2m(δm)+iϵ)
	R2_Tcc(e) = R2(e2m(e)^2, mDˣ⁺,mD⁰, mDˣ⁰,mD⁺, e2m(δm)+iϵ) 
end

# ╔═╡ 38b7f5a1-5021-4387-9deb-fea938b774b8
k1(e) = q(e2m(e)^2,mDˣ⁺,mD⁰)

# ╔═╡ 187aebd1-2fe0-47a5-b960-460e5b3bfeb9
k2(e) = q(e2m(e)^2,mDˣ⁰,mD⁺)

# ╔═╡ b9b3941b-28e8-4245-9ade-5e2e8e31966d
let
	plot(title = "D and R1, R2")
	plot!(e->real(D_Tcc(e+iϵ)),-1, 3, lab="Re D")
	plot!(e->imag(D_Tcc(e+iϵ)),-1, 3, lab="Im D")
	# 
	plot!(e->real(R1_Tcc(e+iϵ)),-1, 3, lab="R1", c=1, ls=:dash)
	plot!(e->imag(R1_Tcc(e+iϵ)),-1, 3, lab="", c=2, ls=:dash)
	# 
	plot!(e->real(R2_Tcc(e+iϵ)),-1, 3, lab="R2", c=1, ls=:dashdot)
	plot!(e->imag(R2_Tcc(e+iϵ)),-1, 3, lab="", c=2, ls=:dashdot)
end

# ╔═╡ 0159c179-a8fb-4342-b60e-7bbd4875d188
md"""
### Numerical computation
with analytic diffs
"""

# ╔═╡ 76b90397-ae34-47db-b317-d9dccb482f5e
N(e) = 2π/e2m(e)

# ╔═╡ b91ddae1-f221-413f-a745-e259c880d251
inva(e) = D_Tcc(e)/N(e)

# ╔═╡ ecca02bb-8a67-4a52-a355-88270d1f260f
ere1 = let # derivative is computed analytically
	N1 = N(eth1)
	inva1 = inva(eth1)
	# 
	dR1 = sympy.diff(R1_Tcc(e), e);
	λdR1 = dR1.replace(e,eth1) |> sympy.N
	# 
	dk1² = sympy.diff(k1(e)^2, e) ;
	λdk1² = dk1².replace(e,eth1) |> sympy.N
	# 
	sr1 = 2/N1 * λdR1 / λdk1²
	r1 = convert(Complex{Float64}, sr1)
	# 
	(; N=N1, a⁻¹=inva1, r=r1)
end

# ╔═╡ cc09971b-1411-47e7-a22d-eaefb1b890a5
md"""
The analytic formula gives an incorrect value because of evaluation of the sqrt.
"""

# ╔═╡ ab4e52bd-7164-422b-9b8d-f84c3dfab5a9
ere2 = let # derivative is computed analytically
	N2 = N(eth2)
	inva2 = inva(eth2)
	# 
	dR2 = sympy.diff(R2_Tcc(e), e);
	λdR2 = dR2.replace(e,eth2) |> sympy.N
	# 
	dk2² = sympy.diff(k2(e)^2, e) ;
	λdk2² = dk2².replace(e,eth2) |> sympy.N
	# 
	sr2 = 2/N2 * λdR2 / λdk2² |> sympy.N
	r2 = convert(Complex{Float64}, sr2)
	# 
	(; N=N2, a⁻¹=inva2, r=r2)
end

# ╔═╡ 85f75140-b103-41e4-86ad-a9680f95fb8c
let
	plot(title = "ERE against the model")
	# 
	plot!(e->real(D_Tcc(e+iϵ)),-1, 2, lab="real")
	plot!(e->imag(D_Tcc(e+iϵ)),-1, 2, lab="imag")
	# 
	plot!(e->real(ere(k1(e+iϵ), ere1)),-1, 2, lab="ERE(1)", c=1, ls=:dash)
	plot!(e->imag(ere(k1(e+iϵ), ere1)),-1, 2, lab="", c=2, ls=:dash)
	# 
	plot!(e->real(ere(k2(e+iϵ), ere2)),-1, 2, lab="ERE(2)", c=1, ls=:dashdot)
	plot!(e->imag(ere(k2(e+iϵ), ere2)),-1, 2, lab="", c=2, ls=:dashdot)
end

# ╔═╡ e20eaee2-ecc1-487e-a5f3-1226bb13f0a1
let
	plot(layout=grid(1,2), size=(800,350),
		title=["R₁" "R₂"], xlab="Δm (MeV)", margin=3mm)
	# 
	plot!(sp=1, e->real(R1_Tcc(e+iϵ)),-1, 3, lab="real")
	plot!(sp=1, e->imag(R1_Tcc(e+iϵ)),-1, 3, lab="imag")
	# 
	plot!(sp=1, e->real(ereR(k1(e+iϵ), ere1)),-1, 3, lab="ERE(1)", c=1, ls=:dash)
	plot!(sp=1, e->imag(ereR(k1(e+iϵ), ere1)),-1, 3, lab="", c=2, ls=:dash)
	# 
	vline!(sp=1, [real(eth1)], lab="", lc=:black)
	# 
	plot!(sp=2, e->real(R2_Tcc(e+iϵ)),-1, 3, lab="real")
	plot!(sp=2, e->imag(R2_Tcc(e+iϵ)),-1, 3, lab="imag")
	# 
	plot!(sp=2, e->real(ereR(k2(e+iϵ), ere2)),-1, 3, lab="ERE(2)", c=1, ls=:dash)
	plot!(sp=2, e->imag(ereR(k2(e+iϵ), ere2)),-1, 3, lab="", c=2, ls=:dash)
	# 
	vline!(sp=2, [real(eth2)], lab="", lc=:black)
end

# ╔═╡ f284efab-8653-452e-9a28-efa83448dbe5
function tophysicsunits(p::NamedTuple)
    @unpack a⁻¹, r = p
    a_fm = 1e-3*fm_times_mev / real(p.a⁻¹)
    r_fm = 1e-3*fm_times_mev * p.r
    (; a_fm, r_fm)
end

# ╔═╡ 187b1315-ff3a-408f-8369-2e0b61f838f8
begin
	df = DataFrame()
	push!(df, (;expantion_point = :(Dˣ⁺+D⁰), ere1..., tophysicsunits(ere1)...))
	push!(df, (;expantion_point = :(Dˣ⁰+D⁺), ere2..., tophysicsunits(ere2)...))
	# 
	# add symbolic
	df.r_ana = 
		[r1_ana(eth1*1e-3, mDˣ⁺,mD⁰, mDˣ⁰,mD⁺),
	 	 r2_ana(eth2*1e-3, mDˣ⁺,mD⁰, mDˣ⁰,mD⁺)]
	#
	select(df, :expantion_point,
		[:a⁻¹, :a_fm, :r_fm] .=> ByRow(x->round(x; digits=3)) .=>
		[:a⁻¹_GeV, "1/Re[a⁻¹]_fm", :r_fm],
		:r_ana => ByRow(x->round(1e-3*fm_times_mev * x; digits=3)) => :r_ana_fm)
end

# ╔═╡ Cell order:
# ╟─0ae0767b-6b8c-490f-8d5c-8dd89b662ade
# ╠═ac869de2-bdbf-11ed-112d-15f08551dee2
# ╠═1936cdb9-70e2-4716-9778-2a4eee9a0103
# ╠═64d893db-8971-45f6-b01e-bf5d61879bba
# ╠═a4031742-03ef-4c96-ac56-60a344af9d17
# ╠═955d81ad-6d32-441b-ac73-9860be3801a5
# ╠═ab4ad8bb-5c14-4165-b73c-0f7027acc850
# ╠═92f9527a-b8c5-453f-bd65-adafb07b7400
# ╠═e1305876-a53e-4221-a6b2-8496ca8c41f4
# ╟─5ee4de99-290d-4b7d-8276-d59fe548560c
# ╠═075ae401-b8ac-4f3f-bc35-11c1f9fe2c01
# ╠═673a2696-8fa4-4fdb-b93d-2b72b0afc593
# ╟─f6182faf-5a3d-4b69-bd2c-590aa2050394
# ╠═f25eec38-79c5-444b-accc-855b3aa07d1a
# ╟─eceebe12-2394-4a57-8d03-d3a91702e8aa
# ╠═74abf809-f627-4ea5-9f31-0b94c82faab2
# ╠═3b55da67-e23c-4a8f-83c3-b5d702b8e205
# ╠═1eea9231-fcd5-4d82-86e0-d188de8810f1
# ╠═46742146-1341-4083-9c4c-31e7c755b839
# ╟─fc3a7973-ed0e-4a8c-aa36-504ceb0917fa
# ╠═5466eb08-17c8-47ec-ac3e-94aea5d6dca5
# ╠═41c17061-0543-4d11-9199-a069dd6a82e6
# ╠═3a0736a0-635a-4a88-aaf4-f1f62d3c43c5
# ╠═38b7f5a1-5021-4387-9deb-fea938b774b8
# ╠═187aebd1-2fe0-47a5-b960-460e5b3bfeb9
# ╠═b9b3941b-28e8-4245-9ade-5e2e8e31966d
# ╟─0159c179-a8fb-4342-b60e-7bbd4875d188
# ╠═76b90397-ae34-47db-b317-d9dccb482f5e
# ╠═b91ddae1-f221-413f-a745-e259c880d251
# ╠═ecca02bb-8a67-4a52-a355-88270d1f260f
# ╟─cc09971b-1411-47e7-a22d-eaefb1b890a5
# ╠═ab4e52bd-7164-422b-9b8d-f84c3dfab5a9
# ╠═85f75140-b103-41e4-86ad-a9680f95fb8c
# ╠═e20eaee2-ecc1-487e-a5f3-1226bb13f0a1
# ╠═187b1315-ff3a-408f-8369-2e0b61f838f8
# ╠═f284efab-8653-452e-9a28-efa83448dbe5
