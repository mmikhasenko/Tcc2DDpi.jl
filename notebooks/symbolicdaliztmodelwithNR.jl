### A Pluto.jl notebook ###
# v0.19.12

using Markdown
using InteractiveUtils

# ╔═╡ 3dbb29a8-1f37-407f-8907-698c80f2bd00
# ╠═╡ show_logs = false
begin
	cd(joinpath(@__DIR__, ".."))
	import Pkg
	Pkg.activate(".")
	Pkg.instantiate()
	# 
	using SymPy
	import PyCall
	PyCall.pyimport_conda("sympy.physics.quantum.spin", "sympy")
	import_from(sympy.physics.quantum.spin, (:CG, ), typ=:Any)
	import_from(sympy.physics.quantum.spin.Rotation, (:d,), typ=:Any)
	#
	using Parameters
	using LaTeXStrings
	using JSON
end

# ╔═╡ 04808f8e-25bf-4c6c-9794-5a9d59b197e8
begin
	wignerd(J,λ1,λ2,θ) = d(Sym(J),λ1,λ2,θ)
	clgd(j1,λ1,j2,λ2,j,λ) = CG(j1,λ1,j2,λ2,j,λ)
end

# ╔═╡ f907394e-a7e0-43d0-8394-328c01dcca7c
begin
	@syms(
		θ12::nonnegative=>"\\theta_{12}",
		θ31::nonnegative=>"\\theta_{31}")
	@syms(
		ζ23_for0::nonnegative=>"\\zeta^0_{2(3)}",
		ζ13_for0::nonnegative=>"\\zeta^0_{1(3)}",
		ζ23_for1::nonnegative=>"\\zeta^1_{2(3)}")
end ;

# ╔═╡ f9f95dab-babb-4190-82af-f1681ce264b4
md"""
### Amplitude $\pi DD$
"""

# ╔═╡ 71f9928c-42c2-43c2-899e-6fa756f6b683
abstract type πDD end

# ╔═╡ 078cad0e-7202-4ec5-b829-1b44f029160e
struct πDD_2xDˣ_3xB{T} <: πDD
	F2::T
	F3::T
	B1::T # DD S-wave
	B2::T # Dπ S-wave
	B3::T # Dπ S-wave
end

# ╔═╡ ca2e840e-df9d-4a27-834d-f6e2f1737378
function amplitude(ch::πDD_2xDˣ_3xB, j0, L, λ)
	@unpack F2, F3 = ch
	@unpack B1, B2, B3 = ch
	A = Sym(0)
	# 
	j_F = 1
	S_F = 1
	# 
	ζ33_for0 = 0 # trivial
	# chain-2
	A += sqrt(Sym(2L+1)) * F2 * sum(
		wignerd(j0,λ,τ,ζ23_for0)*
			wignerd(j_F,τ,0,θ31)*clgd(L,0,S_F,τ,j0,τ)
		for τ in -j_F:j_F)
	# 
	# chain-3
	A += (-1)^j_F * # H_{Dπ} vs H_{πD} in helicity basis
		sqrt(Sym(2L+1)) * F3 * sum(
			wignerd(j0,λ,τ,ζ33_for0)*wignerd(j_F,τ,0,θ12)*clgd(L,0,S_F,τ,j0,τ)
		for τ in -j_F:j_F)
	#
	# j_NR = 0 # of Dpi
	# S_NR = 0 # scalar dimer + pseudoscalar
	# L_NR = j0 # only correct for unnatural parity of X 0-, 1+, 2-
	τ_Dpi = 0
	τ_DD = 0
	# backgr-1
	A += B1 *
		wignerd(j0,λ,τ_DD,ζ13_for0)
			# wignerd(j_NR,τ_DD,0,θ23)
			# clgd(L_NR,0,S_NR,τ_DD,j0,τ_DD)
	# backgr-2
	A += B2 *
		wignerd(j0,λ,τ_Dpi,ζ23_for0)
			# wignerd(j_NR,τ_Dpi,0,θ31)
			# clgd(L_NR,0,S_NR,τ_Dpi,j0,τ_Dpi)
	# backgr-3
	A += B3 *
		wignerd(j0,λ,τ_Dpi,ζ33_for0)
			# wignerd(j_NR,τ_Dpi,0,θ12)
			# clgd(L_NR,0,S_NR,τ_Dpi,j0,τ_Dpi)
	return A
end

# ╔═╡ 7065b936-6ede-4cd7-9290-8e453d903f48
intensity(ch::πDD, j0, L) = sum(abs2, amplitude(ch, j0, L, λ) for λ in -j0:j0)

# ╔═╡ 7608afab-8d37-4578-b1ad-7641de820d79
@syms(
	F2_Dπ::real=>"\\mathcal{F}_2^{D\\pi}",
	F3_Dπ::real=>"\\mathcal{F}_3^{D\\pi}",
	# 
	B2_Dπ::real=>"\\mathcal{B}_2^{D\\pi}",
	B3_Dπ::real=>"\\mathcal{B}_3^{D\\pi}",
	# 
	B1_DD::real=>"\\mathcal{B}_1^{DD}"
)

# ╔═╡ 444842d4-8b3d-4921-9822-4beb2923fa3f
const chains = (F2_Dπ,F3_Dπ,B1_DD,B2_Dπ,B3_Dπ)

# ╔═╡ 5e7042ca-de8e-4347-93bc-9125fe9562b0
_πDD = πDD_2xDˣ_3xB(chains...)

# ╔═╡ 905dbfa9-9d6e-4d2b-b048-a72cb75f687e
I10 = intensity(_πDD, 1, 0) ;

# ╔═╡ 48a50729-332d-4c9d-9b93-b4912f32581d
I10_flat = I10.doit() |> expand |> sympy.trigsimp

# ╔═╡ 1509915d-4e8f-47d5-8db7-091ab21d799c
I10_poly = sympy.poly(I10_flat, chains...) ;

# ╔═╡ 03f5100d-2b80-437e-9257-be90a44808d8
I10_dict = I10_poly.as_dict()

# ╔═╡ 8c23e629-5a24-4c4e-8016-84916badfccf
function decompose(k)
	sum(k) != 2 && error("sum(k)!=2")
	if sum(k .== 2) != 0
		i = findfirst(x->x!=0, k)
		return i, i, 1
	end
	i = findfirst(x->x!=0, k)
	j = findlast(x->x!=0, k)
	return i,j,2
end

# ╔═╡ 025f13d8-af20-4a57-b7dc-f571bc5c89a6
begin # decompose
	M = zeros(Sym,(5,5))
	for (k,v) in I10_dict
		i,j,f = decompose(k)
		M[i,j] = v / f
		M[j,i] = v / f
	end
	M
end ;

# ╔═╡ c878859c-0046-47e4-a541-eac8b83737f3
Markdown.parse(
"""
```math
\\tiny
\\begin{align}
I &=
"""*
sympy.latex([chains...]')*
"\\cdot"*
sympy.latex(M)*
"\\cdot"*
sympy.latex([chains...])*
"""
\\end{align}
```
""")

# ╔═╡ Cell order:
# ╠═3dbb29a8-1f37-407f-8907-698c80f2bd00
# ╠═04808f8e-25bf-4c6c-9794-5a9d59b197e8
# ╠═f907394e-a7e0-43d0-8394-328c01dcca7c
# ╟─f9f95dab-babb-4190-82af-f1681ce264b4
# ╠═71f9928c-42c2-43c2-899e-6fa756f6b683
# ╠═078cad0e-7202-4ec5-b829-1b44f029160e
# ╠═ca2e840e-df9d-4a27-834d-f6e2f1737378
# ╠═7065b936-6ede-4cd7-9290-8e453d903f48
# ╠═7608afab-8d37-4578-b1ad-7641de820d79
# ╠═444842d4-8b3d-4921-9822-4beb2923fa3f
# ╠═5e7042ca-de8e-4347-93bc-9125fe9562b0
# ╠═905dbfa9-9d6e-4d2b-b048-a72cb75f687e
# ╠═48a50729-332d-4c9d-9b93-b4912f32581d
# ╠═1509915d-4e8f-47d5-8db7-091ab21d799c
# ╠═03f5100d-2b80-437e-9257-be90a44808d8
# ╠═8c23e629-5a24-4c4e-8016-84916badfccf
# ╠═025f13d8-af20-4a57-b7dc-f571bc5c89a6
# ╠═c878859c-0046-47e4-a541-eac8b83737f3
