### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ 5ac3c9c0-0ce8-11ed-102b-03e4d63d6c1b
# ╠═╡ show_logs = false
begin
	import Pkg
		Pkg.add([
			Pkg.PackageSpec("LaTeXStrings"),
			Pkg.PackageSpec("SymPy"),
			Pkg.PackageSpec("Parameters"),
			Pkg.PackageSpec("PyCall"),
			Pkg.PackageSpec("JSON")
		])

	using SymPy
	import PyCall
	PyCall.pyimport_conda("sympy.physics.quantum.spin", "sympy")
	import_from(sympy.physics.quantum.spin, (:WignerD, :wignerd), typ=:Any)
	PyCall.pyimport_conda("sympy.physics.wigner",       "sympy")
	import_from(sympy.physics.wigner)
	#
	using Parameters
	using LaTeXStrings
	using JSON
end

# ╔═╡ 04808f8e-25bf-4c6c-9794-5a9d59b197e8
begin
	wignerd(J,λ1,λ2,θ) = WignerD(Sym(J),λ1,λ2,0,θ,0)
	clgd(j1,λ1,j2,λ2,j,λ) =
		clebsch_gordan(Sym(j1),Sym(j2),Sym(j),Sym(λ1),Sym(λ2),Sym(λ))
end

# ╔═╡ f907394e-a7e0-43d0-8394-328c01dcca7c
begin
	@syms(
		θ12::nonnegative=>"\\theta_{12}",
		θ31::nonnegative=>"\\theta_{31}")
	@syms(
		ζ23_for0::nonnegative=>"\\zeta^0_{2(3)}", ζ23_for1::nonnegative=>"\\zeta^1_{2(3)}")
end ;

# ╔═╡ f9f95dab-babb-4190-82af-f1681ce264b4
md"""
### Amplitude $\pi DD$
"""

# ╔═╡ 078cad0e-7202-4ec5-b829-1b44f029160e
struct πDD{T}
	F2::T
	F3::T
end

# ╔═╡ ca2e840e-df9d-4a27-834d-f6e2f1737378
function amplitude(ch::πDD, j0, L, λ)
	@unpack F2, F3 = ch
	j = 1
	S = 1
	sqrt(Sym(2L+1))*sum(
		F2*wignerd(j0,λ,τ,ζ23_for1)*
			wignerd(j,τ,0,θ31)*clgd(L,0,S,τ,j0,τ) + # chain-2
		(-1)^j * # H_{Dπ} vs H_{πD} in helicity basis
			F3*wignerd(j0,λ,τ,0)*wignerd(j,τ,0,θ12)*clgd(L,0,S,τ,j0,τ) # chain-3
		for τ in -j:j)
end

# ╔═╡ 912c2553-6197-4296-b429-429712beadb2
md"""
### Amplitude $\gamma DD$
"""

# ╔═╡ 2b07445f-0710-45b7-9127-6582ea86f2b5
struct γDD{T}
	F2::T
	F3::T
end

# ╔═╡ 2af36530-be75-41ec-8c86-0e86abeb441c
function amplitude(ch::γDD, j0, L, λ, ρ)
	@unpack F2, F3 = ch
	# Tcc → Dˣ D
	j = 1 # Dˣ
	S = 1 # Dˣ⊗D
	# 
	# Dˣ → Dπ P-wave
	jγ = 1 
	lγ = 1
	sγ = 1 # γ⊗D
	# 
	sqrt(Sym(2L+1))*sum(
		# 
		# chain-2: (D₃γ)D₂
			F2*wignerd(j0,λ,τ,ζ23_for0) *
		(-1)^(jγ-ρ′) * # particle-2 phase
		clgd(L,0,S,τ,j0,τ)*wignerd(j,τ,-ρ′,θ31)*clgd(lγ,0,sγ,-ρ′,j,-ρ′) *
			wignerd(jγ,ρ′,ρ,ζ23_for1) - # isospin minus
		# 
		# chain-3: (γD₂)D₃
		(-1)^(lγ+sγ-jγ) * # H_{Dγ} vs H_{γD} in ls basis. Martin-Spearman (5.57)
			F3*wignerd(j0,λ,τ,0)*
		clgd(L,0,S,τ,j0,τ)*wignerd(j,τ,ρ,θ12)*clgd(lγ,0,sγ,ρ′,j,ρ′) * 
			wignerd(jγ,ρ′,ρ,0)
		# 
		for τ in -j:j, ρ′ in -jγ:jγ)
end

# ╔═╡ b1c195bf-de94-45b1-9646-533b47edb7b4
 intensity(ch::πDD, j0, L) = sum(abs2, amplitude(ch, j0, L, λ) for λ in -j0:j0)

# ╔═╡ 81b43385-047a-4d33-83ce-f72beed5b87e
"""
	intensity(ch::γDD, j0, L)

The amplitude is squared and summed over the Tcc and γ helicity states.
The summation is implemented including all states,
```
for λ in -j0:j0, ρ in -1:1
```
howereve, ρ=0 should be excluded.
It is excluded effectively by the Dˣ to D γ decay Clebch-Gardan coefficient,
```
clgd(lγ,0,sγ,ρ′,j,ρ′) = 0 for ρ′=0
```
The Kroniker delta for ρ′ and ρ′ is given by
 - `wignerd(jγ,ρ′,ρ,0)` for the chain 3,
 - `wignerd(jγ,ρ′,ρ,ζ23_for1)` for the chain 2 since `ζ23_for1 = 0 or pi` for the zero mass particle.
"""
intensity(ch::γDD, j0, L) = sum(abs2, amplitude(ch, j0, L, λ, ρ)
	for λ in -j0:j0, ρ in -1:1)

# ╔═╡ e1b3df12-2085-491e-ab7e-9a015178757e
md"""
## Symbolic computation
"""

# ╔═╡ 6b7787e2-d751-422b-b6cd-a5d4532a3b31
const tested_j0L = ((1,0), (0,1), (1,1), (2,1), (2,2))

# ╔═╡ 01c3d283-8e67-4590-9b15-687685b80394
_πDD = πDD(
	SymPy.symbols("\\mathcal{F}_2^{D\\pi}", real = true),
	SymPy.symbols("\\mathcal{F}_3^{D\\pi}", real = true))

# ╔═╡ a3096302-b8e9-4d7b-8ee6-b1059401461e
begin
	IπDD = Dict()
	for (j0,L) in tested_j0L
		p = iseven(L) ? '+' : '-'
		Is = Sym("I_{$(j0)^$(p)}^{(\\pi D D)}")
		IπDD[Is] = intensity(_πDD, j0, L).doit() |> sympy.trigsimp |> expand
	end
end

# ╔═╡ 987d81a4-ccff-476e-9322-a325d1bdf04a
_γDD = γDD(
	SymPy.symbols("\\mathcal{F}_2^{D\\gamma}", real = true),
	SymPy.symbols("\\mathcal{F}_3^{D\\gamma}", real = true))

# ╔═╡ fb0a7eab-7d6d-441e-b17c-b14810e8e569
begin
	IγDD = Dict()
	for (j0,L) in tested_j0L
		p = iseven(L) ? '+' : '-'
		Is = Sym("I_{$(j0)^$(p)}^{(\\gamma D D)}")
		IγDD[Is] = intensity(_γDD, j0, L).doit() |> sympy.trigsimp |> expand
	end
end

# ╔═╡ ecac90d0-1682-4816-94e6-eed9350fbe0b
binomialcoeff(e, F2, F3) = sympy.Poly(e, F2, F3).coeffs()

# ╔═╡ 574d4a03-03cd-4a83-a30e-4fc97c6aeb50
function niceprint(e, F2, F3)
	cv = binomialcoeff(e, F2, F3) .|> sympy.trigsimp .|> sympy.simplify
	"""
	|$(sympy.latex(F2))|^2 \\left[$(sympy.latex(cv[1]))\\right]+\\\\&\\qquad
	\\text{Re}($(sympy.latex(F2))\\overline{$(sympy.latex(F3))}) \\left[$(sympy.latex(cv[2]))\\right]+\\\\&\\qquad
	|$(sympy.latex(F3))|^2 \\left[$(sympy.latex(cv[3]))\\right]
	"""
end ;

# ╔═╡ 69d847b4-c0c6-4f71-ab28-b444bb049602
Markdown.parse(
"""
The unpolarized decay rate to πDD reads:
```math
\\tiny
\\begin{align}
"""*
prod("""
$(sympy.latex(k)) &= $(niceprint(IπDD[k], _πDD.F2, _πDD.F3)) \\\\\\\\\\\\
""" for k in keys(IπDD)) *
"""
\\end{align}
```
""")

# ╔═╡ c8bb2ebb-913d-42f3-b5b0-b5d1a81dd6f3
Markdown.parse(
"""
The unpolarized decay rate into γDD reads:
```math
\\tiny
\\begin{align}
"""*
prod("""
$(sympy.latex(k)) &= $(niceprint(IγDD[k], _γDD.F2, _γDD.F3))) \\\\\\\\\\\\
""" for k in keys(IγDD)) *
"""
\\end{align}
```
""")

# ╔═╡ 88c62736-4631-4521-ba50-85d20b549475
md"""
### Convert sympy to c code
"""

# ╔═╡ 0aa313ea-17f7-4ed3-85a1-cdebeb77b4f1
function expr2ccode(exp, extrasubs)
	exp′ = exp.subs(extrasubs)
	exp′′ = exp′.subs(
			Dict(
				θ12=>Sym("theta12"),
				θ31=>Sym("theta31"),
				ζ23_for1=>Sym("zeta23_for1"),
				ζ23_for0=>Sym("zeta23_for0")))
	exp′′′ = exp′′.subs(Dict(
			Sym("F2")^2 => Sym("F2abs2"),
			Sym("F3")^2 => Sym("F3abs2"),
			Sym("F2")*Sym("F3")=>Sym("ReF2F3x")))
	sympy.ccode(exp′′′)
end

# ╔═╡ 97c44e67-e3c6-498b-9275-2c68aecb5fb3
begin
	codes = Dict()
	for k in keys(IγDD)
		extrasubs = Dict(
				_γDD.F2=>Sym("F2"),
				_γDD.F3=>Sym("F3"))
		codes[k] = expr2ccode(IγDD[k], extrasubs)
	end
	for k in keys(IπDD)
		extrasubs = Dict(
				_πDD.F2=>Sym("F2"),
				_πDD.F3=>Sym("F3"))
		codes[k] = expr2ccode(IπDD[k], extrasubs)
	end
end

# ╔═╡ 9057b18c-40c7-4d99-81a2-50f1af668f59
md"""
### Write to json
"""

# ╔═╡ 9dcd8295-46cb-4916-9ca5-f7f1706154cf
function writejson(path, obj)
    open(path, "w") do io
        JSON.print(io, obj, 4)
    end
end

# ╔═╡ a3ad21c1-e02b-4eed-9f70-10a8d0c7cb58
writejson("code_gDD.json", codes)

# ╔═╡ Cell order:
# ╠═5ac3c9c0-0ce8-11ed-102b-03e4d63d6c1b
# ╠═04808f8e-25bf-4c6c-9794-5a9d59b197e8
# ╠═f907394e-a7e0-43d0-8394-328c01dcca7c
# ╟─f9f95dab-babb-4190-82af-f1681ce264b4
# ╠═078cad0e-7202-4ec5-b829-1b44f029160e
# ╠═ca2e840e-df9d-4a27-834d-f6e2f1737378
# ╠═b1c195bf-de94-45b1-9646-533b47edb7b4
# ╟─912c2553-6197-4296-b429-429712beadb2
# ╠═2b07445f-0710-45b7-9127-6582ea86f2b5
# ╠═2af36530-be75-41ec-8c86-0e86abeb441c
# ╠═81b43385-047a-4d33-83ce-f72beed5b87e
# ╟─e1b3df12-2085-491e-ab7e-9a015178757e
# ╠═6b7787e2-d751-422b-b6cd-a5d4532a3b31
# ╠═01c3d283-8e67-4590-9b15-687685b80394
# ╠═a3096302-b8e9-4d7b-8ee6-b1059401461e
# ╠═987d81a4-ccff-476e-9322-a325d1bdf04a
# ╠═fb0a7eab-7d6d-441e-b17c-b14810e8e569
# ╠═ecac90d0-1682-4816-94e6-eed9350fbe0b
# ╠═574d4a03-03cd-4a83-a30e-4fc97c6aeb50
# ╟─69d847b4-c0c6-4f71-ab28-b444bb049602
# ╟─c8bb2ebb-913d-42f3-b5b0-b5d1a81dd6f3
# ╠═88c62736-4631-4521-ba50-85d20b549475
# ╠═0aa313ea-17f7-4ed3-85a1-cdebeb77b4f1
# ╠═97c44e67-e3c6-498b-9275-2c68aecb5fb3
# ╟─9057b18c-40c7-4d99-81a2-50f1af668f59
# ╠═9dcd8295-46cb-4916-9ca5-f7f1706154cf
# ╠═a3ad21c1-e02b-4eed-9f70-10a8d0c7cb58
