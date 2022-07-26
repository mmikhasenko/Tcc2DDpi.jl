### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ 5ac3c9c0-0ce8-11ed-102b-03e4d63d6c1b
# ╠═╡ show_logs = false
begin
	import Pkg
		Pkg.add([
			Pkg.PackageSpec("Parameters"),
			Pkg.PackageSpec("LaTeXStrings"),
			Pkg.PackageSpec("SymPy"),
			Pkg.PackageSpec("PyCall")
		])

	using SymPy
	import PyCall
	PyCall.pyimport_conda("sympy.physics.quantum.spin", "sympy")
	import_from(sympy.physics.quantum.spin, (:WignerD, :wignerd), typ=:Any)
	PyCall.pyimport_conda("sympy.physics.wigner",       "sympy")
	import_from(sympy.physics.wigner)
	#
	using LaTeXStrings
	using Parameters
end

# ╔═╡ 04808f8e-25bf-4c6c-9794-5a9d59b197e8
begin
	wignerd(J,λ1,λ2,θ) = WignerD(Sym(J),λ1,λ2,0,θ,0)
	clgd(j1,λ1,j2,λ2,j,λ) =
		clebsch_gordan(Sym(j1),Sym(j2),Sym(j),Sym(λ1),Sym(λ2),Sym(λ))
end

# ╔═╡ f907394e-a7e0-43d0-8394-328c01dcca7c
begin
	@syms(F2::real=>"\\mathcal{F}_2", F3::real=>"\\mathcal{F}_3")
	@syms(θ12::nonnegative=>"\\theta_{12}", θ31::nonnegative=>"\\theta_{31}")
	@syms(ζ23::nonnegative=>"\\zeta^0_{2(3)}")
end ;

# ╔═╡ ca2e840e-df9d-4a27-834d-f6e2f1737378
function amplitude(j0, L, λ)
	j = 1
	S = 1
	sqrt(Sym(2L+1))*sum(
		F2*wignerd(j0,λ,τ,ζ23)*wignerd(j,τ,0,θ31)*clgd(L,0,S,τ,j0,τ) + # chain-2
		(-1)^j * # H_{Dπ} vs H_{πD}
			F3*wignerd(j0,λ,τ,0)*wignerd(j,τ,0,θ12)*clgd(L,0,S,τ,j0,τ) # chain-3
		for τ in -j:j)
end

# ╔═╡ b1c195bf-de94-45b1-9646-533b47edb7b4
 intensity(j0, L) = sum(abs2, amplitude(j0, L, τ) for τ in -j0:j0)

# ╔═╡ 95321b67-fda3-47db-85f9-dfe85562e935
@syms I10=>"I_{1^+}" I01=>"I_{0^-}" I11=>"I_{1^-}" I21=>"I_{2^-}" I22=>"I_{2^+}"

# ╔═╡ c73d57e8-885e-4ef0-985f-d6505cf3b830
Iv = Dict(
	I10 => intensity(1, 0).doit() |> sympy.trigsimp |> expand,
	I01 => intensity(0, 1).doit() |> sympy.trigsimp |> expand,
	I11 => intensity(1, 1).doit() |> sympy.trigsimp |> expand,
	I21 => intensity(2, 1).doit() |> sympy.trigsimp |> expand,
	I22 => intensity(2, 2).doit() |> sympy.trigsimp |> expand
) ;

# ╔═╡ ecac90d0-1682-4816-94e6-eed9350fbe0b
binomialcoeff(e) = sympy.Poly(e, F2, F3).coeffs()

# ╔═╡ 574d4a03-03cd-4a83-a30e-4fc97c6aeb50
function niceprint(e)
	cv = binomialcoeff(e) .|> sympy.trigsimp
	"""
	|$(sympy.latex(F2))|^2 \\left[$(sympy.latex(cv[1]))\\right]+\\\\&\\qquad
	\\text{Re}($(sympy.latex(F2))$(sympy.latex(F3))^*) \\left[$(sympy.latex(cv[2]))\\right]+\\\\&\\qquad
	|$(sympy.latex(F3))|^2 \\left[$(sympy.latex(cv[3]))\\right]
	"""
end ;

# ╔═╡ 69d847b4-c0c6-4f71-ab28-b444bb049602
Markdown.parse(
"""
The unpolarized decay rate reads:
```math
\\small
\\begin{align}
"""*
prod("""
$(sympy.latex(k)) &= $(niceprint(Iv[k])) \\\\\\\\\\\\
""" for k in keys(Iv)) *
"""
\\end{align}
```
""")

# ╔═╡ Cell order:
# ╠═5ac3c9c0-0ce8-11ed-102b-03e4d63d6c1b
# ╠═04808f8e-25bf-4c6c-9794-5a9d59b197e8
# ╠═f907394e-a7e0-43d0-8394-328c01dcca7c
# ╠═ca2e840e-df9d-4a27-834d-f6e2f1737378
# ╠═b1c195bf-de94-45b1-9646-533b47edb7b4
# ╠═95321b67-fda3-47db-85f9-dfe85562e935
# ╠═c73d57e8-885e-4ef0-985f-d6505cf3b830
# ╠═ecac90d0-1682-4816-94e6-eed9350fbe0b
# ╟─574d4a03-03cd-4a83-a30e-4fc97c6aeb50
# ╟─69d847b4-c0c6-4f71-ab28-b444bb049602
