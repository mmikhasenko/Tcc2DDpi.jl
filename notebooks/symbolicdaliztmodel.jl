### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# ╔═╡ 3dbb29a8-1f37-407f-8907-698c80f2bd00
begin
	cd(joinpath(@__DIR__, ".."))
	import Pkg
	Pkg.activate(".")
	Pkg.instantiate()
	# 	Pkg.add([
	# 		Pkg.PackageSpec("LaTeXStrings"),
	# 		Pkg.PackageSpec("SymPy"),
	# 		Pkg.PackageSpec("Parameters"),
	# 		Pkg.PackageSpec("PyCall"),
	# 		Pkg.PackageSpec("JSON")
	# 	])

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

# ╔═╡ 71f9928c-42c2-43c2-899e-6fa756f6b683
abstract type πDD end

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
		ζ23_for0::nonnegative=>"\\zeta^0_{2(3)}",
		ζ23_for1::nonnegative=>"\\zeta^1_{2(3)}")
end ;

# ╔═╡ f9f95dab-babb-4190-82af-f1681ce264b4
md"""
### Amplitude $\pi DD$
"""

@with_kw struct Chain{S}
	j::Int
	j0::Int
	L::Int
	ζ_for0::S
	θ::S
	F::S
end

function amplitude(ch::PPP, λ)
	@unpack F, j, ζ_for0, θ = ch
	@unpack j0, L = ch
	S = j
	A = sqrt(Sym(2L+1)) * F * sum(
			wignerd(j0,λ,τ,ζ_for0)*
				wignerd(j_F,τ,0,θ)*clgd(L,0,S,τ,j0,τ)
			for τ in -j:j)
	return A
end

function amplitude(chains::model, λ)
	A = Sym(0)
	for chain in chains
		A += amplitude(chain, λ)
	end
	return A
end

swaporder(A, from, to; j) = A.subs(from, to*(-1)^j)



# application

@syms(
	F2_Dπ::real=>"\\mathcal{F}_2^{D\\pi}",
	F3_πD::real=>"\\mathcal{F}_3^{\\pi D}",
	F3_Dπ::real=>"\\mathcal{F}_3^{D\\pi}",
	# 
	B2_Dπ::real=>"\\mathcal{B}_3^{D\\pi}",
	B3_πD::real=>"\\mathcal{B}_3^{\\pi D}",
	B3_Dπ::real=>"\\mathcal{B}_3^{D\\pi}",
	# 
	B1_DD::real=>"\\mathcal{B}_3^{DD}"
)

model_πDD_2xDˣ_3xB = [
	PPP(; F=F2_Dπ, j=1, j0, L, ζ=ζ23_for0, θ=θ31),
	PPP(; F=F3_πD, j=1, j0, L, ζ=ζ33_for0, θ=θ12),
	PPP(; F=B3_Dπ, j=0, j0, j, ζ=ζ23_for0, θ=θ31),
	PPP(; F=B2_πD, j=0, j0, L, ζ=ζ23_for0, θ=θ12),
	PPP(; F=B2_DD, j=0, j0, L, ζ=ζ13_for0, θ=θ23)
]

let
	A = amplitude(model_πDD_2xDˣ_3xB, j0, L, λ)
	A_Freplaced = swaporder(A, F3_πD, F3_Dπ; j=1)
	A_FBreplaced = swaporder(A_Freplaced, B3_πD, B3_Dπ; j=0)
	# 
	A_FBreplaced
end







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
		for τ in -j:j)
	# 
	# chain-3
	A += (-1)^j * # H_{Dπ} vs H_{πD} in helicity basis
		sqrt(Sym(2L+1)) * F3 * sum(
			wignerd(j0,λ,τ,ζ33_for0)*wignerd(j_F,τ,0,θ12)*clgd(L,0,S_F,τ,j0,τ)
		for τ in -j:j)
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

# ╔═╡ bb7eecb0-450a-4f66-b781-7df2a75f5359
function amplitude_withNR(ch::πDD, λ, F_NR12, F_NR13, F_NR23)
	@unpack F2, F3 = ch
	j0 = 1 # JP of X are 1+
	L = 0 # X->D D* wave if JP of X are 1+
	j = 1 # J of D*
	S = 1 # S of DD*
	j_NR = 0 # of Dpi
	S_NR = 0 # of DDpi then
	L_NR = 1 # X->D {Dpi} wave if JP of X are 1+
	τ_Dpi = 0
	τ_DD = 0
	sqrt(Sym(2L+1))*sum(
		(
			F2*wignerd(j0,λ,τ,ζ23_for0)*
			wignerd(j,τ,0,θ31)*clgd(L,0,S,τ,j0,τ) + # chain-2
		(-1)^j * # H_{Dπ} vs H_{πD} in helicity basis
			F3*wignerd(j0,λ,τ,0)*wignerd(j,τ,0,θ12)
		)
		*clgd(L,0,S,τ,j0,τ)
		# now adding non-resonant component, S-wave in Dpi
		+ (F_NR13/3*
			wignerd(j0,λ,τ_Dpi,ζ23_for0)*wignerd(j_NR,τ_Dpi,0,θ31)
			+ F_NR12/3*wignerd(j0,λ,τ_Dpi,0)*wignerd(j_NR,τ_Dpi,0,θ12)
		)*clgd(L_NR,0,S_NR,τ_Dpi,j0,τ_Dpi)
		# now S-wave in DD
		+ F_NR23/3
		* wignerd(j0,λ,τ_DD,ζ13_for0)*wignerd(j_NR,τ_DD,0,θ23)
		* clgd(L_NR,0,S_NR,τ_DD,j0,τ_DD) 
		# chain-3
		for τ in -j:j)
end

# ╔═╡ 7065b936-6ede-4cd7-9290-8e453d903f48
intensity(ch::πDD, j0, L) = sum(abs2, amplitude(ch, j0, L, λ) for λ in -j0:j0)

# ╔═╡ c2d22258-0654-4ec5-b9f0-6922c597ee37
intensity_withNR(ch::πDD, F_NR12, F_NR13, F_NR23) = sum(abs2, amplitude_withNR(ch, λ, F_NR12, F_NR13, F_NR23) for λ in -1:1)

# ╔═╡ cf359cc5-28f3-4585-ae29-aa1e81efdd91
F_NR12, F_NR13, F_NR23 = Sym("F^{NR}_{12}"), Sym("F^{NR}_{13}"), Sym("F^{NR}_{23}")

# ╔═╡ 783fc548-adc6-4550-9a52-abe4195bffe9
const tested_j0L = ((1,0), (0,1), (1,1), (2,1), (2,2))

# ╔═╡ 5e7042ca-de8e-4347-93bc-9125fe9562b0
_πDD = πDD(
	SymPy.symbols("\\mathcal{F}_2^{D\\pi}", real = true),
	SymPy.symbols("\\mathcal{F}_3^{D\\pi}", real = true))

# ╔═╡ 359156de-388c-4307-ab50-699821e38834
begin
	IπDD = Dict()
	for (j0,L) in tested_j0L
		p = iseven(L) ? '+' : '-'
		Is = Sym("I_{$(j0)^$(p)}^{(\\pi)}")
		IπDD[Is] = intensity(_πDD, j0, L).doit() |> sympy.trigsimp |> expand
	end
end

# ╔═╡ 905dbfa9-9d6e-4d2b-b048-a72cb75f687e
intensity(_πDD, 1, 0).doit()|> sympy.trigsimp

# ╔═╡ fa453dd3-3507-48d3-9c8e-c2f4077f84b8
intensity_withNR(_πDD, F_NR12, F_NR13, F_NR23).doit() |> sympy.trigsimp |> expand

# ╔═╡ 11a388e8-e959-4f36-ad85-537d1d966eb0
intensity_withNR(_πDD, 0, 0, F_NR23).doit() |> sympy.trigsimp |> expand

# ╔═╡ 145e0c53-f855-4d5c-b423-9c2c2bcf4a7a
(intensity_withNR(_πDD, F_NR12, F_NR13, 0) - intensity_withNR(_πDD, 0, 0, 0) ).doit() |> sympy.trigsimp |> expand

# ╔═╡ ecac90d0-1682-4816-94e6-eed9350fbe0b
binomialcoeff(e, F2, F3) = sympy.Poly(e, F2, F3).coeffs()

# ╔═╡ 0a0dd8c4-b392-4757-90ef-06aadb96ba16
function niceprint(e, F2, F3)
	cv = binomialcoeff(e, F2, F3) .|> sympy.trigsimp .|> sympy.simplify
	"""
	|$(sympy.latex(F2))|^2 \\left[$(sympy.latex(cv[1]))\\right]+\\\\&\\qquad
	\\text{Re}($(sympy.latex(F2))\\overline{$(sympy.latex(F3))})
	\\left[$(sympy.latex(cv[2]))\\right]+\\\\&\\qquad
	|$(sympy.latex(F3))|^2 \\left[$(sympy.latex(cv[3]))\\right]
	"""
end;


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
# ╠═3dbb29a8-1f37-407f-8907-698c80f2bd00
# ╠═71f9928c-42c2-43c2-899e-6fa756f6b683
# ╠═04808f8e-25bf-4c6c-9794-5a9d59b197e8
# ╠═f907394e-a7e0-43d0-8394-328c01dcca7c
# ╠═f9f95dab-babb-4190-82af-f1681ce264b4
# ╠═078cad0e-7202-4ec5-b829-1b44f029160e
# ╠═ca2e840e-df9d-4a27-834d-f6e2f1737378
# ╠═bb7eecb0-450a-4f66-b781-7df2a75f5359
# ╠═7065b936-6ede-4cd7-9290-8e453d903f48
# ╠═c2d22258-0654-4ec5-b9f0-6922c597ee37
# ╠═cf359cc5-28f3-4585-ae29-aa1e81efdd91
# ╠═783fc548-adc6-4550-9a52-abe4195bffe9
# ╠═5e7042ca-de8e-4347-93bc-9125fe9562b0
# ╠═359156de-388c-4307-ab50-699821e38834
# ╠═905dbfa9-9d6e-4d2b-b048-a72cb75f687e
# ╠═fa453dd3-3507-48d3-9c8e-c2f4077f84b8
# ╠═11a388e8-e959-4f36-ad85-537d1d966eb0
# ╠═145e0c53-f855-4d5c-b423-9c2c2bcf4a7a
# ╠═ecac90d0-1682-4816-94e6-eed9350fbe0b
# ╠═0a0dd8c4-b392-4757-90ef-06aadb96ba16
# ╠═69d847b4-c0c6-4f71-ab28-b444bb049602
# ╠═c8bb2ebb-913d-42f3-b5b0-b5d1a81dd6f3
# ╠═88c62736-4631-4521-ba50-85d20b549475
# ╠═0aa313ea-17f7-4ed3-85a1-cdebeb77b4f1
# ╠═97c44e67-e3c6-498b-9275-2c68aecb5fb3
# ╠═9057b18c-40c7-4d99-81a2-50f1af668f59
# ╠═9dcd8295-46cb-4916-9ca5-f7f1706154cf
# ╠═a3ad21c1-e02b-4eed-9f70-10a8d0c7cb58
