### A Pluto.jl notebook ###
# v0.19.9

using Markdown
using InteractiveUtils

# ╔═╡ d8489d70-cf9f-11ec-1f26-0b616e422d52
# ╠═╡ show_logs = false
begin
	import Pkg
	cd(joinpath(@__DIR__,".."))
	Pkg.activate(".")

	using ThreeBodyDecay
	using Plots
	using LaTeXStrings
	using Parameters
	using StaticArrays
end

# ╔═╡ 47666396-b2dd-40c7-afe2-02fa25591963
theme(:wong2, frame=:box, grid=false, minorticks=true,
    guidefontvalign=:top, guidefonthalign=:right,
    xlim=(:auto,:auto), ylim=(:auto,:auto),
    lw=1, lab="", colorbar=false)

# ╔═╡ bb9a04d6-8fe9-4295-9a87-74200c485013
md"""
### Lineshape: standard Breit-Wigners
"""

# ╔═╡ cc71b82b-a1fa-43ff-879a-56aadda3b483
function BlattWeiss(z²,l)
	l==0 && return 1.0
	l==1 && return 1/(1+z²)
	l==2 && return 1/(9+3z²+z²^2)
	error("Not implemented, l=$l")
end



# ╔═╡ e357f507-0730-4a4b-9c5e-689adf4734be
begin
	const mπ = 0.13957039
	const mD = 1.86483
	const mDˣ, ΓDˣ = 1.86483+145.4258e-3, 83.4e-6
	const δmTcc = -300e-6 # 300keV
	const mTcc, ΓTcc = mD+mDˣ+δmTcc, 50e-6
end ;

# ╔═╡ 11bd1e32-a909-4971-88c0-ca36d4f7942a
begin
	@with_kw struct RBW
		m::Float64
		mi::Float64
		mj::Float64
		Γ::Float64
		l::Int = 0
		R::Float64 = 1.5
	end
	(lineshape::RBW)(σ::Float64) = amplitude(lineshape, σ)
	# 
	import ThreeBodyDecay: amplitude
	function amplitude(lineshape::RBW, σ::Float64)
		@unpack m, Γ, mi, mj, l, R = lineshape
		sqrtσ = sqrt(σ)
		p = sqrtKallenFact(sqrtσ,mi,mj)/(2sqrtσ)
		p0 = sqrtKallenFact(m,mi,mj)/(2m)
		# 
		ff² = BlattWeiss((p*R)^2,l) / BlattWeiss((p0*R)^2,l)
		Γ_dep = Γ*ff²*(p/p0)^(2l+1)*m/sqrtσ
		return p^l*sqrt(ff²)*m*Γ/(m^2-σ-1im*m*Γ_dep)
	end
	@with_kw struct X3b{T}
		lineshape::T
		m0::Float64
		mk::Float64
		L::Int = 0
		r0::Float64 = 1.5
	end
	function (x::X3b)(σ::Float64)
		@unpack lineshape, m0, mk, r0, L = x
		sqrtσ = sqrt(σ)
		q = sqrtKallenFact(m0,sqrtσ,mk)/(2m0)
		ff² = BlattWeiss((q*r0)^2,L)
		return q^L*sqrt(ff²)*amplitude(lineshape, σ)
	end
	# 
	const Wave = NamedTuple{(:lineshape3b, :j0, :P0, :k, :lineshape, :j, :parity)}
	ms(m0) = ThreeBodyMasses(mπ,mD,mD; m0)
	# 
	function wave2decaychain(wave::Wave, s::Float64)
		@unpack lineshape3b, j0, P0 = wave
		@unpack k, lineshape, j, parity = wave
		# 
		m0 = sqrt(s)
		_ms = ms(m0)
		tbs = ThreeBodySystem(; ms=_ms, two_js=(0,0,0,j0) .|> x2)
		Ps = ('-','-','-',P0)
		d = DecayChainLS(2, lineshape; two_s=2j, parity, Ps, tbs)
		Xlineshape = X3b(; lineshape, m0, mk=_ms[k], L=div(d.HRk.two_ls[1],2))
		return DecayChainLS(k, Xlineshape; two_s=2j, parity, Ps, tbs)
	end
	function amplitude(wave::Wave, s::Float64, σs, λ::Int)
		@unpack lineshape3b = wave
		dX = wave2decaychain(wave, s)
		return Osym(dX,σs,λ)*amplitude(lineshape3b,s)
	end
end

# ╔═╡ 3e672231-c34c-4c99-8d68-2802607a20cc
function Osym(dc::T, σs, λ::Int) where T<:DecayChain
	dc.k == 1 && return amplitude(dc, σs, (0,0,0,2λ))
	dc.k == 2 && return amplitude(dc, σs, (0,0,0,2λ)) +
			minusone()^div(dc.two_s,2) * 
				amplitude(DecayChain(dc, k=3), σs, (0,0,0,2λ))
	error("Unaccounted case, k=$(dc.k)")
end

# ╔═╡ ceb3339d-6fbd-4643-8d46-da2abd6e1af7
BW_Dˣ = RBW(m=mDˣ, Γ=ΓDˣ, mi=mD, mj=mπ, l=1)

# ╔═╡ b40114c5-a04f-42f4-9fe0-8ecff826e4c7
BW_Tcc = RBW(m=mTcc, Γ=ΓTcc, mi=2mD, mj=mπ, l=0)

# ╔═╡ 01ea0651-7afc-40d1-95c5-5e4e91b9058e
begin
	struct model{N,T<:Wave}
		chains::StaticVector{N,T}
		couplings::StaticVector{N,ComplexF64}
	end
	model(dv::AbstractArray, cv::AbstractArray) =
		(N = length(dv); model(SVector{N}(dv), SVector{N}(cv .* (1+0im))))
end

# ╔═╡ 92c48924-1bb9-4eb4-a860-ccb80bf1dbcd
md"""
### An example
"""

# ╔═╡ 3c56835d-16f1-4958-89ed-02ebba3b15b4
wave(_jp::jp) = (lineshape3b=BW_Tcc, j0=_jp.j, P0=_jp.p, k=2, lineshape=BW_Dˣ, j=1, parity='-')

# ╔═╡ 40473fac-21c4-4be7-8a3e-99fbede99b86
const jMax = 1

# ╔═╡ 88b7a1a1-f787-413f-a891-fa28f53b9378
intensity(w::Wave, s, σs) = sum(abs2, 
		amplitude(w, s, σs, λ) for λ in -jMax:jMax)

# ╔═╡ 89debcdf-8e45-4fd1-999c-0cfcde36a224
md"""
### Computation on a grid
"""

# ╔═╡ c4bc1803-12c7-443a-a14a-684a32e6c7a8
plot(ms(mTcc), σs->intensity(wave(str2jp("1+")), mTcc^2, σs), iσx=2, iσy=3)

# ╔═╡ 2512c76b-8875-429c-977b-714b2f1d1405
plot(
	[plot(
	ms(mTcc), σs->intensity(wave(str2jp(_jp)), mTcc^2, σs), iσx=2, iσy=3, title="jp=$_jp")
	for _jp in ["1+", "0-", "2-", "1-"]]..., layout=grid(2,2))

# ╔═╡ Cell order:
# ╠═d8489d70-cf9f-11ec-1f26-0b616e422d52
# ╠═47666396-b2dd-40c7-afe2-02fa25591963
# ╟─bb9a04d6-8fe9-4295-9a87-74200c485013
# ╠═cc71b82b-a1fa-43ff-879a-56aadda3b483
# ╠═11bd1e32-a909-4971-88c0-ca36d4f7942a
# ╠═3e672231-c34c-4c99-8d68-2802607a20cc
# ╠═e357f507-0730-4a4b-9c5e-689adf4734be
# ╠═ceb3339d-6fbd-4643-8d46-da2abd6e1af7
# ╠═b40114c5-a04f-42f4-9fe0-8ecff826e4c7
# ╠═01ea0651-7afc-40d1-95c5-5e4e91b9058e
# ╟─92c48924-1bb9-4eb4-a860-ccb80bf1dbcd
# ╠═3c56835d-16f1-4958-89ed-02ebba3b15b4
# ╠═40473fac-21c4-4be7-8a3e-99fbede99b86
# ╠═88b7a1a1-f787-413f-a891-fa28f53b9378
# ╟─89debcdf-8e45-4fd1-999c-0cfcde36a224
# ╠═c4bc1803-12c7-443a-a14a-684a32e6c7a8
# ╠═2512c76b-8875-429c-977b-714b2f1d1405
