### A Pluto.jl notebook ###
# v0.14.0

using Markdown
using InteractiveUtils

# This Pluto notebook uses @bind for interactivity. When running this notebook outside of Pluto, the following 'mock version' of @bind gives bound variables a default value (instead of an error).
macro bind(def, element)
    quote
        local el = $(esc(element))
        global $(esc(def)) = Core.applicable(Base.get, el) ? Base.get(el) : missing
        el
    end
end

# ╔═╡ d745eee9-7229-4767-8fdc-ea1ecda3e61a
begin
	using PlutoUI
	using Test
	using QuadGK
	using Cuba
	using Parameters # not used
	using Interpolations
	using Optim
	using ThreeBodyDecay
	using Plots
	using Measurements
	theme(:wong2, frame=:box, grid=false, minorticks=true)
end

# ╔═╡ 1919aca2-8900-4736-8119-bc1dd7b8a7e4
using LaTeXStrings

# ╔═╡ b3a386e0-ed2e-4ed0-aa91-725fbad987f4
md"""
# Pole position of three-body model
""" 

# ╔═╡ 6995af97-a29f-4931-ab94-85650d0e48aa
md"""
masses
"""

# ╔═╡ 9052fdd5-92a2-4ca4-bcea-726e9a3df19a
begin
	mπ⁰, mπ⁺ = 0.1349768, 0.13957039
	mD⁰, mD⁺ = 1.86483, 1.86965
	mDˣ⁺ = 2.01026; ΓDˣ⁺ = 83.4e-6
	mDˣ⁰ = 2.00685; ΓDˣ⁰ = 55.2e-6
	mγ = 0.0
end

# ╔═╡ 976f5e93-5c0c-4ba4-9e87-b53e0d87ebee
md"""
couplings
"""

# ╔═╡ 14c58917-3204-43fd-bc62-a54f5481d127
begin
	const h² = 20.13e-3
	const f² = 282.42
	const μ₀ = -3.77
	const μ₊ = 1.0
end

# ╔═╡ 27ffd025-94e5-40aa-9910-7d766a724dd8
e2m(e) = (mD⁰+mDˣ⁺)+e*1e-3

# ╔═╡ f47f2687-3c7a-4127-9675-97be5c22f8c6
m2e(m) = (m-mD⁰-mDˣ⁺)*1e3

# ╔═╡ f825bbfc-a0a4-47e4-9b4e-5da5bc7fe724
md"""
fit results
"""

# ╔═╡ 8a551387-0acb-4aa4-9825-4f46de024492
δm_sw, δm_ΔM = (-363 ± 40)*1e-3, (-375 ± 40)*1e-3;

# ╔═╡ 190be781-8951-4ce8-b06f-6e621c609527
md"""
### $D^*$ propagator
"""

# ╔═╡ b6284e99-49f4-4513-b478-ee4bdad3ae9c
md"""
## General for all channels
"""

# ╔═╡ d808ac86-4893-4927-a120-68cdcd923573
md"""
### Integration over Dalitz
"""

# ╔═╡ c88d9ff1-699a-40ab-b96e-98732de40754
begin
	import Base.^
	^(ms::NamedTuple{(:m1,:m2,:m3),T} where T, n::Int) = Tuple(ms).^n
end

# ╔═╡ 82df0901-7a5c-496b-a090-c18bce6ba755
function J_I(σ,pars)
	m, Γ = pars
	FF = 1.0 # (σ-(mD⁰+mπ⁺)^2) / (m^2-(mD⁰+mπ⁺)^2)
	1/(m^2 - σ + 1im*m*Γ*FF)
end

# ╔═╡ 31268d53-342f-4a29-9219-c13ea17ba418
function J_II(σ,pars)
	m, Γ = pars
	FF = 1.0 # (σ-(mD⁰+mπ⁺)^2) / (m^2-(mD⁰+mπ⁺)^2)
	1/(m^2 - σ - 1im*m*Γ*FF)
end

# ╔═╡ 035f6347-7dc7-45e4-83a0-116359128ef4
pole(R::NamedTuple{(:m,:Γ),T} where T) = R.m^2-1im*R.m*R.Γ

# ╔═╡ 6ae3dc60-696f-489c-ba3d-0969381c3605
md"""
### Interpolation, merging 3-body and 2-body
"""

# ╔═╡ 58176be4-6300-4e5b-9543-73d7af40ebf7
md"""
### Dispersive integrals
"""

# ╔═╡ 5cc68d15-0242-4aa0-ab14-69fe07df5366
md"""
### Pole search of amplitude
"""

# ╔═╡ 29c0d772-6248-4bb0-9c7e-e3d286a92b4f
md"""
For the second sheet - add discontinuity
"""

# ╔═╡ 8c4612df-7d5b-4ac4-a29c-36591225c12e
md"""
to call the denominator on the sum (tuple) of dispersive terms
"""

# ╔═╡ 80b3743a-cb41-4edd-9319-9e16530f2e3a
begin
	denominator_I(ds::Tuple,e,δm) = sum(denominator_I(d,e,δm) for d in ds)
	denominator_II(ds::Tuple,e,δm) = sum(denominator_II(d,e,δm) for d in ds)
end;

# ╔═╡ ba3768a8-17e8-4e9f-b496-a614a028a84e
md"""
## Application to $X\to \pi(\gamma) D D$
three channels
"""

# ╔═╡ 5b9360b1-32ac-4827-9046-ea19ab785f44
md"""
### Scattering length
""" 

# ╔═╡ 36e62ff2-28d5-471d-9a9e-2c0d531ef540
savefig(joinpath("inverse_matching.pdf"))

# ╔═╡ f45a5b27-9c61-4982-a2ef-287282b02d47
md"""
### Pole position for the fit
"""

# ╔═╡ e634d29d-46eb-464a-abea-d4fd30f228aa
δmv = -0.8:0.1:-0.1

# ╔═╡ 4c4cafe9-d0a8-4daa-a475-287756a6014c
md"""
### Second solution
"""

# ╔═╡ d2a66d89-1965-4b76-bee0-f55dbe787722
md"""
### Intensity shape in the $D^0D^0\pi^+$ spectrum
"""

# ╔═╡ 5a685e9f-fa97-479f-a31a-fa7cea7c82c3
@bind δm0 Slider(-0.3:0.1:3, default=2.0, show_value=true)

# ╔═╡ 68caab0a-81ca-4234-9460-bea389ef2ab4
# savefig(joinpath("advanced_pole.pdf"))

# ╔═╡ d12ad72e-274d-43c0-96a2-ade860abebb9
md"""
### Checks
"""

# ╔═╡ 492f23d1-cd3b-46cd-937e-908c56e19674
### Plotting recipies
@recipe function f(f::Function, xlims::Tuple, ylims::Tuple;
		Nx=20, Ny=20, post=log10∘abs2)
	xv = range(xlims..., length=Nx) # MeV
	yv = range(ylims..., length=Ny) # MeV
	ev = xv' .+ 1im .* yv
	zv = post.(f.(ev))
	@series begin
		label := ""
		seriestype := :path
		linecolor := :red
		linewidth := 2
		([xlims...], [0.0, 0.0])
	end
	xlims --> xlims
	ylims --> ylims
	seriestype --> :contour
	colorbar --> false # title="I/II sh."
	(xv,yv,zv)
end

# ╔═╡ e7883cde-4609-4385-91f7-6b408d855bb1
@bind x Slider(0:0.01:1; show_value=true, default=0.3)

# ╔═╡ 54f8c1e0-7cec-4341-8210-338d8874b50c
md"""
### Covariant expressions
"""

# ╔═╡ b895e1e3-a65d-479a-bc8b-f794e8a98cdd
λ(x,y,z) = x^2+y^2+z^2-2x*y-2y*z-2z*x

# ╔═╡ d26c2516-1f71-42e1-b3ff-3004d1d144ab
begin
	abstract type AbstractxDD end
	#
	function ρ_tb(d::AbstractxDD, e)
		sqrt(λ(e2m(e)^2,mDˣ⁺^2,mD⁰^2))/e2m(e)^2
	end	
		
	# πDD
	struct πDD{T1,T2,T3} <: AbstractxDD
		ms::NamedTuple{(:m1,:m2,:m3),T1}
		R12::NamedTuple{T2}
		R13::NamedTuple{T3}
	end
	
	# γDD
	struct γDD{T1,T2,T3} <: AbstractxDD
		ms::NamedTuple{(:m1,:m2,:m3),T1}
		R12::NamedTuple{T2}
		R13::NamedTuple{T3}
	end
end

# ╔═╡ 0355ece1-5558-4213-961f-74e115342ea0
branch_points(d::AbstractxDD) = (
	m2e(d.ms[3] + sqrt(pole(d.R12))),
	m2e(d.ms[2] + sqrt(pole(d.R13))))

# ╔═╡ fbd283a8-f291-481d-98cb-941371717b0f
begin
	channels = [
		πDD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), (m=mDˣ⁺,Γ=ΓDˣ⁺), (m=mDˣ⁺,Γ=ΓDˣ⁺)),
		πDD((m1=mπ⁰,m2=mD⁺,m3=mD⁰), (m=mDˣ⁺,Γ=ΓDˣ⁺), (m=mDˣ⁰,Γ=ΓDˣ⁰)),
		γDD((m1=mγ, m2=mD⁺,m3=mD⁰), (m=mDˣ⁺,Γ=ΓDˣ⁺), (m=mDˣ⁰,Γ=ΓDˣ⁰))]
end;

# ╔═╡ 30748ddd-9757-461f-b95e-c409e1f81718
σ2of3_pm(σ,msq,s) = (msq[1]+msq[3]+
	(σ+msq[1]-msq[2])*(s-σ-msq[3])/(2*σ)) .+ [-1,1] .*
		sqrt(λ(σ,msq[1],msq[2])*λ(s,σ,msq[3]))/(2*σ)

# ╔═╡ 8fe15b8e-4ffb-4aae-bae6-46a71115caed
let ms = channels[1].ms
	msq = ms^2
	s = e2m(-0.12-0.01im)^2
# 	
	σ3(y) = (ms[1]+ms[2])^2 + ((sqrt(s)-ms[3])^2-(ms[1]+ms[2])^2) * y
	σ3v = σ3.(range(0,1, length=100))
	σ2v = σ2of3_pm.(σ3v,Ref(msq),s)
# 	
	plot(σ3v, lab="s12")
	plot!(getindex.(σ2v, 1), lab="s13")
	scatter!(σ2of3_pm(σ3(x),msq,s), lab="s13⁺,s13⁻", mc=:black)
end

# ╔═╡ d93cc18a-8b23-4a0e-9743-febf5e41f176
begin
	struct NonRelBW
	end
	function denominator_I(d::NonRelBW, e, δm)
		ik(e) = 1im*sqrt(λ(e2m(e)^2,mDˣ⁺^2,mD⁰^2))/(2*e2m(e))
		ike = ik(e)
		ik0 = ik(δm+1e-6im)
		ike - real(ik0)
	end
end

# ╔═╡ 218df0f7-36ed-4f12-b285-ee675c0dffeb
s23(v) = v.s+v.msq[1]+v.msq[2]+v.msq[3]-v.s12-v.s13

# ╔═╡ 8b959ab0-adb7-11eb-179d-77286258f875
begin
	p_p12( v) = (v.s+v.s12-v.msq[3])/2
	p12_p2(v) = (v.s12+v.msq[2]-v.msq[1])/2
	p_p2(  v) = (v.s+v.msq[2]-v.s13)/2
	# 
	p_p13( v) = (v.s+v.s13-v.msq[2])/2
	p13_p3(v) = (v.s13+v.msq[3]-v.msq[1])/2
	p_p3(  v) = (v.s+v.msq[3]-v.s12)/2
	#
	p1_p2( v) = (v.s12-v.msq[1]-v.msq[2])/2
	p1_p3( v) = (v.s13-v.msq[1]-v.msq[3])/2
	p2_p3( v) = (s23(v)-v.msq[2]-v.msq[3])/2
	#
	p12_p13(v) = (2*(v.s+v.msq[1])-s23(v)-v.s12-v.s13)/2
# 	
	p12_p3(v) = (v.s-v.s12-v.msq[3])/2
	p13_p2(v) = (v.s-v.s13-v.msq[2])/2
end;

# ╔═╡ 52589f04-a93b-4b61-a6c7-1fd97b4e8a4f
begin
	A(v) = λ(v.s12,v.msq[1],v.msq[2])/(4*v.s12) +
		(p_p12(v) * p12_p2(v)/v.s12 - p_p2(v))^2 / v.s
	B(v) = λ(v.s13,v.msq[1],v.msq[3])/(4*v.s13) +
		(p_p13(v) * p13_p3(v)/v.s13 - p_p3(v))^2 / v.s
	C(v) = D(v)+E(v)
	D(v) = 
		p12_p3(v)*p12_p2(v)/v.s12 +
		p13_p3(v)*p13_p2(v)/v.s13 -
		p12_p13(v)*p12_p2(v)*p13_p3(v)/(v.s12*v.s13) -
		p2_p3(v);
	E(v) =
		p_p12(v)*p_p13(v)*
			p12_p2(v)*p13_p3(v)/(v.s*v.s12*v.s13) +
		p_p2(v)*p_p3(v) / v.s -
		p_p12(v)*p_p3(v)*p12_p2(v)/(v.s12*v.s) -
		p_p13(v)*p_p2(v)*p13_p3(v)/(v.s13*v.s)
	G(v) = (
		2*p1_p2(v)*p1_p3(v)*p2_p3(v) -
		v.msq[2]^2*p1_p3(v)^2 -
		v.msq[3]^2*p1_p2(v)^2) / (2*v.s)
end;

# ╔═╡ 10eccd1f-3956-4cd8-a416-4da8bf384c12
function decay_matrix_element_squared(d::πDD,s,σ3,σ2)
	msq = d.ms^2
	v = (;s,s12=σ3,s13=σ2,msq)
# 	
	J12_I, J12_II = J_I(σ3,d.R12), J_II(σ3,d.R12)
	J13_I, J13_II = J_I(σ2,d.R13), J_II(σ2,d.R13)
# 	
	frakM = A(v) * J12_I * J12_II +
			B(v) * J13_I * J13_II +
			C(v) * (J13_I * J12_II +  J12_I * J13_II)
	f²*frakM/3/4
end

# ╔═╡ 581a83d6-62ea-4d7d-8728-4fddfb533409
function decay_matrix_element_squared(d::γDD,s,σ3,σ2)
	msq = d.ms^2
	v = (;s,s12=σ3,s13=σ2,msq)
# 	
	J12_I, J12_II = J_I(σ3,d.R12), J_II(σ3,d.R12)
	J13_I, J13_II = J_I(σ2,d.R13), J_II(σ2,d.R13)
# 	
	_p1_p2 = p1_p2(v)
	_p1_p3 = p1_p3(v)
	_G = G(v)
# 	
	frakM = μ₊^2 * (_p1_p2^2+_G) * J12_I * J12_II +
			μ₀^2 * (_p1_p3^2+_G) * J13_I * J13_II -
			μ₊*μ₀ * (_p1_p2*_p1_p3 - _G) *
		(J13_I * J12_II +  J12_I * J13_II)
# 	
	h²*frakM/3
end

# ╔═╡ 2b852b33-01da-4724-8294-9ebbcefa03b4
function integrand_mapped_thr(d::AbstractxDD,s,x)
	# 	
	σ3_0, σ3_e = (d.ms[1]+d.ms[2])^2, (√s-d.ms[3])^2
	σ3 = σ3_0 + x[1]*(σ3_e-σ3_0) # straight path
	# 		
	σ2_0, σ2_e = σ2of3_pm(σ3, d.ms^2, s)
	σ2 = σ2_0 + x[2]*(σ2_e-σ2_0) # straight path
	#
	jac = (σ3_e-σ3_0)*(σ2_e-σ2_0)
	decay_matrix_element_squared(d,s,σ3,σ2) / s * jac
end

# ╔═╡ 49058565-1ca9-409c-b3f1-5a5739610d4e
function ρ_thr(d::AbstractxDD, e)
	integrand(x,f) = f[1:2] .= reim(integrand_mapped_thr(d,e2m(e)^2,x))
	v = cuhre(integrand, 2, 2)[1]
	complex(v...)
end

# ╔═╡ fec685a1-db09-4913-bb7f-ac38eeda9942
begin
	struct interpolated{T,V,F1,F2}
		channel::T
		cutoff::F1
		cutoffratio::F2
		itr::V
	end
	function interpolated(channel::AbstractxDD, cutoff::Real)
		eth = m2e(sum(channel.ms))
		ev = eth:0.05:(cutoff+0.1)
		calv = ρ_thr.(Ref(channel),ev)
		itr = interpolate((ev,), real.(calv), Gridded(Linear()))
		cutoffratio = real(ρ_thr(channel,cutoff))/ρ_tb(channel,cutoff)
		interpolated(channel,cutoff,cutoffratio,itr)
	end
	function ρ_thr(d::interpolated, e::Real)
	e < d.cutoff ?
		d.itr(e) :
		ρ_tb( d.channel, e) * d.cutoffratio
	end
	ρ_thr(d::interpolated, e::Complex) = ρ_thr(d.channel,e)
end

# ╔═╡ 22c06ab9-8a74-47e7-aa91-268be952a239
ichannels = interpolated.(channels, 4.0)

# ╔═╡ 9d35f73d-f6a6-43a5-86ff-92eb830f62b4
ρInf = sum(ich.cutoffratio for ich in ichannels)

# ╔═╡ 86555b5d-65f9-48ce-84b4-d3e20058efd9
w_matching = ichannels[1].cutoffratio*3/2 / ρInf * 2/e2m(0)

# ╔═╡ ddef2127-6cee-4361-ba5f-dd774fd0d199
function dispersive(d::interpolated,e)
	s = e2m(e)^2
	eth = m2e(sum(d.channel.ms))
	function integrand(s′)
		e′ = m2e(sqrt(s′))
		ρ_thr(d,e′)/s′/(s′-s-1e-6im)
	end
	s/π * quadgk(integrand, e2m(eth)^2, Inf)[1]
end

# ╔═╡ 436ed12b-870a-4667-981f-4e268a047718
denominator_I(d::interpolated,e,δm) = dispersive(d,e)-real(dispersive(d,δm))

# ╔═╡ ce510797-dcbe-407e-ad86-deab4a1fe3ae
let # π⁺D⁰D⁰
	ev = range(-1,3.0,length=100)
	plot()
	f(e) = denominator_I(Tuple(ichannels), e, δm_ΔM.val) / ρInf
	g(e) = denominator_I(NonRelBW(), e+1e-6im, δm_ΔM.val)*w_matching
	# 	real
	plot!(e->imag(f(e)), ev, lab="", l=(:red,))
	plot!(e->imag(g(e)), ev, lab="", l=(:red,:dash))
	# 	imag
	g0 = real(g(δm_ΔM.val))
	plot!(e->real(f(e)), ev, l=(:black,), lab=L"\mathrm{advanced\,\,BW}")
	plot!(e->real(g(e)-g0), ev, l=(:black,:dash), lab=L"\mathrm{non}\textrm{-}\mathrm{rel}")
	#
	vline!([δm_ΔM.val], lab="", c=2)
	vline!([0 m2e(mDˣ⁰+mD⁺)], lab="", c=[:green :magenta], lw=0.5)
	plot!(xlab=L"\delta' m\,\,[\mathrm{MeV}]",
		  ylab=L"\mathcal{A}^{-1}")
end

# ╔═╡ c0ae7968-e297-4126-8534-68f579cce09a
denominator_I(Tuple(ichannels), 0.0, δm_ΔM) / ρInf / w_matching

# ╔═╡ 7c73adc8-f28e-455d-9ac6-1470796eb665
0.196/real(denominator_I(Tuple(ichannels), 0.0, δm_ΔM) / ρInf / w_matching)

# ╔═╡ b374cd78-1b58-4642-acff-77976d30b395
begin
	plot(e->real(dispersive(ichannels[1],e)), -3, 50,
		xlab=L"\delta'm_0\,\,(\mathrm{MeV})", ylab=L"\Re\,\xi(\delta'm_0)", lab="")
	hline!([34.0], lab="pole")
end

# ╔═╡ 12ccd28b-8d81-4be9-ba99-f81659b6ed30
denominator_II(d::interpolated,e,δm) = imag(e) ≥ 0 ?
	denominator_I(d,e,δm) :
	denominator_I(d,e,δm)+2im*ρ_thr(d,e)

# ╔═╡ 37f832bc-2159-46c1-94ba-097a5babd921
function pole_position(d,δm; init = [δm,-0.1*δm])
	abs2D(x) = abs2(denominator_II(d, complex(x...),δm))
	fitres = optimize(abs2D, init)
	return NamedTuple{(:m_pole, :half_Γ_pole, :invD)}(
		[fitres.minimizer..., fitres.minimum])
end

# ╔═╡ 512681c0-0ec7-49f2-bcf0-995c3f8bfdab
sampledpp = [pole_position(Tuple(ichannels),δm) for δm in δmv]

# ╔═╡ bbb68d19-7b11-4b8a-9e44-a855156bc43a
itr_m, itr_Γ =
	interpolate((δmv,), getproperty.(sampledpp, :m_pole), Gridded(Linear())),
	interpolate((δmv,), 2 .* getproperty.(sampledpp, :half_Γ_pole), Gridded(Linear()))

# ╔═╡ 07b12aa4-e34b-49d9-8587-c05201230ab2
pole_sv = NamedTuple{(:m_pole, :Γ_pole)}([itr_m(δm_sw), itr_Γ(δm_sw)])

# ╔═╡ 55700e24-aa2f-4117-9da8-1cc4c64866b3
itr_m(δm_ΔM), itr_Γ(δm_ΔM)

# ╔═╡ c9ef36b2-72e3-4c0a-9305-f51c441127cf
intensity(e; δm = -0.37) =
	1/abs2(denominator_II(Tuple(ichannels),e,δm))*ρ_thr(ichannels[1],e)

# ╔═╡ 6d6d1a5d-b158-4794-9218-960924aaae33
begin
	xv = collect(range(-1,5,length=100))
	plot(xv, intensity.(xv; δm=δm0), lab="", lw=2, frame=:origin)
end

# ╔═╡ 72a4ef36-af33-4278-886d-63d9f8fe41ab
begin
	plot( e->real(ρ_thr(channels[1],e)), range(-1,4.0,length=40))
	plot!(e->real(ρ_thr(channels[2],e)), range(-1,4.0,length=40))
	plot!(e->real(ρ_thr(channels[3],e)), range(-1,4.0,length=40))
	# vline!([ich2.cutoff])
	# plot!(e->real(ρ_thr(ich2,e)), -1:0.01:1)
end

# ╔═╡ 54328092-fac3-42a0-a861-e322c1077e17
ρ_inf(cutoff) = sum(real(ρ_thr(ch,cutoff))/ρ_tb(ch,cutoff) for ch in channels)

# ╔═╡ 01e7a6d2-2474-4d25-83f0-1d4927da11be
ρ_inf(m2e(3.9))

# ╔═╡ d862aac5-a8ee-464f-a51f-ca479111d78c
begin
	plot( e->real(ρ_thr(ichannels[1],e)), range(-1,10.0,length=100), lab="π⁺D⁰D⁰")
	plot!(e->real(ρ_thr(ichannels[2],e)), range(-1,10.0,length=100), lab="π⁰D⁺D⁰")
	plot!(e->real(ρ_thr(ichannels[3],e)), range(-1,10.0,length=100), lab="γD⁺D⁰")
	vline!(getproperty.(ichannels,:cutoff), lab="cutoff")
	plot!()
end

# ╔═╡ 96065ff9-b6e0-4beb-9bf9-0d81de0fdc57
itr(ρ,xv) = interpolate((xv,), real.(ρ.(xv)), Gridded(Linear()))

# ╔═╡ 880fa92f-4c7f-4f9a-9de1-49683efacf5d
md"""
## Delendences and Tests
"""

# ╔═╡ 9bd7390e-a21c-44c1-bb24-216a0d50650d
import Plots.PlotMeasures.mm

# ╔═╡ d23a416b-7edc-4bc7-828b-a0fb6425dca8

begin
	plot(layout=grid(2,1,heights=(0.30,0.70)), size=(600,500), link=:x)
	#
	plot!(sp=1, x->log(intensity(x; δm=δm_sw.val)), -1, 1.8, lab="", lw=2,
		ticks=nothing, bottom_margin=-5mm, xaxis=false,
		ylab=L"\left|\!\!\mathcal{A}^{U}(s)\,\right|^2 \varrho(s)")
	# 	
	plot!(sp=2, getproperty.(sampledpp, :m_pole),
		 getproperty.(sampledpp, :half_Γ_pole), l=(1,:gray), lab="")
	massv = range(δm_sw.val-δm_sw.err, δm_sw.val+δm_sw.err, length=10)
	plot!(sp=2, itr_m.(massv), itr_Γ.(massv)/2, l=(2,:black), lab="")
	scatter!(sp=2, [itr_m(δm_sw.val)], [itr_Γ(δm_sw.val)/2], m=(5,:black), lab="")
	for bp in vcat(collect.(branch_points.(channels))...)
		plot!(sp=2, [bp, complex(1.8,imag(bp))], lc=:red, lab="")
		scatter!(sp=2, [bp], m=(:o,3,:red), lab="")
	end
	vline!(sp=2, [δm_sw.val], lab="", lc=2, frame=:origin)
	vline!(sp=1, [δm_sw.val], lab="", lc=2, frame=:origin)
# 	
	plot!(sp=2, xlab=L"\Re\,(\delta' m)\,\,[\mathrm{MeV}]",
		  ylab=L"\Im\,(\delta' m)\,\,[\mathrm{MeV}]")
end

# ╔═╡ 56f2ef6b-fcab-433b-b6e6-fbe8444ed0b6
@testset "scalar production" begin
	let
		s,s12,s13,msq = 1,2,3,(2,3,5)
		v = (;s,s12,s13,msq)
		@test p12_p3(v)+v.s12 == p_p12(v)
		@test p13_p2(v)+v.s13 == p_p13(v)
		@test p12_p2(v) + p2_p3(v) == p_p2(v)
		@test p13_p2(v) + msq[2] == p_p2(v)
	end
end

# ╔═╡ b5ecdd32-3d72-4666-8e81-5ad917ac0b98
@testset "ABCG" begin
	let
		s,s12,s13,msq = 1,2,3,(2,3,3)
		v = (;s,s12,s13,msq)
		@test A(v) != 0
		@test B(v) != 0
		@test C(v) != 0
		@test C((;s,s12,s13,msq)) ≈ C((;s12=s13,s13=s12,s,msq))
		@test D(v) != 0
		@test E(v) != 0
		@test G(v) != 0
	end
end

# ╔═╡ Cell order:
# ╟─b3a386e0-ed2e-4ed0-aa91-725fbad987f4
# ╟─6995af97-a29f-4931-ab94-85650d0e48aa
# ╠═9052fdd5-92a2-4ca4-bcea-726e9a3df19a
# ╟─976f5e93-5c0c-4ba4-9e87-b53e0d87ebee
# ╠═14c58917-3204-43fd-bc62-a54f5481d127
# ╠═27ffd025-94e5-40aa-9910-7d766a724dd8
# ╠═f47f2687-3c7a-4127-9675-97be5c22f8c6
# ╟─f825bbfc-a0a4-47e4-9b4e-5da5bc7fe724
# ╠═8a551387-0acb-4aa4-9825-4f46de024492
# ╟─190be781-8951-4ce8-b06f-6e621c609527
# ╠═82df0901-7a5c-496b-a090-c18bce6ba755
# ╠═31268d53-342f-4a29-9219-c13ea17ba418
# ╠═035f6347-7dc7-45e4-83a0-116359128ef4
# ╟─b6284e99-49f4-4513-b478-ee4bdad3ae9c
# ╠═d26c2516-1f71-42e1-b3ff-3004d1d144ab
# ╠═0355ece1-5558-4213-961f-74e115342ea0
# ╠═10eccd1f-3956-4cd8-a416-4da8bf384c12
# ╠═581a83d6-62ea-4d7d-8728-4fddfb533409
# ╟─d808ac86-4893-4927-a120-68cdcd923573
# ╠═30748ddd-9757-461f-b95e-c409e1f81718
# ╠═2b852b33-01da-4724-8294-9ebbcefa03b4
# ╠═49058565-1ca9-409c-b3f1-5a5739610d4e
# ╠═c88d9ff1-699a-40ab-b96e-98732de40754
# ╟─6ae3dc60-696f-489c-ba3d-0969381c3605
# ╠═fec685a1-db09-4913-bb7f-ac38eeda9942
# ╠═d93cc18a-8b23-4a0e-9743-febf5e41f176
# ╟─58176be4-6300-4e5b-9543-73d7af40ebf7
# ╠═ddef2127-6cee-4361-ba5f-dd774fd0d199
# ╟─5cc68d15-0242-4aa0-ab14-69fe07df5366
# ╠═436ed12b-870a-4667-981f-4e268a047718
# ╟─29c0d772-6248-4bb0-9c7e-e3d286a92b4f
# ╠═12ccd28b-8d81-4be9-ba99-f81659b6ed30
# ╟─8c4612df-7d5b-4ac4-a29c-36591225c12e
# ╠═80b3743a-cb41-4edd-9319-9e16530f2e3a
# ╠═c9ef36b2-72e3-4c0a-9305-f51c441127cf
# ╠═37f832bc-2159-46c1-94ba-097a5babd921
# ╟─ba3768a8-17e8-4e9f-b496-a614a028a84e
# ╠═fbd283a8-f291-481d-98cb-941371717b0f
# ╟─72a4ef36-af33-4278-886d-63d9f8fe41ab
# ╠═54328092-fac3-42a0-a861-e322c1077e17
# ╠═01e7a6d2-2474-4d25-83f0-1d4927da11be
# ╠═22c06ab9-8a74-47e7-aa91-268be952a239
# ╠═d862aac5-a8ee-464f-a51f-ca479111d78c
# ╟─5b9360b1-32ac-4827-9046-ea19ab785f44
# ╠═9d35f73d-f6a6-43a5-86ff-92eb830f62b4
# ╠═86555b5d-65f9-48ce-84b4-d3e20058efd9
# ╟─ce510797-dcbe-407e-ad86-deab4a1fe3ae
# ╠═36e62ff2-28d5-471d-9a9e-2c0d531ef540
# ╠═c0ae7968-e297-4126-8534-68f579cce09a
# ╠═7c73adc8-f28e-455d-9ac6-1470796eb665
# ╟─f45a5b27-9c61-4982-a2ef-287282b02d47
# ╠═e634d29d-46eb-464a-abea-d4fd30f228aa
# ╠═512681c0-0ec7-49f2-bcf0-995c3f8bfdab
# ╠═bbb68d19-7b11-4b8a-9e44-a855156bc43a
# ╠═07b12aa4-e34b-49d9-8587-c05201230ab2
# ╠═55700e24-aa2f-4117-9da8-1cc4c64866b3
# ╟─4c4cafe9-d0a8-4daa-a475-287756a6014c
# ╠═b374cd78-1b58-4642-acff-77976d30b395
# ╟─d2a66d89-1965-4b76-bee0-f55dbe787722
# ╠═5a685e9f-fa97-479f-a31a-fa7cea7c82c3
# ╠═6d6d1a5d-b158-4794-9218-960924aaae33
# ╟─d23a416b-7edc-4bc7-828b-a0fb6425dca8
# ╟─68caab0a-81ca-4234-9460-bea389ef2ab4
# ╟─d12ad72e-274d-43c0-96a2-ade860abebb9
# ╟─492f23d1-cd3b-46cd-937e-908c56e19674
# ╠═e7883cde-4609-4385-91f7-6b408d855bb1
# ╟─8fe15b8e-4ffb-4aae-bae6-46a71115caed
# ╟─54f8c1e0-7cec-4341-8210-338d8874b50c
# ╠═8b959ab0-adb7-11eb-179d-77286258f875
# ╠═52589f04-a93b-4b61-a6c7-1fd97b4e8a4f
# ╠═b895e1e3-a65d-479a-bc8b-f794e8a98cdd
# ╠═218df0f7-36ed-4f12-b285-ee675c0dffeb
# ╠═96065ff9-b6e0-4beb-9bf9-0d81de0fdc57
# ╟─880fa92f-4c7f-4f9a-9de1-49683efacf5d
# ╠═d745eee9-7229-4767-8fdc-ea1ecda3e61a
# ╠═1919aca2-8900-4736-8119-bc1dd7b8a7e4
# ╠═9bd7390e-a21c-44c1-bb24-216a0d50650d
# ╠═56f2ef6b-fcab-433b-b6e6-fbe8444ed0b6
# ╠═b5ecdd32-3d72-4666-8e81-5ad917ac0b98
