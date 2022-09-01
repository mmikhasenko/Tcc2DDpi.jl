### A Pluto.jl notebook ###
# v0.19.5

using Markdown
using InteractiveUtils

# ╔═╡ c26b8580-e57a-11ec-37cd-f9e7a240c1de
begin
	cd(joinpath(@__DIR__, ".."))
	import Pkg
	Pkg.activate(".")
	Pkg.instantiate()
	# 
	using X2DDpi
	using Plots
	using LaTeXStrings
	using Parameters
	using Measurements
end

# ╔═╡ fb0a593c-35b7-4248-853a-4c82fc5be7c6
theme(:wong2, size=(500,350), minorticks=true, grid=false, frame=:box,
    guidefontvalign=:top, guidefonthalign=:right,
    foreground_color_legend = nothing,
    legendfontsize=9, legend =:topright,
    xlim=(:auto,:auto), ylim=(:auto,:auto))

# ╔═╡ 483d68b3-de26-417f-a1af-b0b19d6b3cdd
settings = transformdictrecursively!(
	readjson("settings.json"), ifstringgivemeasurement)

# ╔═╡ 4ce11392-a99e-40ff-84e3-8b045a8387eb
begin
	@unpack δm0 = settings["fitresults"]
	@unpack Γ0_90CL_syst, Γ0_95CL_syst = settings["fitresults"]
	@unpack Γ0_90CL_stat, Γ0_95CL_stat = settings["fitresults"]
end

# ╔═╡ 2b22abef-3d44-4650-9365-275a1a855b45
const δm0_val = δm0.val

# ╔═╡ a0aea720-d9b9-4b80-b176-ce1837b04c3a
begin
	modelDict = readjson(joinpath("results","nominal","model.json"))
	const ichannels = interpolated.(d2nt.(modelDict["ichannels"]))
	channels = getproperty.(ichannels, :channel)
end

# ╔═╡ 12f1d2f3-5003-4a44-9dee-d7222ea038bc
const ampX0 = Amplitude(Tuple(ichannels), zero) # zero is 1/Γ

# ╔═╡ 2ce57cb3-069d-4adb-a212-d9b5216585a5
D(e) = denominator_I(ampX0, e, δm0_val)

# ╔═╡ c31e7105-4fe1-4303-80bd-d2067f0982b9
I(e) = 1/abs2(D(e))*ρ_thr(ampX0.ik[1], e)

# ╔═╡ 0fe42f23-ece0-4256-aaf2-6825b27633f8
let
	plot(I, -1, 2, yscale=:log10, lab="")
	plot!(xlab=L"p_\mathrm{rel}\,\,(\mathrm{MeV})",
		ylab=L"\mathrm{decay\,\,rate}\,\,(\mathrm{a.u.})")
	vline!([0 1e3*(mDˣ⁰+mD⁺-mDˣ⁺-mD⁰)],
		lab=[L"D^{*+}D^0" L"D^{*0}D^+"])
end

# ╔═╡ 90f76e7d-ec62-43cb-bcfd-1e48a340bf52
function p2e(p; m1=mDˣ⁺, m2=mD⁰) # D*+D0
	_p = p<0 ? 1im*p/1e3 : p/1e3
	E(p,m) = sqrt(p^2+m^2)
	m² = m1^2 + m2^2 + 2*E(_p,m1)*E(_p, m2) + 2*(_p)^2
	return (sqrt(m²) - m1 - m2)*1e3 |> real
end

# ╔═╡ 828074ac-ce86-4f30-935a-77be82bdabdc
e2p(e; m1=mDˣ⁺, m2=mD⁰) =
	sqrt(X2DDpi.λ((e/1e3+m1+m2)^2,m1^2,m2^2))/2/(e/1e3+m1+m2)*1e3

# ╔═╡ c8044366-8dfe-45be-8739-dc372ef827d7
begin
	plot(p->I(p2e(p)), 0, 90, lab="", lw=2, size=(500,380))
	plot!(p->I(p2e(p)), -50, 0, lab="", lw=2, lc=:lightgray)
	# plot!(yscale=:log10)
	plot!(ylims=(0,1.5e3))
	plot!(xlab=L"p_\mathrm{rel}\,\,(\mathrm{MeV})",
		ylab=L"\mathrm{decay\,\,rate}\,\,(\mathrm{a.u.})")
	vline!(e2p.([0im 1e3*(mDˣ⁰+mD⁺-mDˣ⁺-mD⁰)]) .|> real,
		lab=[L"D^{*+}D^0" L"D^{*0}D^+"])
end

# ╔═╡ a417f197-7655-483b-8c5c-ca286703e72e
savefig(joinpath("plots", "coorelation_function.pdf"))

# ╔═╡ Cell order:
# ╠═c26b8580-e57a-11ec-37cd-f9e7a240c1de
# ╠═fb0a593c-35b7-4248-853a-4c82fc5be7c6
# ╠═483d68b3-de26-417f-a1af-b0b19d6b3cdd
# ╠═4ce11392-a99e-40ff-84e3-8b045a8387eb
# ╠═2b22abef-3d44-4650-9365-275a1a855b45
# ╠═a0aea720-d9b9-4b80-b176-ce1837b04c3a
# ╠═12f1d2f3-5003-4a44-9dee-d7222ea038bc
# ╠═2ce57cb3-069d-4adb-a212-d9b5216585a5
# ╠═c31e7105-4fe1-4303-80bd-d2067f0982b9
# ╠═0fe42f23-ece0-4256-aaf2-6825b27633f8
# ╠═90f76e7d-ec62-43cb-bcfd-1e48a340bf52
# ╠═828074ac-ce86-4f30-935a-77be82bdabdc
# ╠═c8044366-8dfe-45be-8739-dc372ef827d7
# ╠═a417f197-7655-483b-8c5c-ca286703e72e
