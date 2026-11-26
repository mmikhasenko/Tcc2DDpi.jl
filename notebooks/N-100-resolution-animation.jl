### A Pluto.jl notebook ###
# v0.20.8

using Markdown
using InteractiveUtils

# ╔═╡ 61b01d0b-1f8f-4e9d-805a-b533cd39410c
begin
    import Pkg
    Pkg.activate(@__DIR__)
end

# ╔═╡ 96463f8e-8c97-11f0-3099-2dfd09730179
# ╠═╡ show_logs = false
begin
	cd(joinpath(@__DIR__, ".."))
	using Plots
	using LaTeXStrings
	using Random
	using FHist
	using X2DDpi
	using Parameters
	using NumericalDistributions: NumericallyIntegrable, fft_convolve
	using NumericalDistributions.Distributions
	import NumericalDistributions
end

# ╔═╡ 5f534c6c-14d5-4665-a39e-e9a4aa198e38
md"""
# Experimental resolution

The notebook demonstrates how experimental resolution impacts the mass spectrum of particles and sensitivity to resonance physics
"""

# ╔═╡ 0135ec0d-da27-45b5-aefe-a5fcdef0a22c
theme(:wong2, size = (500, 350), dpi=300,
	minorticks = true, grid = false, frame = :box,
	guidefontvalign = :top, guidefonthalign = :right,
    foreground_color_legend = nothing,
	legendfontsize = 9, legend = :topright,
    xlim = (:auto, :auto), ylim = (:auto, :auto),
	xlab = L"\delta m_{D^0D^0\pi^+}\,\,(\mathrm{MeV})",
    ylab = "Cross section",
	fontfamily = "Computer Modern")

# ╔═╡ a60f4ac1-aac7-4165-9643-d3975d508962
begin
	settings = transformdictrecursively!(
	    readjson("settings.json"),
	    ifstringgivemeasurement)
	@unpack δm0 = settings["fitresults"]
	#
	# retrieve the model
	modelDict = readjson(joinpath("results", "nominal", "model.json"))
	ichannels = X2DDpi.interpolated.(d2nt.(modelDict["ichannels"]))
	const ampX0 = Amplitude(Tuple(ichannels)) # zero is 1/Γ
	const ich_D⁰D⁰π⁺ = ichannels[1]
end;

# ╔═╡ 8467b165-92d4-4c0b-b5a1-84ab33674ed6
@time intensity_dist =
    NumericalDistributions.interpolated(range(-2.0, 3.0, 500)) do e
        abs2(1 / denominator_I(ampX0, e, δm0.val)) * ρ_thr(ich_D⁰D⁰π⁺, e2m(e))
    end;

# ╔═╡ b9458d32-ea2a-4ea3-a940-07b257278c7b
let
    plot(size=(300,250))
    plot!(ylims=(0,:auto)) # , yscale=:log10
    plot!(-1, 3, l = (:black, 2), lab = "LHCb model", fill=0, fillalpha=0.3, lc=1, c=2) do x
		pdf(intensity_dist, x)
	end
    #
    vline!([0], lab = "", l = (:green, :dash), lw = 1)
	savefig(joinpath(get(ENV, "TCC_EXPORT_OUTPUTS", joinpath(@__DIR__, "..", "plots")), "no-resolution.png"))
	plot!()
end

# ╔═╡ 30228da0-0da7-497c-be39-f56817528259
function gaus_conv_dist(σ, model)
	resolution_dist = truncated(Normal(0, σ), -5σ, 5σ)
    cov_dist = fft_convolve(model, resolution_dist)
end;

# ╔═╡ 5e992883-030e-4980-9b8c-be43f39d1f97
function plot_data_styled(data; σR::Number)
	n = length(data)
	h = Hist1D(data, binedges = range(-2, 4, 21))
	#
    plot(ylim = (0, 0.8 * n))
    scatter!(bincenters(h), bincounts(h), yerr = sqrt.(bincounts(h)), lab = "Simulation")
    plot!(h, seriestype = :stepbins, fill = 0, alpha = 0.3, lab = "")
    vline!([0.0], lab = "", l = (:green, :dash))
	annotate!([
        (0, 0.3n, text("Dˣ⁺D⁰",10, :left, :bottom, rotation=90)),
		(2, 0.2n, text("\$\\sigma\$ = $(round(σR; digits=1)) MeV", 15, :bottom))
    ])
end;

# ╔═╡ 0b15c750-69f6-494b-9b04-279f4ef2dd85
let
    σR = 0.1
    cov_dist = gaus_conv_dist(σR, intensity_dist)
	#
    Random.seed!(1223)
    nEvents = 300
    data = rand(cov_dist, nEvents)
	plot_data_styled(data; σR)
	plot!(title="Test plot")
end

# ╔═╡ 93003e04-e906-4f53-a750-52e91731b1cb
scan_σR_values = let
	_v = [0.01, 0.1:0.02:0.3..., 0.3:0.2:3.0...]
	# _v = [0.01, 0.1:0.2:0.3..., 0.5:0.8:3.0...] # debug
	vcat(_v, reverse(_v))
end;

# ╔═╡ 6be6cba8-a5fb-4f62-adb5-6402c3f69c91
# ╠═╡ show_logs = false
let
	anim = Plots.Animation()
	for σR in scan_σR_values
	    cov_dist = gaus_conv_dist(σR, intensity_dist)
		#
	    Random.seed!(1223)
	    nEvents = 300
	    data = rand(cov_dist, nEvents)
		plot_data_styled(data; σR)
		#
		Plots.frame(anim)
	end
	_gif = gif(anim, joinpath(get(ENV, "TCC_EXPORT_OUTPUTS", joinpath(@__DIR__, "..", "plots")), "resolution.gif"); fps=5, loop=0)
end

# ╔═╡ Cell order:
# ╠═61b01d0b-1f8f-4e9d-805a-b533cd39410c
# ╟─5f534c6c-14d5-4665-a39e-e9a4aa198e38
# ╠═96463f8e-8c97-11f0-3099-2dfd09730179
# ╠═0135ec0d-da27-45b5-aefe-a5fcdef0a22c
# ╠═a60f4ac1-aac7-4165-9643-d3975d508962
# ╠═8467b165-92d4-4c0b-b5a1-84ab33674ed6
# ╠═b9458d32-ea2a-4ea3-a940-07b257278c7b
# ╠═30228da0-0da7-497c-be39-f56817528259
# ╠═5e992883-030e-4980-9b8c-be43f39d1f97
# ╠═0b15c750-69f6-494b-9b04-279f4ef2dd85
# ╠═93003e04-e906-4f53-a750-52e91731b1cb
# ╠═6be6cba8-a5fb-4f62-adb5-6402c3f69c91
