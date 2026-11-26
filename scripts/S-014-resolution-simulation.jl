using Pkg
cd(joinpath(@__DIR__, ".."))
Pkg.activate(".")

using Plots
using LaTeXStrings
theme(:wong2, size = (500, 350), minorticks = true, grid = false, frame = :box,
    guidefontvalign = :top, guidefonthalign = :right,
    foreground_color_legend = nothing,
    legendfontsize = 9, legend = :topright,
    xlim = (:auto, :auto), ylim = (:auto, :auto),
    xlab = L"\delta' m\,\,(\mathrm{MeV})",
    ylab = "Cross section",
    fontfamily = "Computer Modern")

using Random
using FHist
using X2DDpi
using Parameters
using NumericalDistributions: NumericallyIntegrable, fft_convolve
using NumericalDistributions.Distributions
using NumericalDistributions: NumericalDistributions

settings = transformdictrecursively!(
    readjson("settings.json"),
    ifstringgivemeasurement)
@unpack δm0 = settings["fitresults"]

#
# retrieve the model
modelDict = readjson(joinpath("results", "nominal", "model.json"))
ichannels = interpolated.(d2nt.(modelDict["ichannels"]))
const ampX0 = Amplitude(Tuple(ichannels), zero) # zero is 1/Γ
const ich_D⁰D⁰π⁺ = ichannels[1]

@time intensity_dist = let
    NumericalDistributions.interpolated(range(-2.0, 4.0, 200)) do e
        abs2(1 / denominator_I(ampX0, e, δm0.val)) * ρ_thr(ich_D⁰D⁰π⁺, e2m(e))
    end
end

let
    plot()
    plot!(x -> pdf(intensity_dist, x), -1, 2, l = (:black,), lab = "LHCb model")
    #
    vline!([0], lab = "", c = :green, lw = 1, yscale = :log10)
end


let
    σR = 0.3
    resolution_dist = truncated(Normal(0, σR), -5σR, 5σR)
    cov_dist = fft_convolve(intensity_dist, resolution_dist)

    Random.seed!(1223)
    nEvents = 300
    data = rand(cov_dist, nEvents)
    h = Hist1D(data, binedges = range(-2, 4, 21))
    plot(ylim = (0, 0.8 * nEvents))
    scatter!(bincenters(h), bincounts(h), yerr = sqrt.(bincounts(h)), lab = "Simulations")
    plot!(h, seriestype = :stepbins, fill = 0, alpha = 0.1, lab = "")
    vline!([0.0], lab = "", c = :purple)
end
