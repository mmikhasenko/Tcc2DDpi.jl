using Pkg
cd(joinpath(@__DIR__,".."))
Pkg.activate(".")

using Plots
using LaTeXStrings
theme(:wong2, size=(500,350), minorticks=true, grid=false, frame=:box,
    guidefontvalign=:top, guidefonthalign=:right,
    foreground_color_legend = nothing,
    legendfontsize=9, legend =:topright,
    xlim=(:auto,:auto), ylim=(:auto,:auto))

using X2DDpi
using Parameters

settings = transformdictrecursively!(
    readjson("settings.json"),
    ifstringgivemeasurement)
@unpack δm0 = settings["fitresults"]

@unpack mgrid, invA_nonrel_real, invA_nonrel_imag, invA_advans_real, invA_advans_imag =
    readjson(joinpath("results","nominal","inverse_amplitude.json"))["inverse_amplitude"]

let
    ev = range(-1,3.0,length=100)
    plot(xlab=L"\delta' m\,\,(\mathrm{MeV})", ylab=L"\mathrm{Re}\,\mathcal{A}^{-1}(black),\,\,\mathrm{Im}\,\mathcal{A}^{-1}(red)")
    f(e) = denominator_I(Tuple(ichannels), e, δm0.val) / (w_matching*ρInf)
    g(e) = denominator_I(NonRelBW(), e+1e-6im, δm0.val)
    # 	real
    plot!(mgrid, invA_advans_real, l=(:black,), lab="LHCb model")
    plot!(mgrid, invA_nonrel_real, l=(:black,:dash), lab="Non-rel. eff.-range exp.")
    # 	imag
    plot!(mgrid, invA_advans_imag, l=(:red,), lab="")
    plot!(mgrid, invA_nonrel_imag, l=(:red,:dash), lab="")
    #
    # vline!([δm0.val], lab="", c=2)
    vline!([0], lab="", c=:green, lw=1)
end
savefig(joinpath("plots","nominal","effectiverangenonrel.pdf"))

