
using Plots
using LaTeXStrings
theme(:wong2, frame=:box, grid=false, minorticks=true, 
    guidefontvalign=:top, guidefonthalign=:right)

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
    plot()
    f(e) = denominator_I(Tuple(ichannels), e, δm0.val) / ρInf / (w_matching*1e3)
    g(e) = denominator_I(NonRelBW(), e+1e-6im, δm0.val)
    # 	real
    plot!(mgrid, invA_advans_real, lab="", l=(:red,))
    plot!(mgrid, invA_nonrel_real, lab="", l=(:red,:dash))
    # 	imag
    plot!(mgrid, invA_advans_imag, l=(:black,), lab=L"\mathrm{advanced\,\,BW}")
    plot!(mgrid, invA_nonrel_imag, l=(:black,:dash), lab=L"\mathrm{non}\textrm{-}\mathrm{rel}")
    #
    vline!([δm0.val], lab="", c=2)
    vline!([0 m2e(mDˣ⁰+mD⁺)], lab="", c=[:green :magenta], lw=0.5)
    plot!(xlab=L"\delta' m\,\,[\mathrm{MeV}]",
        ylab=L"\mathcal{A}^{-1}", leg=:topright)
end
savefig(joinpath("plots","nominal","inverse_amplitude.pdf"))

# latex
pgfplotsx()
savefig(joinpath("plots","latex","inverse_amplitude.tex"))
