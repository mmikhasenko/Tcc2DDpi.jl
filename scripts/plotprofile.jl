using Pkg
Pkg.activate(joinpath(@__DIR__,".."))
# 
using X2DDpi
using Parameters
# 
using Plots
using LaTeXStrings
theme(:wong2, size=(500,350), minorticks=true, grid=false, frame=:box,
    guidefontvalign=:top, guidefonthalign=:right,
    foreground_color_legend = nothing,
    legendfontsize=9, legend =:topright,
    xlim=(:auto,:auto), ylim=(:auto,:auto))


using AlgebraPDF
# 


settings = transformdictrecursively!(readjson("settings.json"), ifstringgivemeasurement)
#
@unpack δm0 = settings["fitresults"]

# retrieve the model
modelDict = readjson(joinpath("results","nominal","model.json"))
ichannels = interpolated.(d2nt.(modelDict["ichannels"]))
channels = getproperty.(ichannels, :channel)

ampX0 = Amplitude(Tuple(ichannels), zero) # zero is 1/Γ
# 
ρInf = sum(ich.cutoffratio for ich in ichannels)
ampX(Γ) = Amplitude(Tuple(ichannels),
    e->(e2m(δm0.val)^2-e2m(e)^2)/(e2m(δm0.val)*Γ/1e3/ρInf))
# 

FI = FunctionWithParameters(
    (e;p)->1e-9*abs2(1/denominator_I(ampX(p.Γ),e,p.δm)), (δm=δm0.val, Γ=300.0))
# 
let
    plot(xlab=L"\delta\prime m\,\,[\mathrm{MeV}]")
    plot!(updatepar(FI, :Γ, 10000.0), -0.5, 0.1, lab=L"\Gamma=\mathrm{Inf}")
    plot!(updatepar(FI, :Γ, 135.0), -0.5, 0.1, lab=L"\Gamma=135\,\,\mathrm{MeV}")
    plot!(updatepar(FI, :Γ,  95.0), -0.5, 0.1,  lab=L"\Gamma=95\,\,\mathrm{MeV}")
end
savefig(joinpath("plots","nominal","profile.pdf"))
