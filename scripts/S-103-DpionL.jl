using Pkg
cd(joinpath(@__DIR__,".."))
Pkg.activate(".")

using X2DDpi
using AlgebraPDF
# 
using Plots
using LaTeXStrings
theme(:wong2, size=(500,350), minorticks=true, grid=false, frame=:box,
    guidefontvalign=:top, guidefonthalign=:right,
    foreground_color_legend = nothing,
    legendfontsize=9, legend =:topright,
    xlim=(:auto,:auto), ylim=(:auto,:auto),
    xlab=L"m(D^0\pi^+)\,[\mathrm{GeV}]")

# 
modelDict = readjson(joinpath("results","nominal","model.json"))
ichannels = interpolated.(d2nt.(modelDict["ichannels"]))
channels = getproperty.(ichannels, :channel)

const ms1 = channels[1].ms
const nt3m = NamedTuple{(:m1,:m2,:m3)}
const parsDˣ = channels[1].R12
# 
const λ = X2DDpi.λ
momp(m,M, ms::nt3m) = ms.m1+ms.m2 < m ? sqrt(λ(m^2,ms.m1^2,ms.m2^2))/m : 0.0
momq(m,M, ms::nt3m) = m < M-ms.m3 ? sqrt(λ(m^2,M^2,ms.m3^2))/M : 0.0

#
rho_pqL(m,M, ms::nt3m; p) = momp(m,M, ms)^(2p.Lp+1)*momq(m,M, ms)^(2p.Lq+1)

d0 = FunctionWithParameters((m;p)->rho_pqL(m,p.M, ms1; p); p=(Lp=1,Lq=0,M=e2m(-0.3)))
nd0 = Normalized(d0, (ms1.m1+ms1.m2, e2m(0)-ms1.m3)) 
dD = FunctionWithParameters((x;p)->J_I(x^2,parsDˣ), ∅)
nd1 = Normalized(d0*abs2(dD), (ms1.m1+ms1.m2, e2m(0)-ms1.m3))

p1 = let
    plot(leg=:topleft, title=L"\mathrm{}\delta'm_0 =-0.3\,\mathrm{MeV}")
    plot!(updatepars(nd1,(Lp=1,Lq=0)),1,300, lab="P-wave + S-wave")
    plot!(updatepars(nd1,(Lp=1,Lq=1)),1,300, lab="P-wave + P-wave")
    plot!(updatepars(nd1,(Lp=1,Lq=2)),1,300, lab="P-wave + D-wave")
end

p2 = let
    plot(leg=:topleft, title=L"\mathrm{with\,\,no\,\,}D^{*+}\mathrm{\,\,propagator}")
    plot!(updatepars(nd0,(Lp=1,Lq=0)), lab="P-wave + S-wave")
    plot!(updatepars(nd0,(Lp=1,Lq=1)), lab="P-wave + P-wave")
    plot!(updatepars(nd0,(Lp=1,Lq=2)), lab="P-wave + D-wave")
end

plot(p1,p2, layout=grid(2,1), size=(400,700))
savefig(joinpath("plots","nominal","Dpi_spectrum.pdf"))


