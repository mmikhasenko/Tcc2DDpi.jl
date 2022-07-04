using Pkg
cd(joinpath(@__DIR__,".."))
Pkg.activate(".")
# 
using X2DDpi
using Parameters
using Measurements
using Interpolations
using Statistics
using Optim

using Plots
import Plots.PlotMeasures: mm
using LaTeXStrings
theme(:wong2, frame=:box, grid=false, minorticks=true, 
    guidefontvalign=:top, guidefonthalign=:right,
    xlim=(:auto,:auto), ylim=(:auto,:auto),
    xlab=L"\delta m_{D^0D^0\pi^+}\,\,[\mathrm{MeV}]", lw=1.2)


settings = transformdictrecursively!(readjson("settings.json"), ifstringgivemeasurement)
#
@unpack δm0 = settings["fitresults"]
const δm0_val = δm0.val

# construct the model
sqrtpih(z) = sqrt(z*cis(-π/2))*cis(+π/4)
sqrt_λ_fact(x,y,z) = sqrtpih(x-(sqrt(y)+sqrt(z))^2) * sqrtpih(x-(sqrt(y)-sqrt(z))^2)
k1(s) = sqrt_λ_fact(s,(mDˣ⁺^2-1im*mDˣ⁺*ΓDˣ⁺),mD⁰^2)/s
k2(s) = sqrt_λ_fact(s,(mDˣ⁰^2-1im*mDˣ⁰*ΓDˣ⁰),mD⁺^2)/s
#
D(s) =  - 1im*k1(s) - 1im*k2(s) +
    real(1im*k1(e2m(δm0_val)^2) + 1im*k2(e2m(δm0_val)^2))
#
# 
@assert real(D(e2m(δm0_val)^2)) == 0.0

# plot spectrrum
plot(e->abs(1/D(e2m(e)^2+1e-6im)), -0.9,1.5)
# 

# plot complex plane
sur(x,y) = abs2(D(e2m(x+1im*y)^2))
heatmap(-0.5:0.01:0.3, -0.1:0.01:0.1, log ∘ sur)

# seacrch for a pole
pole_position = let
    fr = optimize(x->sur(x[1],x[2]), [δm0_val,-56.0e-3/2], BFGS())
    NamedTuple{(:m_pole, :half_Γ_pole, :invDabs2)}(
            [fr.minimizer..., fr.minimum])
end

writejson(joinpath("results","nominal","pole_Dstarcomplexmass.json"), transformdictrecursively!(
        Dict{Symbol,Any}(
            :pole_position => pole_position
        ), ifmeasurementgivestring)
    )
