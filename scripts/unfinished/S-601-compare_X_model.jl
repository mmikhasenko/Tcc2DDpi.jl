using Pkg
cd(joinpath(@__DIR__, "..", ".."))
Pkg.activate(".")
# 
using X2DDpi
using Parameters
using Measurements
using Interpolations
using Statistics
using Plots

settings = transformdictrecursively!(readjson("settings.json"), ifstringgivemeasurement)
#
@unpack δm0 = settings["fitresults"]

# retrieve the model
modelDict = readjson(joinpath("results", "nominal", "X_model.json"))
ichannels = interpolated.(d2nt.(modelDict["ichannels"]))
channels = getproperty.(ichannels, :channel)

ampX0 = Amplitude(Tuple(ichannels), zero) # zero is 1/Γ

let
    plot(e -> ρ_thr(ichannels[1], e) / ichannels[1].cutoffratio, -5.0, 0.0, lab = "π⁰")
    plot!(e -> ρ_thr(ichannels[2], e) / ichannels[2].cutoffratio, -5.0, 0.0, lab = "γ")
end

let
    plot(e -> ρ_thr(ichannels[1].channel, e) / ichannels[1].cutoffratio, -5.0, 0.0, lab = "π⁰")
    plot!(e -> ρ_thr(ichannels[2].channel, e) / ichannels[2].cutoffratio, -5.0, 0.0, lab = "γ")
end
