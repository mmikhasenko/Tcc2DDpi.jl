using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))
cd(joinpath(@__DIR__, "..", ".."))
# 
using X2DDpi
using Parameters

settings = transformdictrecursively!(readjson("settings.json"), ifstringgivemeasurement)
@unpack cutoff, estep = settings["phspmatching"]
#
channels = [
    πDD((m1 = mπ⁰, m2 = mD⁰, m3 = mD⁰), BW(m = mDˣ⁰, Γ = ΓDˣ⁰), BW(m = mDˣ⁰, Γ = ΓDˣ⁰)),
    πDD((m1 = mγ, m2 = mD⁰, m3 = mD⁰), BW(m = mDˣ⁰, Γ = ΓDˣ⁰), BW(m = mDˣ⁰, Γ = ΓDˣ⁰))]
#
ichannels = interpolated.(channels, cutoff; estep = estep) # cutoff

# save
writejson(joinpath("results", "nominal", "X_model.json"), Dict(
    :channels => X2DDpi.obj2nt.(channels),
    :ichannels => X2DDpi.obj2nt.(ichannels),
)
)
