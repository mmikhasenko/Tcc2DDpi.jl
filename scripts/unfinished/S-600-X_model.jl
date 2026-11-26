using Pkg
Pkg.activate(joinpath(@__DIR__, "..", ".."))
cd(joinpath(@__DIR__, "..", ".."))
# 
using X2DDpi
using Parameters

settings = transformdictrecursively!(readjson("settings.json"), ifstringgivemeasurement)
@unpack e_cutoff, e_step = settings["phspmatching"]
m_cutoff = e2m(e_cutoff)
#
channels = [
    πDD((m1 = mπ⁰, m2 = mD⁰, m3 = mD⁰), BW(m = mDˣ⁰, Γ = ΓDˣ⁰), BW(m = mDˣ⁰, Γ = ΓDˣ⁰)),
    πDD((m1 = mγ, m2 = mD⁰, m3 = mD⁰), BW(m = mDˣ⁰, Γ = ΓDˣ⁰), BW(m = mDˣ⁰, Γ = ΓDˣ⁰))]
#
ichannels = interpolated.(channels, m_cutoff; mstep = e_step / 1e3) # cutoff

# save
writejson(joinpath("results", "nominal", "X_model.json"), Dict(
    :channels => X2DDpi.obj2nt.(channels),
    :ichannels => X2DDpi.obj2nt.(ichannels),
)
)
