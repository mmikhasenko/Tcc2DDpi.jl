using Pkg
cd(joinpath(@__DIR__,".."))
Pkg.activate(".")
# 
using X2DDpi
using Parameters

settings = transformdictrecursively!(readjson("settings.json"), ifstringgivemeasurement)
@unpack cutoff, estep = settings["phspmatching"]
#
channels = [
    πDD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), BW(m=mDˣ⁺,Γ=ΓDˣ⁺), BW(m=mDˣ⁺,Γ=ΓDˣ⁺)),
    πDD((m1=mπ⁰,m2=mD⁺,m3=mD⁰), BW(m=mDˣ⁺,Γ=ΓDˣ⁺), BW(m=mDˣ⁰,Γ=ΓDˣ⁰)),
    γDD((m1=mγ, m2=mD⁺,m3=mD⁰), BW(m=mDˣ⁺,Γ=ΓDˣ⁺), BW(m=mDˣ⁰,Γ=ΓDˣ⁰))]
#
ichannels = interpolated.(channels, cutoff; estep=estep) # cutoff

# save
writejson(joinpath("results","nominal","model.json"), Dict(
        :channels => obj2nt.(channels),
        :ichannels => obj2nt.(ichannels)
        )
    )
#
