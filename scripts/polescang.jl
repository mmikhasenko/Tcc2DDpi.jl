using Pkg
cd(joinpath(@__DIR__,".."))
Pkg.activate(".")
# 
using X2DDpi
using Parameters
using Measurements
using CSV
using DataFrames

settings = transformdictrecursively!(readjson("settings.json"), ifstringgivemeasurement)
#
@unpack δm0 = settings["fitresults"]

# retrieve the model
modelDict = readjson(joinpath("results","nominal","model.json"))
ichannels = interpolated.(d2nt.(modelDict["ichannels"]))
channels = getproperty.(ichannels, :channel)
# 
import X2DDpi:pole_position
pole_position(δm,g) = pole_position(
    Amplitude(Tuple(ichannels),
        e->(e2m(δm)^2-e2m(e)^2)/(g^2/2)),
    δm)
    #
gdm = CSV.File(joinpath("results", "scang.txt"))
nt_gdm = collect(NamedTuple.(gdm))

nt_gdm = [(g=1e4, dm = δm0.val, δdm = δm0.err), nt_gdm...]
nt_gdm[end] = (; nt_gdm[end]..., δdm = 0.043)
# push!(nt_gdm) = (; nt_gdm[end]..., g=0.01, δdm = 0.043)

pole_mddl = [pole_position( row[2]        /1000, row[1]) for row in gdm]
pole_left = [pole_position((row[2]-row[3])/1000, row[1]) for row in gdm]
pole_rght = [pole_position((row[2]+row[3])/1000, row[1]) for row in gdm]

writejson(joinpath("results","nominal","scang.json"),
        Dict(
            :middle => pole_mddl,
            :middle_plus_error => pole_left,
            :middle_minus_error => pole_rght,
        )
    )
#