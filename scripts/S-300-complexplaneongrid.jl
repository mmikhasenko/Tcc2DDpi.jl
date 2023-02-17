using Pkg
cd(joinpath(@__DIR__,".."))
Pkg.activate(".")
# 
using X2DDpi
using Parameters
using Measurements
using Interpolations
using Statistics

settings = transformdictrecursively!(readjson("settings.json"), ifstringgivemeasurement)
#
@unpack δm0 = settings["fitresults"]
const δm0_val = δm0.val


# retrieve the model
modelDict = readjson(joinpath("results","nominal","model.json"))
ichannels = interpolated.(d2nt.(modelDict["ichannels"]))
channels = getproperty.(ichannels, :channel)

ampX = Amplitude(Tuple(ichannels), zero)

dA(x,y) = y<0 ? 
    denominator_II(ampX, x+1im*y, δm0_val) :
    denominator_I(ampX, x+1im*y, δm0_val)

ff = log ∘ abs2 ∘ dA
xv = -0.5:0.003:0.10
yv = -0.2:0.003:0.14
@time zv = -ff.(xv', yv)

# contour(xv, yv, calv, levels=200)

writejson(joinpath("results","nominal","logA_2d_grid.json"), Dict(
    :x => xv,
    :y => yv,
    :z => zv,
))
