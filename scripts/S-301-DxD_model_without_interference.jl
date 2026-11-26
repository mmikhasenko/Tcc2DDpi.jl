using Pkg
cd(joinpath(@__DIR__,".."))
Pkg.activate(".")
# 
using X2DDpi
using Parameters
using QuadGK

settings = transformdictrecursively!(readjson("settings.json"), ifstringgivemeasurement)
@unpack e_cutoff, e_step = settings["phspmatching"]
m_cutoff = e2m(e_cutoff)
@unpack δm0 = settings["fitresults"]
const δm0_val = δm0.val
#

channel = DˣD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), (m=mDˣ⁺,Γ=ΓDˣ⁺))

ichannel = interpolated(channel, m_cutoff; mstep = e_step / 1e3) # cutoff

const ampX = Amplitude(ichannel, zero)

f(x,y) = y<0 ? 
    denominator_II(ampX, x+1im*y, δm0_val) :
    denominator_I(ampX, x+1im*y, δm0_val)
# 
ff = log ∘ abs2 ∘ f
xv = -0.5:0.003:0.10
yv = -0.2:0.003:0.14
@time calv = -ff.(xv', yv)

# contour(xv, yv, calv, levels=200)

writejson(joinpath("results","nominal","complexplane_model_without_interference.json"), Dict(
    :x => xv,
    :y => yv,
    :z => zv,
))
