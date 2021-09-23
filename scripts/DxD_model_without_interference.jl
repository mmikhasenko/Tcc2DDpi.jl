using Pkg
cd(joinpath(@__DIR__,".."))
Pkg.activate(".")
# 
using X2DDpi
using Parameters
using QuadGK

settings = transformdictrecursively!(readjson("settings.json"), ifstringgivemeasurement)
@unpack cutoff, estep = settings["phspmatching"]
@unpack δm0 = settings["fitresults"]
const δm0_val = δm0.val
#

import X2DDpi: AbstractxDD
struct DˣD{T1,T2} <: AbstractxDD
    ms::NamedTuple{(:m1,:m2,:m3),T1}
    R::NamedTuple{T2}
end

import X2DDpi: decay_matrix_element, f²
decay_matrix_element(d::DˣD,s,σ3,σ2) = J_I(σ3,d.R) * J_II(σ3,d.R) * f²/3/4

channel = DˣD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), (m=mDˣ⁺,Γ=ΓDˣ⁺))

import X2DDpi: integrand_mapped_thr
function integrand_mapped_thr(d::DˣD,s,x)
	# 	
	σ3_0, σ3_e = (d.ms[1]+d.ms[2])^2, (√s-d.ms[3])^2
    # 
    # cut straight
    σ3 = σ3_0 + x[1]*(σ3_e-σ3_0) # straight path
	jac = (σ3_e-σ3_0)
    # 
    # cut down
    # σ3_i = real((√s-d.ms[3])^2)+sign(imag(s))*1e-6im
    # σ3 = x[1] < 0.5 ? 
    #         σ3_0 + ( x[1]     / 0.5)*(σ3_i-σ3_0) : # straight path
	#         σ3_i + ((x[1]-0.5)/ 0.5)*(σ3_e-σ3_i) # straight path
	# #
	# jac = x[1] < 0.5 ?
    #         (σ3_e-σ3_0) / 0.5 :
    #         (σ3_e-σ3_i) / 0.5
    # 
    σ2 = 0.0
	decay_matrix_element(d,s,σ3,σ2) / (2π*s) * jac
end

import X2DDpi: ρ_thr
function ρ_thr(d::DˣD, e::Complex)
	integrand(x) = integrand_mapped_thr(d,e2m(e)^2,[x])
	v = quadgk(integrand, 0, 1)[1]
	complex(v...) / (8π)^2
end

function ρ_thr(d::DˣD, e::Real)
	integrand(x) = real(integrand_mapped_thr(d,e2m(e)^2,[x]))
	v = quadgk(integrand, 0, 1)[1]
	v[1] / (8π)^2
end


import X2DDpi: ρ_tb, λ
function ρ_tb(d::DˣD, e::Real)
	M,m = d.R.m, d.ms.m3
	sqrts = e2m(e)
	sqrts < M+m ? 0.0 :
    	sqrt(λ(e2m(e)^2,M^2,m^2))/e2m(e)^2
end	

ichannel = interpolated(channel, cutoff; estep=estep) # cutoff

const ampX = Amplitude(ichannel, zero)

# using Plots
# theme(:wong2, grid=false, frame=:box)

# using Plots
# plot(ρ_thr.(Ref(ichannel), -1:0.01:20))

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
