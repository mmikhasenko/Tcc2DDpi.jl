module xDDPhaseSpace

using Cuba
using Parameters

# Reexport all methods from submodules
export AbstractxDD
export mapdalitzmethod, masses
export integrand_mapped_thr, ρ_thr
include("abstractxdd.jl")


export e2m
export WithThrE
include("withthr.jl")


import Base.^
^(ms::NamedTuple{(:m1, :m2, :m3), T} where {T}, n::Int) = Tuple(ms) .^ n

export σ2of3_pm, σ3of1_pm, σ3of2_pm
export σ3of1, σ2of1
export AbstractDalitzMapping
export LinearDalitzMapping, HookSqrtDalitzMapping
export mapdalitz
include("integrationpaths.jl")

export BW, BW_norm, BW_Swave, ZeroBW
export Jᴵ, Jᴵᴵ
include("isobar.jl")

export πDD
export branch_points
export constructchannel
include("covariant_exptessions.jl")
include("mainmodel-pi.jl")

export obj2nt
include("serialization.jl")



end # xDDPhaseSpace
