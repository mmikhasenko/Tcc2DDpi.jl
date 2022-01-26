module X2DDpi

using QuadGK
using Cuba
using Interpolations
using Optim
using ThreeBodyDecay
using Measurements
using JSON
using Parameters
using RecipesBase
using LaTeXStrings


export mγ
export mπ⁰, mπ⁺
export mD⁰, mD⁺
export mDˣ⁺, ΓDˣ⁺
export mDˣ⁰, ΓDˣ⁰
# 
export μDˣ⁺D⁰
# 
export e2m, m2e
export fm_times_mev
# 
export ΔNLL_90CL, ΔNLL_95CL, ΔNLL_68CL
# 
export tophysicsunits
# 
export k3b, kNR
include("constants.jl")


import Base.^
^(ms::NamedTuple{(:m1,:m2,:m3),T} where T, n::Int) = Tuple(ms).^n

include("covariant_exptessions.jl")

export BW, BW_norm, BW_Swave
export Jᴵ, Jᴵᴵ
include("isobar.jl")

export σ2of3_pm, σ3of1_pm, σ3of2_pm
export σ3of1, σ2of1
export AbstractDalitzMapping
export LinearDalitzMapping, HookSqrtDalitzMapping
export mapdalitz
include("integrationpaths.jl")

export ρ_thr
export masses
export decay_matrix_element_squared
include("abstractxdd.jl")

export ρ_tb
export πDD, γDD
export branch_points
export constructchannel
include("mainmodel.jl")

export DˣD
include("singlechannel.jl")

export ChannelWithIntegrationMethod
include("leftrightmodel.jl")

export interpolated
include("intepolate_merge.jl")

export dispersive
include("dispersion.jl")

export Amplitude
export pole_position
export intensity
export denominator_I, denominator_II
export NonRelBW
include("denominator.jl")

export projecttocos23
export projectto1, projectto3
include("projections.jl")


# serilization
# export obj2nt
# 
export writejson, readjson
export writetoml, readtoml
export transformdictrecursively!
export ifstringgivemeasurement
export ifmeasurementgivestring
export d2nt
include("io.jl")


export dalitzplot
export dpdzspectum, dzdzspectum
export dzdzpipspectum 
include("plotrecipes.jl")


export circleintegral
export cauchy, cauchy′, cauchy′′
export effectiverangeexpansion
include("cauchyintegrals.jl")

end