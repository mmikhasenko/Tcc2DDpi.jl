module X2DDpi

using QuadGK
using Cuba
using Interpolations
using Optim
# using ThreeBodyDecay
using Measurements
using JSON


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
export ΔNLL_90CL, ΔNLL_95CL, ΔNLL_66CL
include("constants.jl")


import Base.^
^(ms::NamedTuple{(:m1,:m2,:m3),T} where T, n::Int) = Tuple(ms).^n

include("covariant_exptessions.jl")

export J_I, J_II
include("isobar.jl")
export πDD, γDD
export branch_points
include("abstractxdd.jl")

export ρ_thr
include("dalitz_integral.jl")

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

export writejson, readjson
export writetoml, readtoml
export transformdictrecursively!
export ifstringgivemeasurement
export ifmeasurementgivestring
export d2nt
include("io.jl")

end