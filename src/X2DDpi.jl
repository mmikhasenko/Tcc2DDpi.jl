module X2DDpi

using QuadGK
using Cuba
using Interpolations
using Optim
# using ThreeBodyDecay
using Measurements
using JSON


export mπ⁰, mπ⁺
export mD⁰, mD⁺
export mDˣ⁺, ΓDˣ⁺
export mDˣ⁰, ΓDˣ⁰
export mγ

export e2m, m2e
include("constants.jl")


import Base.^
^(ms::NamedTuple{(:m1,:m2,:m3),T} where T, n::Int) = Tuple(ms).^n

include("covariant_exptessions.jl")

export J_I, J_II
include("isobar.jl")
export πDD, γDD
include("abstractxdd.jl")

export ρ_thr
include("dalitz_integral.jl")

export interpolated
include("intepolate_merge.jl")

export dispersive
include("dispersion.jl")

export pole_position
export intensity
export denominator_I, denominator_II
include("denominator.jl")

export writejson, readjson
export writetoml, readtoml
export transformdictrecursively!
export ifstringgivemeasurement
export ifmeasurementgivestring
include("io.jl")

end