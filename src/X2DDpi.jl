module X2DDpi

using QuadGK
using Cuba
using Interpolations
using Optim
using ThreeBodyDecays
using Measurements
using JSON
using Parameters
using RecipesBase
using LaTeXStrings

# Include the xDDPhaseSpace module
include("xDDPhaseSpace/xDDPhaseSpace.jl")

# Import specific methods from xDDPhaseSpace module
import .xDDPhaseSpace: AbstractxDD, masses
import .xDDPhaseSpace: decay_matrix_element_squared
import .xDDPhaseSpace: integrand_mapped_thr, ρ_thr
import .xDDPhaseSpace: mapdalitz, mapdalitzmethod

# serialization
import .xDDPhaseSpace: obj2nt

# model
import .xDDPhaseSpace: πDD, branch_points, constructchannel

# analytic continuation
import .xDDPhaseSpace: Jᴵ, Jᴵᴵ

# kinematics
import .xDDPhaseSpace: e2m
import .xDDPhaseSpace: σ2of3_pm, σ3of1_pm, σ3of2_pm
import .xDDPhaseSpace: σ3of1, σ2of1

# 
export mγ
export mπ⁰, mπ⁺
export mD⁰, mD⁺
export mDˣ⁺, ΓDˣ⁺
export mDˣ⁰, ΓDˣ⁰
export thr_D⁰Dˣ⁺
# 
export μDˣ⁺D⁰, Eᵦˣ⁺
# 
export e2m, m2e
export WithThrE
export fm_times_mev
# 
export ΔNLL_90CL, ΔNLL_95CL, ΔNLL_68CL
# 
export tophysicsunits
# 
export k3b, kNR
include("constants.jl")

# Reexport xDDPhaseSpace functionality
export AbstractxDD, ρ_thr, masses, interpolated, dispersive
export πDD, branch_points, constructchannel

# model 1: main
export γDD
include("mainmodel-gamma.jl")

# model 2: single channel
export DˣD
include("singlechannel.jl")

# model 3: leftright model
# integration paths
import .xDDPhaseSpace: AbstractDalitzMapping
import .xDDPhaseSpace: LinearDalitzMapping, HookSqrtDalitzMapping
# 
export ChannelWithIntegrationMethod
include("leftrightmodel.jl")


# advanced functions
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
export cauchyintegral, cauchyintegral′
export cauchyintegral′′, cauchyintegral′′′
export highertermexpansion
export effectiverangeexpansion
export ComplexBranchPointExpansion
export CircularIntegral, CircularSum
include("cauchyintegrals.jl")


export ere, hte
export EffectiveRangeFit
include("fiteffectiverange.jl")


export RegressionResult
export linear_regression
export effective_range_scan
export linear_regression_in_sqrt_ϵf
export effective_range_extrapolation
include("extrapolation.jl")

end