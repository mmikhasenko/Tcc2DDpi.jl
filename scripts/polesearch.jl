using Pkg
Pkg.activate(joinpath(@__DIR__,".."))
# 
using X2DDpi
using Parameters
using Measurements
using Interpolations
using Statistics

settings = transformdictrecursively!(readjson("settings.json"), ifstringgivemeasurement)
#
@unpack δm0 = settings["fitresults"]
@unpack dm_min, dm_max, dm_N = settings["polepositiongrid"]

@unpack Γ0_68CL_syst, Γ0_90CL_syst, Γ0_95CL_syst = settings["fitresults"]
@unpack Γ0_68CL_stat, Γ0_90CL_stat, Γ0_95CL_stat = settings["fitresults"]

# retrieve the model
modelDict = readjson(joinpath("results","nominal","model.json"))
ichannels = interpolated.(d2nt.(modelDict["ichannels"]))
channels = getproperty.(ichannels, :channel)

# 
ρInf = sum(ich.cutoffratio for ich in ichannels)

ampX0 = Amplitude(Tuple(ichannels), zero)
#
# pole
δmv = range(dm_min, dm_max, length=20) # dm_N 
ppsampled = [(@show δm; pole_position(ampX0,δm)) for δm in δmv]

itr_m, itr_Γ =
	interpolate((δmv,), getproperty.(ppsampled, :m_pole), Gridded(Linear())),
	interpolate((δmv,), 2 .* getproperty.(ppsampled, :half_Γ_pole), Gridded(Linear()))
# 
pole_sv = NamedTuple{(:m_pole, :Γ_pole)}([itr_m(δm0), itr_Γ(δm0)])

writejson(joinpath("results","nominal","pole_default.json"), transformdictrecursively!(
        Dict{Symbol,Any}(
            :pole_position => Dict{Symbol,Any}(
                :pole_sv => pole_sv
            ),
            :pole_interpolation => Dict{Symbol,Any}(
                :mgrid => δmv,
                :ppvalues => ppsampled
            ),
        ), ifmeasurementgivestring)
    )
# 

#######################################################################
ampX(Γ) = Amplitude(Tuple(ichannels),
    e->(e2m(δm0.val)^2-e2m(e)^2)/(e2m(δm0.val)*Γ/1e3/ρInf))

ppsampled_stat_syst = [(@show δm; pole_position(ampX(Γ),δm))
    for δm in δmv, Γ in [Γ0_68CL_stat, Γ0_90CL_stat, Γ0_95CL_stat,
                         Γ0_68CL_syst, Γ0_90CL_syst, Γ0_95CL_syst]]
# 
writejson(joinpath("results","nominal","pole_interpolation_stat_syst.json"),
        Dict(
            :pole_interpolation => Dict{Symbol,Any}(
                :mgrid => δmv,
                :ppvalues_68CL_stat => ppsampled_stat_syst[:,1],
                :ppvalues_90CL_stat => ppsampled_stat_syst[:,2],
                :ppvalues_95CL_stat => ppsampled_stat_syst[:,3],
                :ppvalues_68CL_syst => ppsampled_stat_syst[:,4],
                :ppvalues_90CL_syst => ppsampled_stat_syst[:,5],
                :ppvalues_95CL_syst => ppsampled_stat_syst[:,6],
            ),
        )
    )
# 