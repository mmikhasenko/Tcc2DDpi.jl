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

writejson(joinpath("results","nominal","pole_interpolation.json"), transformdictrecursively!(
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
#######################################################################

# save results to a file
writejson(joinpath("results","nominal","pole.json"), transformdictrecursively!(
        Dict{Symbol,Any}(
            :size_estimate => Dict{Symbol,Any}(
                :momentum => γ,
                :R_eff_ΔE => R_ΔE,
            ),
            :effective_range_parameters => Dict{Symbol,Any}(
                :Re_inv_scatt_length => real(inverse_scattering_length),
                :Im_inv_scatt_length => imag(inverse_scattering_length),
                :Re_scatt_length => real(scattering_length),
                :Im_scatt_length => imag(scattering_length),
                :effective_range => (; r_90CL, r_95CL),
                :g_coupling => (; g_90CL, g_95CL),
                :technical => Dict{Symbol,Any}(
                    :rho_inf => ρInf,
                    :w_matching => w_matching)
                ),
            :pole_position => Dict{Symbol,Any}(
                :pole_sv => pole_sv),
        ), ifmeasurementgivestring)
    )
#

#######################################################################

#
ev = range(-1,3.0,length=100)
D_advans_δm0(e) = denominator_I(ampX0, e, δm0.val) / (w_matching * ρInf)
D_nonrel_δm0(e) = denominator_I(NonRelBW(), e+1e-6im, δm0.val)
# 
invA_nonrel, invA_advans = D_nonrel_δm0.(ev), D_advans_δm0.(ev)

writejson(joinpath("results","nominal","inverse_amplitude.json"),
        Dict(
            :inverse_amplitude => Dict{Symbol,Any}(
                :mgrid => ev,
                :invA_nonrel_real => real.(invA_nonrel),
                :invA_nonrel_imag => imag.(invA_nonrel),
                :invA_advans_real => real.(invA_advans),
                :invA_advans_imag => imag.(invA_advans),
            ),
        )
    )
#
