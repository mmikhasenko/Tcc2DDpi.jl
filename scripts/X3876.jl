using X2DDpi
using Parameters
using Measurements
using Interpolations
# using Plots

settings = transformdictrecursively!(readjson("settings.json"), ifstringgivemeasurement)

@unpack cutoff, estep = settings["phspmatching"]
@unpack dm_min, dm_max, dm_N = settings["polepositiongrid"]
@unpack δm_sw, δm_ΔM = settings["fitresults"]
# 
channels = [
    πDD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), (m=mDˣ⁺,Γ=ΓDˣ⁺), (m=mDˣ⁺,Γ=ΓDˣ⁺)),
    πDD((m1=mπ⁰,m2=mD⁺,m3=mD⁰), (m=mDˣ⁺,Γ=ΓDˣ⁺), (m=mDˣ⁰,Γ=ΓDˣ⁰)),
    γDD((m1=mγ, m2=mD⁺,m3=mD⁰), (m=mDˣ⁺,Γ=ΓDˣ⁺), (m=mDˣ⁰,Γ=ΓDˣ⁰))]
#
ichannels = interpolated.(channels, 4.0) # cutoff

# pole
δmv = range(dm_min, dm_max, length=dm_N)
sampledpp = [pole_position(Tuple(ichannels),δm) for δm in δmv]

itr_m, itr_Γ =
	interpolate((δmv,), getproperty.(sampledpp, :m_pole), Gridded(Linear())),
	interpolate((δmv,), 2 .* getproperty.(sampledpp, :half_Γ_pole), Gridded(Linear()))
# 
pole_sv = NamedTuple{(:m_pole, :Γ_pole)}([itr_m(δm_sw), itr_Γ(δm_sw)])
pole_ΔM = NamedTuple{(:m_pole, :Γ_pole)}([itr_m(δm_ΔM), itr_Γ(δm_ΔM)])


# scattering length

ρInf = sum(ich.cutoffratio for ich in ichannels)
w_matching = ichannels[1].cutoffratio*3/2 / ρInf * 2/e2m(0)
scattering_length = denominator_I(Tuple(ichannels), 0.0, δm_ΔM) / ρInf / w_matching
# 0.196/real(scattering_length)

# effective range

# save
writejson(joinpath("results","nominal","pole.json"), transformdictrecursively!(
        Dict{Symbol,Any}(
            :pole_position => Dict{Symbol,Any}(
                :pole_sv => pole_sv,
                :pole_ΔM => pole_ΔM),
            :pole_interpolation => Dict{Symbol,Any}(
                :grid => δmv,
                :values => sampledpp
            )
        ), ifmeasurementgivestring)
    )
# 