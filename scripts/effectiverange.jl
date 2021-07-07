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

@unpack Γ0_90CL_syst, Γ0_95CL_syst = settings["fitresults"]
@unpack Γ0_90CL_stat, Γ0_95CL_stat = settings["fitresults"]

# retrieve the model
modelDict = readjson(joinpath("results","nominal","model.json"))
ichannels = interpolated.(d2nt.(modelDict["ichannels"]))
channels = getproperty.(ichannels, :channel)

ampX0 = Amplitude(Tuple(ichannels), zero) # zero is 1/Γ

# naive size estimation
γ = sqrt(-2μDˣ⁺D⁰*1e3*δm0)
R_ΔE = fm_times_mev/γ

# scattering length
ρInf = sum(ich.cutoffratio for ich in ichannels)
w_matching = ichannels[1].cutoffratio*3/2 / ρInf * 2/e2m(0) # 1/GeV
inverse_scattering_length = denominator_I(ampX0, 0.0, δm0) / (w_matching * ρInf) * 1e3 # MeV
scattering_length = fm_times_mev/inverse_scattering_length

# effective range
#          1 / (GeV * MeV) / (1/GeV)      = 1/MeV
reff(Γ) = 8 / (e2m(0) * Γ) / (w_matching) # 1/MeV

g_90CL, g_95CL = sqrt.((2*e2m(0)).*[Γ0_90CL_syst, Γ0_95CL_syst] .* 1e-3 /ρInf)  # GeV
r_90CL, r_95CL = round.(reff.([Γ0_90CL_syst, Γ0_95CL_syst]) * fm_times_mev; digits=3)  # fm

# save results to a file
writejson(joinpath("results","nominal","effective_range.json"), transformdictrecursively!(
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
        ), ifmeasurementgivestring)
    )
#

#######################################################################

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
