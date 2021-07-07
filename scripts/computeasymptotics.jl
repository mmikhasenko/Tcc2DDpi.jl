using X2DDpi
using Parameters
using Measurements

settings = transformdictrecursively!(readjson("settings.json"), ifstringgivemeasurement)
@unpack δm0 = settings["fitresults"]

channels = [
    πDD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), (m=mDˣ⁺,Γ=ΓDˣ⁺), (m=mDˣ⁺,Γ=ΓDˣ⁺)),
    πDD((m1=mπ⁰,m2=mD⁺,m3=mD⁰), (m=mDˣ⁺,Γ=ΓDˣ⁺), (m=mDˣ⁰,Γ=ΓDˣ⁰)),
    γDD((m1=mγ, m2=mD⁺,m3=mD⁰), (m=mDˣ⁺,Γ=ΓDˣ⁺), (m=mDˣ⁰,Γ=ΓDˣ⁰))]

# Table D.2
ρ_thr_m0     = sum(ρ_thr.(channels, δm0.val))
ρ_thr_zero   = sum(ρ_thr.(channels, 0))
ρ_thr_second = sum(ρ_thr.(channels, m2e(mDˣ⁰+mD⁺))) 


Φ2_p0(e) = sqrt(X2DDpi.λ(e2m(e)^2,mDˣ⁺^2,mD⁰^2))/e2m(e)^2
Φ2_0p(e) = sqrt(X2DDpi.λ(e2m(e)^2,mDˣ⁺^2,mD⁺^2))/e2m(e)^2


m_matching = 3.9 # GeV
ρ_thr_mpoint = sum(ρ_thr.(channels, m2e(m_matching))) 
ρInf = ρ_thr(channels[1], m2e(m_matching)) / Φ2_p0(m2e(m_matching)) +
       ρ_thr(channels[2], m2e(m_matching)) / Φ2_0p(m2e(m_matching)) +
       ρ_thr(channels[3], m2e(m_matching)) / Φ2_p0(m2e(m_matching))
#

# save results to a file
writejson(joinpath("results","nominal","values_of_rho.json"),
        Dict{Symbol,Any}(
            :rho_thr => Dict{Symbol,Any}(
                :at_m0               => ρ_thr_m0,
                :at_zero             => ρ_thr_zero,
                :at_second_threshold => ρ_thr_second,
                :at_matching_point => ρ_thr_mpoint,
                :forced_asymptotics => ρInf,
                #
                :by_channel_at_m0 => ρ_thr.(channels, 0)
            ),
            :additional => Dict{Symbol,Any}(
                :matching_point => m_matching,
                :matching_point => m_matching,
                :Φ2_p0 => Φ2_p0(m2e(m_matching)),
                :Φ2_0p => Φ2_0p(m2e(m_matching)),
                :forced_asymptotics_p0 => ρ_thr_mpoint / Φ2_p0(m2e(m_matching))
            ),
        )
    )
#