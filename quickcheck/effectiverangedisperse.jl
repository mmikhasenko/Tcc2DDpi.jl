using X2DDpi
using Parameters


settings = transformdictrecursively!(readjson("settings.json"), ifstringgivemeasurement)
@unpack δm0 = settings["fitresults"]
const δm0_val = δm0.val
@unpack cutoff, estep = settings["phspmatching"]


ch = DˣD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), BW(m=mDˣ⁺, Γ=ΓDˣ⁺))
ich = interpolated(ch, cutoff; estep=estep)
const A₀ = Amplitude(ich)


const Eᵦ = X2DDpi.Eᵦˣ⁺
@time effrangepars = # 12s
    effectiverangeexpansion(
        # Δe->-1im*ρ_thr(channel(ΓDˣ⁺), Eᵦ+Δe),
        Δe->denominator_II(A₀, Eᵦ+Δe, δm0_val),
        Δe->k3b(Eᵦ+Δe),
        abs(imag(Eᵦ))/20)
# 
tophysicsunits(effrangepars)
# (a_fm = -7.552369007612229,
#  r_fm = 0.16551794481511717 - 0.021681725310022524im)
