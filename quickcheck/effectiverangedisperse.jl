using X2DDpi
using Parameters



using Plots
using LaTeXStrings
theme(:wong2, frame=:box, grid=false, minorticks=true, 
    guidefontvalign=:top, guidefonthalign=:right,
    xlim=(:auto,:auto), ylim=(:auto,:auto), lw=1.2)


settings = transformdictrecursively!(readjson("settings.json"), ifstringgivemeasurement)
@unpack δm0 = settings["fitresults"]
const δm0_val = δm0.val
@unpack cutoff, estep = settings["phspmatching"]



# DˣD
const A₀_DˣD = Amplitude(
    interpolated(
        DˣD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), BW(m=mDˣ⁺, Γ=ΓDˣ⁺)),
        cutoff; estep=estep))
# 
plot(e->abs2(denominator_II(A₀_DˣD, e, δm0_val)), -1, 1)
plot(Δe->abs2(denominator_II(A₀_DˣD, Eᵦˣ⁺ + Δe, δm0_val)), -1, 1)

let
    @time effrangepars = # 12s
        effectiverangeexpansion(
            Δe->denominator_II(A₀_DˣD, Eᵦˣ⁺+Δe, δm0_val),
            Δe->k3b(Eᵦˣ⁺+Δe),
            ComplexBranchPointExpansion(CircularSum(abs(imag(Eᵦˣ⁺))/20, 50)))
    # 
    (; tophysicsunits(effrangepars)..., effrangepars...)
end

# exact:
# (a_fm = -7.552369007612229,
#  r_fm = 0.16551794481511717 - 0.021681725310022524im)

# N=50
# (a_fm = -7.552369007612229,
#  r_fm = 0.16551794481511717 - 0.021681725310022524im)



# πDD
const A₀_πDD = let ch = πDD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), BW(m=mDˣ⁺,Γ=ΓDˣ⁺), BW(m=mDˣ⁺,Γ=ΓDˣ⁺))

    set3 = ChannelWithIntegrationMethod(ch, HookSqrtDalitzMapping{3}())
    set2 = ChannelWithIntegrationMethod(ch, HookSqrtDalitzMapping{2}())
    # 
    ich2 = interpolated(set2, cutoff; estep=estep)
    ich3 = interpolated(set3, cutoff; estep=estep)
    # 
    Amplitude((ich2,ich3))
end

plot(e->abs2(denominator_II(A₀_πDD, e, δm0_val)), -1, 1)
plot(Δe->abs2(denominator_II(A₀_πDD, Eᵦˣ⁺ + Δe, δm0_val)), -1, 1)


# D₀ = denominator_I(A₀_πDD, 0.0, δm0_val)
# ich = interpolated(ch, cutoff; estep=estep)
# Aₓ = Amplitude(ich)
# Dₓ = denominator_I(Aₓ, 0.0, δm0_val)
# # test
# abs((Dₓ - D₀) / Dₓ) < 1e-4


let
    @time effrangepars = # 12s
        effectiverangeexpansion(
            Δe->denominator_II(A₀_πDD, Eᵦˣ⁺+Δe, δm0_val),
            Δe->k3b(Eᵦˣ⁺+Δe),
            ComplexBranchPointExpansion(CircularSum(abs(imag(Eᵦˣ⁺))/20, 150)))
    # 
    (; tophysicsunits(effrangepars)..., effrangepars...)
end

# N = 250
# (a_fm = -7.3176188107221085, r_fm = -0.7794936309216776 + 0.2597180191937492im)

# N = 150
# (a_fm = -7.317634027931594, r_fm = -0.779907328004435 + 0.260605536433528im)

# N = 50
# (a_fm = -7.317640051885012, r_fm = -0.7797722395543267 + 0.26066081944521835im)

let
    Δev = range(-0.1, 0.1, length=100)
    calv1 = [denominator_II(A₀_πDD, Eᵦˣ⁺+Δe, δm0_val) for Δe in Δev]
    calv2 = [denominator_II(A₀_DˣD, Eᵦˣ⁺+Δe, δm0_val) for Δe in Δev]
    calv3 = [-1im*k3b(Eᵦˣ⁺+Δe) for Δe in Δev]
    # 
    calv1, calv2, calv3 = imag.(calv1), imag.(calv2), imag.(calv3)
    # 
    plot(Δev, [calv1 ./ (calv1[end])  calv2 ./ (calv2[end])  calv3 ./ (calv3[end])],
        lab=["πDD" "DˣD" "k"])
end




# πDD+πDD+γDD
@time const A₀_full = let 
    ch1 = πDD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), BW(m=mDˣ⁺,Γ=ΓDˣ⁺), BW(m=mDˣ⁺,Γ=ΓDˣ⁺))
    ch2 = πDD((m1=mπ⁰,m2=mD⁺,m3=mD⁰), BW(m=mDˣ⁺,Γ=ΓDˣ⁺), BW(m=mDˣ⁰,Γ=ΓDˣ⁰))
    ch3 = γDD((m1=mγ, m2=mD⁺,m3=mD⁰), BW(m=mDˣ⁺,Γ=ΓDˣ⁺), BW(m=mDˣ⁰,Γ=ΓDˣ⁰))
    # 
    iπDD2 = interpolated(
        ChannelWithIntegrationMethod(ch1, HookSqrtDalitzMapping{2}()),
        cutoff; estep=estep)
    iπDD3 = interpolated(
        ChannelWithIntegrationMethod(ch1, HookSqrtDalitzMapping{3}()),
        cutoff; estep=estep)
    # 
    iπDD2′ = interpolated(
        ChannelWithIntegrationMethod(ch2, HookSqrtDalitzMapping{2}()),
        cutoff; estep=estep)
    iπDD3′ = interpolated(
        ChannelWithIntegrationMethod(ch2, HookSqrtDalitzMapping{3}()),
        cutoff; estep=estep)
    # 
    iγDD2 = interpolated(
        ChannelWithIntegrationMethod(ch3, HookSqrtDalitzMapping{2}()),
        cutoff; estep=estep)
    iγDD3 = interpolated(
        ChannelWithIntegrationMethod(ch3, HookSqrtDalitzMapping{3}()),
        cutoff; estep=estep)
    # 
    Amplitude((iπDD2,iπDD3,iπDD2′,iπDD3′,iγDD2,iγDD3))
end

plot(e->abs2(denominator_II(A₀_full, e, δm0_val)), -1, 1)
plot(Δe->abs2(denominator_II(A₀_full, Eᵦˣ⁺ + Δe, δm0_val)), -1, 1)

# D₀ = denominator_I(A₀_πDD, 0.0, δm0_val)
# ich = interpolated(ch, cutoff; estep=estep)
# Aₓ = Amplitude(ich)
# Dₓ = denominator_I(Aₓ, 0.0, δm0_val)
# # test
# abs((Dₓ - D₀) / Dₓ) < 1e-4

let
    @time effrangepars = # 12s
        effectiverangeexpansion(
            Δe->denominator_II(A₀_full, Eᵦˣ⁺+Δe, δm0_val),
            Δe->k3b(Eᵦˣ⁺+Δe),
            ComplexBranchPointExpansion(CircularSum(abs(imag(Eᵦˣ⁺))/20, 50)))
    # 
    (; tophysicsunits(effrangepars)..., effrangepars...)
end

# N = 50
# (a_fm = -5.959338885798581, r_fm = -4.361863371414783 + 0.460836323360997im)

# N = 150
# (a_fm = -5.959335691531121, r_fm = -4.36195566780003 + 0.4608003964193325im)
