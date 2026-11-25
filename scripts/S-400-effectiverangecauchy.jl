using Pkg
cd(joinpath(@__DIR__, ".."))
Pkg.activate(".")
Pkg.instantiate()
# 
using X2DDpi
using Parameters
using Measurements
using Interpolations
using Statistics
using LaTeXStrings
using DataFrames
using Markdown

using Plots
theme(:wong2, size=(500, 350), minorticks=true, grid=false, frame=:box,
    guidefontvalign=:top, guidefonthalign=:right,
    foreground_color_legend=nothing,
    legendfontsize=9, legend=:topright,
    xlim=(:auto, :auto), ylim=(:auto, :auto))


#        _|              _|                
#    _|_|_|    _|_|_|  _|_|_|_|    _|_|_|  
#  _|    _|  _|    _|    _|      _|    _|  
#  _|    _|  _|    _|    _|      _|    _|  
#    _|_|_|    _|_|_|      _|_|    _|_|_|  

settings = transformdictrecursively!(readjson("settings.json"), ifstringgivemeasurement)
#
@unpack δm0 = settings["fitresults"]
const δm0_val = δm0.val
@unpack e_cutoff, e_step = settings["phspmatching"]
m_cutoff = e2m(e_cutoff)


#                            _|            
#    _|_|_|    _|_|      _|_|_|    _|_|    
#  _|        _|    _|  _|    _|  _|_|_|_|  
#  _|        _|    _|  _|    _|  _|        
#    _|_|_|    _|_|      _|_|_|    _|_|_|  


# πDD + πDD + γDD: Dˣ⁺ + Dˣ⁺ + Dˣ⁰
"""
The default full model.
In includes three three-body cuts, π⁺D⁰D⁰ + π⁰D⁰D⁺ + γD⁰D⁺ and all possible resonances:
 - two Dˣ⁺ in π⁺D⁰D⁰
 - Dˣ⁺ and Dˣ⁰ in π⁰D⁰D⁺
 - Dˣ⁺ and Dˣ⁰ in γD⁰D⁺
"""
const model0 = let
    ch1 = πDD((m1=mπ⁺, m2=mD⁰, m3=mD⁰), BW(m=mDˣ⁺, Γ=ΓDˣ⁺), BW(m=mDˣ⁺, Γ=ΓDˣ⁺))
    ch2 = πDD((m1=mπ⁰, m2=mD⁺, m3=mD⁰), BW(m=mDˣ⁺, Γ=ΓDˣ⁺), BW(m=mDˣ⁰, Γ=ΓDˣ⁰))
    ch3 = γDD((m1=mγ, m2=mD⁺, m3=mD⁰), BW(m=mDˣ⁺, Γ=ΓDˣ⁺), BW(m=mDˣ⁰, Γ=ΓDˣ⁰))
    # 
    @time iπDD2 = interpolated(
        ChannelWithIntegrationMethod(ch1, HookSqrtDalitzMapping{2}()),
        m_cutoff; mstep = e_step / 1e3)
    @time iπDD3 = interpolated(
        ChannelWithIntegrationMethod(ch1, HookSqrtDalitzMapping{3}()),
        m_cutoff; mstep = e_step / 1e3)
    # 
    @time iπDD2′ = interpolated(
        ChannelWithIntegrationMethod(ch2, HookSqrtDalitzMapping{2}()),
        m_cutoff; mstep = e_step / 1e3)
    @time iπDD3′ = interpolated(
        ChannelWithIntegrationMethod(ch2, HookSqrtDalitzMapping{3}()),
        m_cutoff; mstep = e_step / 1e3)
    # 
    @time iγDD2 = interpolated(
        ChannelWithIntegrationMethod(ch3, HookSqrtDalitzMapping{2}()),
        m_cutoff; mstep = e_step / 1e3)
    @time iγDD3 = interpolated(
        ChannelWithIntegrationMethod(ch3, HookSqrtDalitzMapping{3}()),
        m_cutoff; mstep = e_step / 1e3)
    # 
    Amplitude((iπDD2, iπDD3, iπDD2′, iπDD3′, iγDD2, iγDD3))
end


# πDD + πDD + γDD: Dˣ⁺
"""
The model with the Dˣ⁺ resonance only.
However, it includes three three-body cuts, π⁺D⁰D⁰ + π⁰D⁰D⁺ + γD⁰D⁺:
 - two Dˣ⁺ in π⁺D⁰D⁰
 - one Dˣ⁺ in π⁰D⁰D⁺
 - one Dˣ⁺ in γD⁰D⁺
"""
const model1 = let
    # ch1 = πDD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), BW(m=mDˣ⁺,Γ=ΓDˣ⁺), BW(m=mDˣ⁺,Γ=ΓDˣ⁺))
    ch2 = πDD((m1=mπ⁰, m2=mD⁺, m3=mD⁰), BW(m=mDˣ⁺, Γ=ΓDˣ⁺), ZeroBW(m=mDˣ⁰, Γ=ΓDˣ⁰))
    ch3 = γDD((m1=mγ, m2=mD⁺, m3=mD⁰), BW(m=mDˣ⁺, Γ=ΓDˣ⁺), ZeroBW(m=mDˣ⁰, Γ=ΓDˣ⁰))
    # 
    # iπDD2 = interpolated(
    #     ChannelWithIntegrationMethod(ch1, HookSqrtDalitzMapping{2}()),
    #     m_cutoff; mstep = e_step / 1e3)
    # iπDD3 = interpolated(
    #     ChannelWithIntegrationMethod(ch1, HookSqrtDalitzMapping{3}()),
    #     m_cutoff; mstep = e_step / 1e3)
    iπDD2, iπDD3 = model0.ik[1], model0.ik[2]
    # 
    @time iπDD3′ = interpolated(
        ChannelWithIntegrationMethod(ch2, HookSqrtDalitzMapping{3}()),
        m_cutoff; mstep = e_step / 1e3)
    # 
    @time iγDD3 = interpolated(
        ChannelWithIntegrationMethod(ch3, HookSqrtDalitzMapping{3}()),
        m_cutoff; mstep = e_step / 1e3)
    # 
    Amplitude((iπDD2, iπDD3, iπDD3′, iγDD3))
end

# πDD: Dˣ⁺ + Dˣ⁺
"""
The model with the π⁺D⁰D⁰ system with two Dˣ⁺ resonances.
"""
const model2 = let
    # ch1 = πDD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), BW(m=mDˣ⁺,Γ=ΓDˣ⁺), BW(m=mDˣ⁺,Γ=ΓDˣ⁺))
    iπDD2, iπDD3 = model0.ik[1], model0.ik[2]
    # 
    # set3 = ChannelWithIntegrationMethod(ch, HookSqrtDalitzMapping{3}())
    # ich3 = interpolated(set3, m_cutoff; mstep = e_step / 1e3)
    # 
    Amplitude((iπDD2, iπDD3))
end

# πDD: Dˣ⁺
"""
The model with the π⁺D⁰D⁰ system with a single, unsymmetrized Dˣ⁺ resonance.
"""
const model3 = let
    ch = πDD((m1=mπ⁺, m2=mD⁰, m3=mD⁰), BW(m=mDˣ⁺, Γ=ΓDˣ⁺), ZeroBW(m=mDˣ⁺, Γ=ΓDˣ⁺))
    # 
    set3 = ChannelWithIntegrationMethod(ch, HookSqrtDalitzMapping{3}())
    ich3 = interpolated(set3, m_cutoff; mstep = e_step / 1e3)
    # 
    Amplitude((ich3,))
end

# πDD: Dˣ⁺
const model4 = Amplitude(
    interpolated(
        DˣD((m1=mπ⁺, m2=mD⁰, m3=mD⁰), BW(m=mDˣ⁺, Γ=ΓDˣ⁺)),
        m_cutoff; mstep = e_step / 1e3))
# 
# 
modelnames = [:model0, :model1, :model2, :model3, :model4]
const Nm = length(modelnames)



# let
#     @time effrangepars = # 12s
#         effectiverangeexpansion(
#             Δe -> denominator_II(model0, Eᵦˣ⁺ + Δe, δm0_val),
#             Δe -> k3b(Eᵦˣ⁺ + Δe),
#             ComplexBranchPointExpansion(CircularSum(abs(imag(Eᵦˣ⁺)) / 20, 50)))
#     # 
#     (; tophysicsunits(effrangepars)..., effrangepars...)
# end


df = DataFrame(; modelnames)
df.model = eval.(df.modelnames)
# 
df1 = DataFrame(zeros(Complex{Float64}, Nm, 4 * 3),
    [:a⁻¹, :r, :r_fm, :N,
        :a⁻¹⁺ᵟ, :r⁺ᵟ, :r_fm⁺ᵟ, :N⁺ᵟ,
        :a⁻¹⁻ᵟ, :r⁻ᵟ, :r_fm⁻ᵟ, :N⁻ᵟ])
# 
df2 = DataFrame(
    zeros(Float64, Nm, 1 * 3),
    [:a_fm, :a_fm⁺ᵟ, :a_fm⁻ᵟ,])
# 
df = hcat(df, df1, df2)

let Nappr = 150, Rexp = abs(imag(Eᵦˣ⁺)) / 20
    # 
    for i in 1:Nm
        𝒜 = df.model[i]
        @time effrangepars = # 12s
            effectiverangeexpansion(
                Δe -> denominator_II(𝒜, Eᵦˣ⁺ + Δe, δm0_val),
                Δe -> k3b(Eᵦˣ⁺ + Δe),
                ComplexBranchPointExpansion(CircularSum(Rexp, Nappr)))
        # 
        efrpars = (; tophysicsunits(effrangepars)..., effrangepars...)
        # 
        df.a_fm[i] = efrpars.a_fm
        df.r_fm[i] = efrpars.r_fm
        df.a⁻¹[i] = efrpars.a⁻¹
        df.r[i] = efrpars.r
        df.N[i] = efrpars.N
    end
end

let Nappr = 150, Rexp = abs(imag(Eᵦˣ⁺)) / 20
    # 
    for i in 1:Nm
        𝒜 = df.model[i]
        # 
        # dm0+δdm0
        @time effrangepars⁺ᵟ = # 12s
            effectiverangeexpansion(
                Δe -> denominator_II(𝒜, Eᵦˣ⁺ + Δe, δm0_val + δm0.err),
                Δe -> k3b(Eᵦˣ⁺ + Δe),
                ComplexBranchPointExpansion(CircularSum(Rexp, Nappr)))
        # 
        efrpars⁺ᵟ = (; tophysicsunits(effrangepars⁺ᵟ)..., effrangepars⁺ᵟ...)
        # 
        df.a_fm⁺ᵟ[i] = efrpars⁺ᵟ.a_fm
        df.r_fm⁺ᵟ[i] = efrpars⁺ᵟ.r_fm
        df.a⁻¹⁺ᵟ[i] = efrpars⁺ᵟ.a⁻¹
        df.r⁺ᵟ[i] = efrpars⁺ᵟ.r
        df.N⁺ᵟ[i] = efrpars⁺ᵟ.N
        #
        # dm0-δdm0
        @time effrangepars⁻ᵟ = # 12s
            effectiverangeexpansion(
                Δe -> denominator_II(𝒜, Eᵦˣ⁺ + Δe, δm0_val - δm0.err),
                Δe -> k3b(Eᵦˣ⁺ + Δe),
                ComplexBranchPointExpansion(CircularSum(Rexp, Nappr)))
        # 
        efrpars⁻ᵟ = (; tophysicsunits(effrangepars⁻ᵟ)..., effrangepars⁻ᵟ...)
        # 
        df.a_fm⁻ᵟ[i] = efrpars⁻ᵟ.a_fm
        df.r_fm⁻ᵟ[i] = efrpars⁻ᵟ.r_fm
        df.a⁻¹⁻ᵟ[i] = efrpars⁻ᵟ.a⁻¹
        df.r⁻ᵟ[i] = efrpars⁻ᵟ.r
        df.N⁻ᵟ[i] = efrpars⁻ᵟ.N
    end
end


# test
summarized_columns = vcat((["a_fm" "r_fm"] .* ["", "⁺ᵟ", "⁻ᵟ"])...)
df_summary = select(df, :modelnames,
    summarized_columns .=> ByRow(x -> round(x; digits=2)) .=> summarized_columns .* "_round")
print(df_summary)

md"""
- The effective range is positive +0.15 fm for the model with one Dˣ⁺ resonance (model3).
The positive slope comes from the two-body loop integral. 
- Once the OPE is included (model2), | Dˣ⁺ + Dˣ⁺ |², the effective range turns negative, -0.78 fm
- The next step is to add two other three-body thresholds, D⁰D⁺π⁰ and D⁰D⁺γ, where the Dˣ⁺ can decay, the effective range changes unsignificantly, -0.63 fm
- The full model is obtained by adding the Dˣ⁰D⁺ threshold, the value of the effective range is shifted to -4.36 fm
"""

writejson(joinpath("results", "nominal", "effective-range-table.json"),
    hcat(df_summary, select(df, Not([:modelnames, :model]))))
#

