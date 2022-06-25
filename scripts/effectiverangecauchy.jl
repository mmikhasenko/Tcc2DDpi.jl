using Pkg
cd(joinpath(@__DIR__,".."))
Pkg.activate(".")
# 
using X2DDpi
using Parameters
using Measurements
using Interpolations
using Statistics
using LaTeXStrings

using Plots
theme(:wong2, size=(500,350), minorticks=true, grid=false, frame=:box,
    guidefontvalign=:top, guidefonthalign=:right,
    foreground_color_legend = nothing,
    legendfontsize=9, legend =:topright,
    xlim=(:auto,:auto), ylim=(:auto,:auto))

@with_kw struct ZeroBW <: X2DDpi.AbstractLinesShape
    m::Float64
    Γ::Float64
end

X2DDpi.Jᴵ(σ::Number,pars::ZeroBW) = 0.0
X2DDpi.Jᴵᴵ(σ::Number,pars::ZeroBW) = 0.0

#        _|              _|                
#    _|_|_|    _|_|_|  _|_|_|_|    _|_|_|  
#  _|    _|  _|    _|    _|      _|    _|  
#  _|    _|  _|    _|    _|      _|    _|  
#    _|_|_|    _|_|_|      _|_|    _|_|_|  

settings = transformdictrecursively!(readjson("settings.json"), ifstringgivemeasurement)
#
@unpack δm0 = settings["fitresults"]
const δm0_val = δm0.val
@unpack cutoff, estep = settings["phspmatching"]


#                            _|            
#    _|_|_|    _|_|      _|_|_|    _|_|    
#  _|        _|    _|  _|    _|  _|_|_|_|  
#  _|        _|    _|  _|    _|  _|        
#    _|_|_|    _|_|      _|_|_|    _|_|_|  


# πDD + πDD + γDD: Dˣ⁺ + Dˣ⁺ + Dˣ⁰
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


# πDD + πDD + γDD: Dˣ⁺
@time const A₀_full_DˣD = let 
    # ch1 = πDD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), BW(m=mDˣ⁺,Γ=ΓDˣ⁺), BW(m=mDˣ⁺,Γ=ΓDˣ⁺))
    ch2 = πDD((m1=mπ⁰,m2=mD⁺,m3=mD⁰), BW(m=mDˣ⁺,Γ=ΓDˣ⁺), ZeroBW(m=mDˣ⁰,Γ=ΓDˣ⁰))
    ch3 = γDD((m1=mγ, m2=mD⁺,m3=mD⁰), BW(m=mDˣ⁺,Γ=ΓDˣ⁺), ZeroBW(m=mDˣ⁰,Γ=ΓDˣ⁰))
    # 
    # iπDD2 = interpolated(
    #     ChannelWithIntegrationMethod(ch1, HookSqrtDalitzMapping{2}()),
    #     cutoff; estep=estep)
    # iπDD3 = interpolated(
    #     ChannelWithIntegrationMethod(ch1, HookSqrtDalitzMapping{3}()),
    #     cutoff; estep=estep)
    iπDD2, iπDD3 = A₀_full.ik[1], A₀_full.ik[2]
    # 
    iπDD3′ = interpolated(
        ChannelWithIntegrationMethod(ch2, HookSqrtDalitzMapping{3}()),
        cutoff; estep=estep)
    # 
    iγDD3 = interpolated(
        ChannelWithIntegrationMethod(ch3, HookSqrtDalitzMapping{3}()),
        cutoff; estep=estep)
    # 
    Amplitude((iπDD2,iπDD3,iπDD3′,iγDD3))
end

# πDD: Dˣ⁺ + Dˣ⁺
const A₀_π⁺D⁰D⁰ = let
    # ch1 = πDD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), BW(m=mDˣ⁺,Γ=ΓDˣ⁺), BW(m=mDˣ⁺,Γ=ΓDˣ⁺))
    iπDD2, iπDD3 = A₀_full.ik[1], A₀_full.ik[2]
    # 
    # set3 = ChannelWithIntegrationMethod(ch, HookSqrtDalitzMapping{3}())
    # ich3 = interpolated(set3, cutoff; estep=estep)
    # 
    Amplitude((iπDD2,iπDD3))
end

# πDD: Dˣ⁺
const A₀_DˣD = let
    ch = πDD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), BW(m=mDˣ⁺,Γ=ΓDˣ⁺), ZeroBW(m=mDˣ⁺,Γ=ΓDˣ⁺))
    # 
    set3 = ChannelWithIntegrationMethod(ch, HookSqrtDalitzMapping{3}())
    ich3 = interpolated(set3, cutoff; estep=estep)
    # 
    Amplitude((ich3,))
end

# πDD: Dˣ⁺
const A₀′_DˣD = Amplitude(
    interpolated(
        DˣD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), BW(m=mDˣ⁺, Γ=ΓDˣ⁺)),
        cutoff; estep=estep))
# 
# 
modelnames = [:A₀_full, :A₀_full_DˣD, :A₀_π⁺D⁰D⁰, :A₀_DˣD, :A₀′_DˣD]
const Nm = length(modelnames)




using DataFrames
df = DataFrame(; modelnames)
df.model = eval.(df.modelnames)

let 
    @time effrangepars = # 12s
        effectiverangeexpansion(
            Δe->denominator_II(A₀_full, Eᵦˣ⁺+Δe, δm0_val),
            Δe->k3b(Eᵦˣ⁺+Δe),
            ComplexBranchPointExpansion(CircularSum(abs(imag(Eᵦˣ⁺))/20, 50)))
    # 
    (; tophysicsunits(effrangepars)..., effrangepars...)
end

df.a_fm = ones(Nm)
df.r_fm = ones(Nm).*1im
df.a⁻¹ = ones(Nm).*1im
df.r = ones(Nm).*1im
df.N = ones(Nm).*1im


let Nappr = 150, Rexp = abs(imag(Eᵦˣ⁺))/20
    # 
    for i in 1:Nm
        𝒜 = df.model[i]
        @time effrangepars = # 12s
            effectiverangeexpansion(
                Δe->denominator_II(𝒜, Eᵦˣ⁺+Δe, δm0_val),
                Δe->k3b(Eᵦˣ⁺+Δe),
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


# uncertainties
df.a_fm⁺ᵟ = ones(Nm)
df.r_fm⁺ᵟ = ones(Nm).*1im
df.a⁻¹⁺ᵟ = ones(Nm).*1im
df.r⁺ᵟ = ones(Nm).*1im
df.N⁺ᵟ = ones(Nm).*1im
df.a_fm⁻ᵟ = ones(Nm)
df.r_fm⁻ᵟ = ones(Nm).*1im
df.a⁻¹⁻ᵟ = ones(Nm).*1im
df.r⁻ᵟ = ones(Nm).*1im
df.N⁻ᵟ = ones(Nm).*1im

let Nappr = 150, Rexp = abs(imag(Eᵦˣ⁺))/20
    # 
    for i in 1:Nm
        𝒜 = df.model[i]
        # 
        # dm0+δdm0
        @time effrangepars⁺ᵟ = # 12s
            effectiverangeexpansion(
                Δe->denominator_II(𝒜, Eᵦˣ⁺+Δe, δm0_val+δm0.err),
                Δe->k3b(Eᵦˣ⁺+Δe),
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
                Δe->denominator_II(𝒜, Eᵦˣ⁺+Δe, δm0_val-δm0.err),
                Δe->k3b(Eᵦˣ⁺+Δe),
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



ps((a,b,c)) = "$(a)^{+$(b)}_{-$(c)}"

# test
print(select(df[[4,5],:], [:modelnames, :a_fm, :r_fm, :N]))
# modelnames  a_fm      r_fm                  N
#─────────────────────────────────────────────────────────────────────
# A₀_DˣD      -7.55077  0.159869-0.0238077im  0.00347571-3.76763e-5im
# A₀′_DˣD     -7.55237  0.165517-0.0216824im   0.0034758-3.76788e-5im

# interference term
print(select(df[[3,4],:], [:modelnames, :a_fm, :r_fm, :N]))
# modelnames  a_fm      r_fm                   N
#───────────────────────────────────────────────────────────────────────
# A₀_π⁺D⁰D⁰   -7.31763  -0.779907+0.260606im   0.00702464+0.000589854im
# A₀_DˣD      -7.55077   0.159869-0.0238077im  0.00347571-3.76763e-5im

# higher threshold
print(select(df[[1,2],:], [:modelnames, :a_fm, :r_fm, :N]))
# modelnames  a_fm      r_fm                   N
#──────────────────────────────────────────────────────────────────────
# A₀_full      -5.95934   -4.36196+0.4608im    0.0103207+0.000671301im
# A₀_full_DˣD  -7.32118  -0.626768+0.159036im   0.010316+0.000554445im


function printuncertainty(nt, s, processing=identity)
    v = nt[s]
    v⁺ᵟ = nt[Symbol(s,"⁺ᵟ")]
    v⁻ᵟ = nt[Symbol(s,"⁻ᵟ")]
    Δv⁺ = v⁺ᵟ-v
    Δv⁻ = v-v⁻ᵟ
    # 
    abc = [v,Δv⁺,Δv⁻]
    return ps(processing(x) for x in abc)
end



Xc(r::Real,a⁻¹::Real) = 1/sqrt(1+2r*a⁻¹/X2DDpi.fm_times_mev)
Xc(r::Complex,a⁻¹::Complex) = Xc(real(r),real(a⁻¹))
# 
Z_polosa_90 = map(r->Xc(r, df.a⁻¹[1]*1e3), (-11.9,0) .+ df.r_fm[1])
Z_polosa_95 = map(r->Xc(r, df.a⁻¹[1]*1e3), (-16.9,0) .+ df.r_fm[1])
# 
Z_hanhart_90 = map(r->Xc(r, df.a⁻¹[1]*1e3), (-11.9,0) .+ df.r_fm[2])
Z_hanhart_95 = map(r->Xc(r, df.a⁻¹[1]*1e3), (-16.9,0) .+ df.r_fm[2])




writejson(joinpath("results","nominal","effectiverangecauchy.json"),
        Dict(
            :allmodelcomputation => [
                Symbol(df.modelnames[i]) => Dict{Symbol,Any}(
                    :a_fm => df.a_fm[i],
                    :r_fm => df.r_fm[i],
                    :inv_a => df.a⁻¹[i],
                    :r => df.r[i],
                    :N => df.N[i]
                ) for i in 1:size(df,1)],
# 
            :investigationsummary => Dict(
                :inv_scatt_length_MeV => round(df.a⁻¹[1]*1e3, digits=2),
                :effective_range_fm => round(df.r_fm[1], digits=2),
                :effective_range_disp_fm => round(df.r_fm[2], digits=2),
                :effective_range_high_fm => round(df.r_fm[1]-df.r_fm[2], digits=2),
            ),
# 
            :withuncertainties => Dict(
                :inv_scatt_length_MeV_Re => printuncertainty(df[1,:], :a⁻¹, x->round(real(1e3*x),digits=2)),
                :inv_scatt_length_MeV_Im => printuncertainty(df[1,:], :a⁻¹, x->round(imag(1e3*x),digits=2)),
                :effective_range_fm_Re => printuncertainty(df[1,:], :r_fm, x->round(real(x),digits=2)),
                :effective_range_fm_Im => printuncertainty(df[1,:], :r_fm, x->round(imag(x),digits=2)),
                :effective_range_disp_fm_Re => printuncertainty(df[2,:], :r_fm, x->round(real(x),digits=2)),
                :effective_range_high_fm_Im => printuncertainty(df[2,:], :r_fm, x->round(real(x),digits=2))
            ),
# 
            :compositeness => Dict(
                :Z_hanhart_90 => Z_hanhart_90,
                :Z_hanhart_95 => Z_hanhart_95,
                :Z_polosa_90 => Z_polosa_90,
                :Z_polosa_95 => Z_polosa_95
            )
    ))
#


struct ERA
    a⁻¹::Complex{Float64}
    r::Complex{Float64}
    N::Complex{Float64}
end
# 
struct EFRMatching
    𝒜::Amplitude
    era::ERA
end
EFRMatching(𝒜, nt) = EFRMatching(𝒜, ERA(nt.a⁻¹, nt.r, nt.N))

function X2DDpi.denominator_II(era::ERA,e)
    @unpack a⁻¹, r, N = era
    return N*(a⁻¹ + r/2*k3b(e)^2-1im*k3b(e))
end

using RecipesBase
@recipe function f(m::EFRMatching)
    Δev = range(-0.01, 0.0251, length=30)
    @series begin
        calv = [denominator_II(m.era, Eᵦˣ⁺+Δe) for Δe in Δev]
        label := "effective range"
        (Δev, real(calv)) 
    end
    calv = [denominator_II(m.𝒜, Eᵦˣ⁺+Δe, δm0_val) for Δe in Δev]
    (Δev, real(calv))
end

let
    plot(size=(200*5,150), layout=grid(1,5), yaxis=false, yticks=false)
    plot!(sp=1,EFRMatching(df.model[1], df[1,:]), lab=string(df.modelnames[1]))
    plot!(sp=2,EFRMatching(df.model[2], df[2,:]), lab=string(df.modelnames[2]))
    plot!(sp=3,EFRMatching(df.model[3], df[3,:]), lab=string(df.modelnames[3]))
    plot!(sp=4,EFRMatching(df.model[4], df[4,:]), lab=string(df.modelnames[4]))
    plot!(sp=5,EFRMatching(df.model[5], df[5,:]), lab=string(df.modelnames[5]))
end
savefig(joinpath("plots", "testmatchefr.pdf"))

# 
@unpack w_matching, rho_inf =
    readjson(joinpath("results", "nominal", "effective_range.json"))["effective_range_parameters"]["technical"]
# 
N′ = 1 / (w_matching * rho_inf)
let i = 1
    plot(xlab=L"\delta' m\,\,(\mathrm{MeV})", ylab=L"\mathrm{Re}\,\mathcal{A}^{-1}(black),\,\,\mathrm{Im}\,\mathcal{A}^{-1}(red)")
    plot!(Δe->N′*real(denominator_II(df.model[i], Δe, δm0_val)), -1, 3, lab="LHCb Model", lc=:black)
    plot!(Δe->N′*imag(denominator_II(df.model[i], Δe, δm0_val)), -1, 3, lab="", lc=:red)
    # 
    @unpack a⁻¹, r, N = df[i,:]
    era = ERA( a⁻¹, r, N)
    plot!(Δe->N′*real(denominator_II(era, Δe)), -1, 3, lab="Eff.-range exp.", ls=:dash, lc=:black)
    plot!(Δe->N′*imag(denominator_II(era, Δe)), -1, 3, lab="", ls=:dash, lc=:red)
    # plot!(Δe->imag(2*denominator_II(df.model[4], Eᵦˣ⁺+Δe, δm0_val)),-0.01, 0.1, lab=string(df.modelnames[4]))
    vline!([0], lab="", c=:green, lw=1)
    # scatter!([0.0], [N′*real(denominator_II(era, 0))], ms=5, mc=:black, lab="exp. point")
    # scatter!([0.0], [N′*imag(denominator_II(era, 0))], ms=5, mc=:red, lab="")
end
savefig(joinpath("plots","nominal", "effectiverangecauchy.pdf"))


