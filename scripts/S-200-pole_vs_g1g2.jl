using Pkg
cd(joinpath(@__DIR__,".."))
Pkg.activate(".")
# 
using X2DDpi
using Parameters

using Plots
using LaTeXStrings
theme(:wong2, size=(500,350), minorticks=true, grid=false, frame=:box,
    guidefontvalign=:top, guidefonthalign=:right,
    foreground_color_legend = nothing,
    legendfontsize=9, legend =:topright,
    xlim=(:auto,:auto), ylim=(:auto,:auto),
    xlab=L"m(D^0\pi^+)\,[\mathrm{GeV}]")
# 


settings = transformdictrecursively!(readjson("settings.json"), ifstringgivemeasurement)
@unpack cutoff, estep = settings["phspmatching"]
@unpack δm0 = settings["fitresults"]
const δm0_val = δm0.val
#

function rotate((x,y), ϕ)
    x′ = x*cos(ϕ) - y*sin(ϕ)
    y′ = x*sin(ϕ) + y*cos(ϕ)
    return (x′,y′)
end


function precalculate_channels_g1g2(g1,g2)
    channels = [
        πDD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), BW_norm(m=mDˣ⁺,Γ=ΓDˣ⁺,n=g1), BW_norm(m=mDˣ⁺,Γ=ΓDˣ⁺,n=g1)),
        πDD((m1=mπ⁰,m2=mD⁺,m3=mD⁰), BW_norm(m=mDˣ⁺,Γ=ΓDˣ⁺,n=g1), BW_norm(m=mDˣ⁰,Γ=ΓDˣ⁰,n=g2)),
        γDD((m1=mγ, m2=mD⁺,m3=mD⁰), BW_norm(m=mDˣ⁺,Γ=ΓDˣ⁺,n=g1), BW_norm(m=mDˣ⁰,Γ=ΓDˣ⁰,n=g2))]
    #
    ichannels = interpolated.(channels, cutoff; estep=estep)
    # 
    (; g1, g2, channels, ichannels)
end

ϕv = range(-π/2,π/2,length=11)
ichsv = [precalculate_channels_g1g2(rotate((1,1),ϕ)...) for ϕ in ϕv]

import X2DDpi: obj2nt
obj2nt(nt::NamedTuple{(:g1, :g2, :channels, :ichannels)}) = (;
    g1=nt.g1, g2=nt.g2,
    channels=obj2nt.(nt.channels), ichannels=obj2nt.(nt.ichannels))

writejson(joinpath("results","nominal","setofmodelg1g2.json"), Dict(
        :angles => ϕv,
        :sermodels => obj2nt.(ichsv))
    )
# 


#                                      _|  
#  _|  _|_|    _|_|      _|_|_|    _|_|_|  
#  _|_|      _|_|_|_|  _|    _|  _|    _|  
#  _|        _|        _|    _|  _|    _|  
#  _|          _|_|_|    _|_|_|    _|_|_|  



setofmodelg1g2 = readjson(joinpath("results","nominal","setofmodelg1g2.json"))
@unpack angles, sermodels = setofmodelg1g2


function dict2model(d::Dict)
    ichannels = interpolated.(d2nt.(d["ichannels"]))
    return Amplitude(Tuple(ichannels))
end

models = [dict2model(m) for m in sermodels]

# takes 20 secs
polev = [pole_position(m,δm0_val) for m in models]

let
    plot(ylim=(0,70),
        xlab=L"\phi:\,\,[g_1, g_2]^T = R(\phi)[1,-1]^T", ylab=L"\Gamma_{T_{cc}^+}\,\,[\mathrm{keV}]")
    plot!(ϕv, -2e3 .* getproperty.(polev, :half_Γ_pole), lab="")
    vline!([0.0], lab="", lc=:red,
        ann=(0,10, text(L"\mathrm{nominal}:\,\,\,\,I=1",10,rotation=90,:left,:bottom)))
    # 
    vline!([-π/4], lab="", lc=:orange,
        ann=(-π/4,10,text(L"g(D^{*+}D^+)=0",10,rotation=90,:left,:bottom)))
    vline!([π/4], lab="", lc=:orange,
        ann=( π/4,40,text(L"g(D^{*+}D^0)=0",10,rotation=90,:left,:bottom)))
end
# 
savefig(joinpath("plots","nominal","pole_vs_g1g2.pdf"))


ps1 = [let ch = getproperty.(m.ik, :channel)
    Nbins = 100
    xv = range(3.734, 3.758, length=Nbins+1)
    m = e2m(δm0_val)
    # 
    y2v = map(e1->projectto1(ch[2], m^2, e1^2), xv)
    y3v = map(e1->projectto1(ch[3], m^2, e1^2), xv)

    dpdzspectrum(xv, y2v, y3v)
end for m in models]

# ps2 = [let #m = models[6]
#     ch = getproperty.(m.ik, :channel)
#     Nbins = 100
#     xv = range(3.734, 3.758, length=Nbins+1)
#     m = e2m(δm0_val)
#     # 
#     y2v = map(e1->projectto1(ch[2], m^2, e1^2), xv)
#     y3v = map(e1->projectto1(ch[3], m^2, e1^2), xv)

#     dpdzspectrum(xv, y2v, y3v)
# end for m in models]

# DzDzpip
ps2 = [let #m = models[6]
    ich = m.ik
    ch = getproperty.(m.ik, :channel)
    # 
    Nbins = 200
    ev = range(-1, 2.5, length=Nbins+1)
    calv = map(e->1/abs2(denominator_I(m, e, δm0_val)), ev)/1e3
    phsps = map(e->ρ_thr(ich[1],e), ev)

    plot(xlab=L"e\,\,[MeV]", ylab=L"|A|^2 \rho")
    plot!(ev, calv .* phsps, lab="", fill=0, α=0.5, c=:red)
    plot!(ev, calv .* phsps, lab="", l=(:black, 1.5))
end for m in models]

# DzDzpip
ps3 = [let #m = models[6]
    ich = m.ik
    ch = getproperty.(m.ik, :channel)
    # 
    Nbins = 200
    xv = range(3.729, 3.740, length=Nbins+1)
    m = e2m(δm0_val)
    # 
    y1v = map(e1->projectto1(ch[1], m^2, e1^2), xv)
    # 
    plot(xlab=L"m(D^0D^0)\,\,[\mathrm{GeV}]")
    plot!(xv, y1v, lab="", fill=0, α=0.5, c=:red)
    plot!(xv, y1v, lab="", l=(:black, 1.5))
end for m in models]

# Dzpip
ps4 = [let #m = models[6]
    ich = m.ik
    ch = getproperty.(m.ik, :channel)
    # 
    Nbins = 200
    xv = range(2.004, 2.0105, length=Nbins+1)
    m = e2m(δm0_val)
    # 
    y1v = map(e1->projectto3(ch[1], m^2, e1^2), xv)
    # 
    plot(xlab=L"m(D^0\pi^+)\,\,[\mathrm{GeV}]")
    plot!(xv, y1v, lab="", fill=0, α=0.5, c=:red)
    plot!(xv, y1v, lab="", l=(:black, 1.5))
end for m in models]

let
    plot(permutedims([[plot() for i in 1:length(ps2)] ps2 ps3 ps1 ps4])...,
        layout=grid(11,5, widths=(0.08,0.23,0.23,0.23,0.23)),
        size=(1000,1000),
        xlab="", ylab="", 
        xaxis=nothing, yaxis=nothing, lab="", frame=:nothing)
    # 
    plot!(sp=1, title=L"g_1/g_2")
    plot!(sp=2, title=L"m(D^0D^0\pi^+)")
    plot!(sp=3, title=L"m(D^0D^0)")
    plot!(sp=4, title=L"m(D^0D^+)")
    plot!(sp=5, title=L"m(D^0\pi^+)")
    # 
    for i in 1:11
        sm = sermodels[i]
        @unpack g1, g2 = sm
        plot!(sp=1+5*(i-1),
            ann=(0.5,0.5,text("$(round(g1, digits=1))/$(round(-g2, digits=1))", 13)),
            xaxis=false, yaxis=false)
    end
    plot!()
end
savefig(joinpath("plots","nominal","summary_g1g2.pdf"))
