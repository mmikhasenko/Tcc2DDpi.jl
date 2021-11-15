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


function pole_with_g1g2(g1,g2)
    channels = [
        πDD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), BW_norm(m=mDˣ⁺,Γ=ΓDˣ⁺,c=g1), BW_norm(m=mDˣ⁺,Γ=ΓDˣ⁺,c=g1)),
        πDD((m1=mπ⁰,m2=mD⁺,m3=mD⁰), BW_norm(m=mDˣ⁺,Γ=ΓDˣ⁺,c=g1), BW_norm(m=mDˣ⁰,Γ=ΓDˣ⁰,c=g2)),
        γDD((m1=mγ, m2=mD⁺,m3=mD⁰), BW_norm(m=mDˣ⁺,Γ=ΓDˣ⁺,c=g1), BW_norm(m=mDˣ⁰,Γ=ΓDˣ⁰,c=g2))]
    #
    ichannels = interpolated.(channels, cutoff; estep=estep) # cutoff
    #
    ampX0 = Amplitude(Tuple(ichannels), zero)

    pole_position(ampX0,δm0_val)
end

# test 
pole_nominal  = pole_with_g1g2(1,1)
pole_nominal′ = pole_with_g1g2(2,2)
@assert abs(1 - 
        pole_nominal.half_Γ_pole /
        pole_nominal′.half_Γ_pole) < 1e-2 
# limiting cases
pole_Dˣ⁺D⁰ = pole_with_g1g2(1,0)
pole_Dˣ⁰D⁺ = pole_with_g1g2(0,1)


ϕv = range(-π/2,π/2,length=11)
calv = [pole_with_g1g2(rotate((1,1),ϕ)...) for ϕ in ϕv]

let
    plot(ylim=(0,70),
        xlab=L"\phi:\,\,[g_1, g_2]^T = R(\phi)[1,-1]^T", ylab=L"\Gamma_{T_{cc}^+}\,\,[\mathrm{keV}]")
    plot!(ϕv, -2e3 .* getproperty.(calv, :half_Γ_pole), lab="")
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

