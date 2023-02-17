using Pkg
cd(joinpath(@__DIR__,".."))
Pkg.activate(".")
# 
using X2DDpi
using Parameters
using QuadGK

using Plots
using LaTeXStrings
theme(:wong2, size=(500,350), minorticks=true, grid=false, frame=:box,
    guidefontvalign=:top, guidefonthalign=:right,
    foreground_color_legend = nothing,
    legendfontsize=9, legend =:topright,
    xlim=(:auto,:auto), ylim=(:auto,:auto),
    xlab=L"m(D^+D^0)\,[\mathrm{GeV}]")
# 

settings = transformdictrecursively!(readjson("settings.json"), ifstringgivemeasurement)
@unpack δm0 = settings["fitresults"]
const δm0_val = δm0.val
const Γ0 = 48e-3 # MeV
# 
const srange = (e2m.(δm0_val .+ (-1, 1) .* 2Γ0)) .^ 2
# 

# 
channels = [
    πDD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), BW(m=mDˣ⁺,Γ=ΓDˣ⁺), BW(m=mDˣ⁺,Γ=ΓDˣ⁺)),
    πDD((m1=mπ⁰,m2=mD⁺,m3=mD⁰), BW(m=mDˣ⁺,Γ=ΓDˣ⁺), BW(m=mDˣ⁰,Γ=ΓDˣ⁰)),
    γDD((m1=mγ, m2=mD⁺,m3=mD⁰), BW(m=mDˣ⁺,Γ=ΓDˣ⁺), BW(m=mDˣ⁰,Γ=ΓDˣ⁰))]
# 


data_2 = let 
    m = e2m(δm0_val)
    xv = range(3.734, 3.758, length=100)
    yv = map(e1->projectto1(channels[2], m^2, e1^2), xv)
    (; xv, yv, m)
end

data_3 = let 
    @unpack xv, m = data_2
    yv = map(e1->projectto1(channels[3], m^2, e1^2), xv)
    (; xv, yv, m)
end

writejson(joinpath("results","nominal","DpDz_spectrum_atthepeak.json"), Dict(
    :pi0D0Dp => data_2,
    :gD0Dp => data_3))
#


#  _|_|_|    _|          _|  
#  _|    _|  _|          _|  
#  _|_|_|    _|    _|    _|  
#  _|    _|    _|  _|  _|    
#  _|_|_|        _|  _|      


# BW model

# 
bw_e(s,e0,Γ) = J_I(s,
    BW_norm(m=e2m(e0), Γ=Γ*1e-3, n=sqrt(e2m(e0)*(Γ*1e-3) / π)))
# 
const bwn = quadgk(s->abs2(bw_e(s,δm0_val,Γ0)), e2m(δm0_val-2Γ0)^2, e2m(δm0_val+2Γ0)^2)[1]

# 
import X2DDpi: projectto1
projectto1_bw3(d, σ1) = projectto1(d, σ1,
    s->bw_e(s,δm0_val,Γ0),
    srange...)
# 


#      _|_|            _|  _|  
#    _|      _|    _|  _|  _|  
#  _|_|_|_|  _|    _|  _|  _|  
#    _|      _|    _|  _|  _|  
#    _|        _|_|_|  _|  _|  


# full model


modelDict = readjson(joinpath("results","nominal","model.json"))
const ichannels = interpolated.(d2nt.(modelDict["ichannels"]))
const ampX0 = Amplitude(Tuple(ichannels))

ampTcc(s) = 1e-3/denominator_I(ampX0, m2e(sqrt(s)), δm0_val)

let
    plot()
    #
    n = abs2(ampTcc(e2m(δm0_val)^2))
    plot!(e->abs2(ampTcc(e2m(e)^2))/n,-0.5, 1.5, yscale=:log10)
    #
    n = abs2(bw_e(e2m(δm0_val)^2,δm0_val,Γ0))
    plot!(e->abs2(bw_e(e2m(e)^2,δm0_val,Γ0))/n,-0.5, 1.5, yscale=:log10)
end

projectto1_Tcc(d, σ1) = projectto1(d, σ1,
    ampTcc,
    srange[1], (mDˣ⁰+mD⁺)^2)
# 
# takes 30 mins
data_2′ = let 
    @unpack xv = data_2
    yv = map(e1->projectto1_Tcc(channels[2], e1^2), xv)
    (; xv, yv)
end

# takes 30 mins
data_3′ = let 
    @unpack xv = data_2
    yv = map(e1->projectto1_Tcc(channels[3], e1^2), xv)
    (; xv, yv)
end

writejson(joinpath("results","nominal","DpDz_spectrum_withDxzDpcutoff.json"), Dict(
    :pi0D0Dp => data_2′,
    :gD0Dp => data_3′))
# 





let 
    plot()
    n23 = sum(data_2.yv+data_3.yv) / 1e3
    plot!(data_3.xv, data_3.yv ./ n23, lab=L"\gamma D^+D^0", fill=0, c=3, α=0.5)
    plot!(data_3.xv, data_3.yv ./ n23, lab="", l=(:black, 1.5))
    # 
    plot!(data_3.xv, (data_2.yv+data_3.yv) ./ n23, lab=L"\pi^0 D^+D^0", fill_between=data_3.yv./ n23, c=:red, α=0.3)
    plot!(data_3.xv, (data_2.yv+data_3.yv) ./ n23, lab="", l=(:black, 1.5))
end
# 
savefig(joinpath("plots","nominal","DpDz_spectrum.pdf"))

let 
    n23′ = sum(data_2′.yv+data_3′.yv) / 1e3
    plot!(data_3′.xv, data_3′.yv ./ n23′, lab="", l=(:black, 1.5), ls=:dash)
    plot!(data_3′.xv, (data_2′.yv+data_3′.yv)  ./ n23′, lab="", l=(:black, 1.5), ls=:dash)
end
# 
savefig(joinpath("plots","nominal","DpDz_spectrum_with_full.pdf"))


# let
#     plot()
#     n23 = sum(data_2.yv+data_3.yv) / 1e3
#     plot!(data_3.xv, data_3.yv ./ n23, lab=L"\gamma D^+D^0", fill=0, c=3, α=0.5)
#     plot!(data_3.xv, data_3.yv ./ n23, lab="", l=(:black, 1.5))
#     # 
#     plot!(data_3.xv, (data_2.yv+data_3.yv) ./ n23, lab=L"\pi^0 D^+D^0", fill_between=data_3.yv./ n23, c=:red, α=0.3)
#     plot!(data_3.xv, (data_2.yv+data_3.yv) ./ n23, lab="", l=(:black, 1.5))
# end

@time data_DpD0 = let 
    ch = channels
    Nbins = 100
    xv = range(3.734, 3.758, length=Nbins+1)
    m = e2m(δm0_val)
    # 
    y2v = map(e1->projectto1(ch[2], m^2, e1^2), xv)
    y3v = map(e1->projectto1(ch[3], m^2, e1^2), xv)
    # 
    (; xv, y2v, y3v)
end

dpdzspectrum(data_DpD0.xv, data_DpD0.y2v, data_DpD0.y3v)

