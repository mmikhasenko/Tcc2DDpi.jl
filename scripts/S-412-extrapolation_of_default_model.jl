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
@unpack cutoff, estep = settings["phspmatching"]


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
        cutoff; estep=estep)
    @time iπDD3 = interpolated(
        ChannelWithIntegrationMethod(ch1, HookSqrtDalitzMapping{3}()),
        cutoff; estep=estep)
    # 
    @time iπDD2′ = interpolated(
        ChannelWithIntegrationMethod(ch2, HookSqrtDalitzMapping{2}()),
        cutoff; estep=estep)
    @time iπDD3′ = interpolated(
        ChannelWithIntegrationMethod(ch2, HookSqrtDalitzMapping{3}()),
        cutoff; estep=estep)
    # 
    @time iγDD2 = interpolated(
        ChannelWithIntegrationMethod(ch3, HookSqrtDalitzMapping{2}()),
        cutoff; estep=estep)
    @time iγDD3 = interpolated(
        ChannelWithIntegrationMethod(ch3, HookSqrtDalitzMapping{3}()),
        cutoff; estep=estep)
    # 
    Amplitude((iπDD2, iπDD3, iπDD2′, iπDD3′, iγDD2, iγDD3))
end


effective_range_intepolation_settings = ( 
    ϵ0 = abs(imag(Eᵦˣ⁺)),
    ϵf_grid = [0.02, 0.05, 0.1, 0.2, 0.5],
    circular_sum_n = 100,
    model_δm0 = δm0_val)

ere_scan = effective_range_scan(model0,
    effective_range_intepolation_settings)
# 
ssc = effective_range_extrapolation(ere_scan)

@show ssc.r_regression;
@show ssc.inva_regression * 1000;

writejson(
    joinpath("results",
    "nominal",
    "effective-range-extrapolation.json"),
    ssc)


let
    df = DataFrame(ere_scan)
    # 
    ϵfv = df.ϵf
    xv = sqrt.(ϵfv)
    yv = df.r
    regression = linear_regression(xv, yv)
    @unpack α, β = regression
    # 
    plot(layout=grid(1,2), size=(900,400),
        title = "Regression",
        ylab=["Re r" "Im r"] .* " [fm]", xlab="scan scale, √ϵf",
        bottom_margin = 3mm, left_margin = 4mm)
    # 
    plot!(sp=1, x-> real(β + α * x), range(xv[[1,end]]..., 30))
    scatter!(sp=1, xv, yv .|> real, m=(7, :d))
    #
    plot!(sp=2, x-> imag(β + α * x), range(xv[[1,end]]..., 30))
    scatter!(sp=2, xv, yv .|> imag, m=(7, :d))
    #
    plot!()
end

