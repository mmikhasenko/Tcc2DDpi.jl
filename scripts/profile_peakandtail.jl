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

# 
ch⁺⁰⁰_⁺⁺ = πDD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), BW(m=mDˣ⁺,Γ=ΓDˣ⁺), BW(m=mDˣ⁺,Γ=ΓDˣ⁺))
ch⁺⁰⁰_⁺X = πDD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), BW(m=mDˣ⁺,Γ=ΓDˣ⁺), ZeroBW(m=mDˣ⁺,Γ=ΓDˣ⁺))
# 
ch⁰⁺⁰_⁺⁰ = πDD((m1=mπ⁰,m2=mD⁺,m3=mD⁰), BW(m=mDˣ⁺,Γ=ΓDˣ⁺), BW(m=mDˣ⁰,Γ=ΓDˣ⁰))
ch⁰⁺⁰_⁺X = πDD((m1=mπ⁰,m2=mD⁺,m3=mD⁰), BW(m=mDˣ⁺,Γ=ΓDˣ⁺), ZeroBW(m=mDˣ⁰,Γ=ΓDˣ⁰))
#
chγ⁺⁰_⁺⁰ = γDD((m1=mγ, m2=mD⁺,m3=mD⁰), BW(m=mDˣ⁺,Γ=ΓDˣ⁺), BW(m=mDˣ⁰,Γ=ΓDˣ⁰))
chγ⁺⁰_⁺X = γDD((m1=mγ, m2=mD⁺,m3=mD⁰), BW(m=mDˣ⁺,Γ=ΓDˣ⁺), ZeroBW(m=mDˣ⁰,Γ=ΓDˣ⁰))
# 

@time for qch in [:ch⁺⁰⁰_⁺⁺, :ch⁺⁰⁰_⁺X, :ch⁰⁺⁰_⁺⁰, :ch⁰⁺⁰_⁺X, :chγ⁺⁰_⁺⁰, :chγ⁺⁰_⁺X]
    eval(quote
       $(Symbol("i",qch)) = interpolated($(qch), cutoff; estep=estep)
    end)
end


# πDD + πDD + γDD: Dˣ⁺ + Dˣ⁺ + Dˣ⁰
const A₀_full = (ich⁺⁰⁰_⁺⁺, ich⁰⁺⁰_⁺⁰, ichγ⁺⁰_⁺⁰) |> Amplitude

# πDD + πDD + γDD: Dˣ⁺
const A₀_full_DˣD = (ich⁺⁰⁰_⁺⁺, ich⁰⁺⁰_⁺X, ichγ⁺⁰_⁺X) |> Amplitude

# πDD: Dˣ⁺ + Dˣ⁺
const A₀_π⁺D⁰D⁰ =  (ich⁺⁰⁰_⁺⁺,) |> Amplitude

# πDD: Dˣ⁺
const A₀_DˣD = (ich⁺⁰⁰_⁺X,) |> Amplitude

# πDD: Dˣ⁺
const A₀′_DˣD = (interpolated(
        DˣD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), BW(m=mDˣ⁺, Γ=ΓDˣ⁺)),
        cutoff; estep=estep),) |> Amplitude
# 
# 
modelnames = [:A₀_full, :A₀_full_DˣD, :A₀_π⁺D⁰D⁰, :A₀_DˣD, :A₀′_DˣD]
const Nm = length(modelnames)



struct TccLineshape
    xv::Vector{Float64}
    yv::Vector{Float64}
end

function lineshape(a::Amplitude, ev=-1:0.003:2.5)
    yv = broadcast(e->1/abs2(denominator_I(a, e, δm0_val)), ev)
    yv .*= broadcast(e->ρ_thr(a.ik[1],e), ev)
    return TccLineshape(ev, yv)
end


@recipe function f(lsh::TccLineshape)
    # 
    @unpack xv, yv = lsh
    # 
    xguide --> "e [MeV]"
    yguide --> "|A|²ρᵢ"
    # 
    (xv, yv ./ sum(yv))
end


lsv = TccLineshape[]
for modelname in modelnames
    model = eval(modelname)
    push!(lsv, lineshape(model))
end


function ftail(lsh::TccLineshape)
    @unpack xv, yv = lsh
    xyv = collect(zip(xv, yv))
    Ipeak = sum(filter(xy->xy[1]<0, xyv)) do xy
        xy[2]
    end
    Itot = sum(yv)
    return (;ftail = (Itot-Ipeak) / Itot, Itot, Ipeak)
end


let 
    plot()
    for (modelname, lsh) in zip(modelnames, lsv)
        ft = ftail(lsh)
        plot!(lsh, label="$(modelname): fₜ = $(round(100*ft.ftail; digits=1))%")
    end
    plot!()
end
plot!(yscale=:log10)
savefig(joinpath("plots", "profile_peakandtail.pdf"))