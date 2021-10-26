using Pkg
cd(joinpath(@__DIR__,".."))
Pkg.activate(".")
# 
using X2DDpi
using Parameters
using Measurements
using Interpolations
using Statistics
using Plots
using CSV
using DataFrames
using Plots.StatsBase
using QuadGK
using ThreeBodyDecay
using AlgebraPDF
using LaTeXStrings


theme(:wong2, size=(500,350), minorticks=true, grid=false, frame=:box,
    guidefontvalign=:top, guidefonthalign=:right,
    foreground_color_legend = nothing,
    legendfontsize=9, legend =:topright,
    xlim=(:auto,:auto), ylim=(:auto,:auto))
# 


function projectto3(d, s, σ3)

    σ3_0, σ3_e = (d.ms[1]+d.ms[2])^2, (√s-d.ms[3])^2
    !(σ3_0 < σ3 < σ3_e) && return 0.0
    # 
	σ2_0, σ2_e = X2DDpi.σ2of3_pm(σ3, d.ms^2, s)
    
    M²of2(σ2) = abs2(X2DDpi.decay_matrix_element(d,s,σ3,σ2))
    return quadgk(M²of2, σ2_0, σ2_e)[1] / (2π*s)
end
# 

polyakovfactor(m) = max(0, tanh((mDˣ⁺ - m ) / 350e-6))


#        _|              _|                
#    _|_|_|    _|_|_|  _|_|_|_|    _|_|_|  
#  _|    _|  _|    _|    _|      _|    _|  
#  _|    _|  _|    _|    _|      _|    _|  
#    _|_|_|    _|_|_|      _|_|    _|_|_|  

data_fig3 = CSV.read(joinpath("quickcheck", "Fig3.txt"), DataFrame)
rename!(data_fig3, map(s->filter(x -> !isspace(x), s), names(data_fig3)))


#                            _|            
#    _|_|_|    _|_|      _|_|_|    _|_|    
#  _|        _|    _|  _|    _|  _|_|_|_|  
#  _|        _|    _|  _|    _|  _|        
#    _|_|_|    _|_|      _|_|_|    _|_|_|  


settings = transformdictrecursively!(readjson("settings.json"), ifstringgivemeasurement)
#
@unpack δm0 = settings["fitresults"]
const δm0_val = δm0.val

π⁺D⁰D⁰ = πDD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), (m=mDˣ⁺,Γ=ΓDˣ⁺), (m=mDˣ⁺,Γ=ΓDˣ⁺))
# 
m_fig3 = e2m(δm0_val) # δm0_val

# # text
# let d = π⁺D⁰D⁰, m = m_fig3
#     plot(e3->projectto3(d, m^2, e3^2),
#         (d.ms[1]+d.ms[2]), (m-d.ms[3]), lab="")
# end
# let d = π⁺D⁰D⁰, m = m_fig3
#     plot(polyakovfactor, (d.ms[1]+d.ms[2]), (m-d.ms[3]),
#         lc=:red, lab="Polyakov factor", leg=:topleft)
# end

lims0 = let d = π⁺D⁰D⁰, m = m_fig3
    (d.ms[1]+d.ms[2]), (m-d.ms[3])
end
lims_fig3 = (2.004, 2.0105)



function project_convolv(d, e; lims)
    proj(e3) = projectto3(d, e2m(e)^2, e3^2)
    pdf = Normalized(FunctionWithParameters((x;p)->proj(x); p=∅), lims)
    c = convGauss(pdf, 363e-6)
    return c
end

c3      = project_convolv(π⁺D⁰D⁰,      δm0_val; lims=lims_fig3)
# 
data_nominal = let c = c3
    intensity(m3) = func(c, m3)*polyakovfactor(m3)
    xv = range(lims_fig3..., length=100)
    yv = intensity.(xv)
    (; xv,yv)
end

# plot
let scale = 0.0027
    plot(xlab=L"m(D^0\pi^+) [\mathrm{GeV}]")
    plot!(data_nominal.xv, data_nominal.yv*scale, lab="",
        fill=0, c=:red, α=0.3)
    plot!(data_nominal.xv, data_nominal.yv*scale, lab="", l=(:black,1))
end
# overlap with the data
h = Plots.StatsBase.fit(Histogram, data_fig3.mass .* 1e-3,
    weights(data_fig3.weight), range(lims_fig3..., length=30))
let ed = h.edges[1]
    xv = (ed[1:end-1]+ed[2:end]) ./ 2
    scatter!(xv, h.weights, yerr=sqrt.(h.weights), lab="Data", c=:black, leg=:topleft)
end

                                                                                                     
#                                                                    _|                                
#    _|_|_|    _|_|    _|_|_|  _|_|    _|_|_|      _|_|_|  _|  _|_|        _|_|_|    _|_|    _|_|_|    
#  _|        _|    _|  _|    _|    _|  _|    _|  _|    _|  _|_|      _|  _|_|      _|    _|  _|    _|  
#  _|        _|    _|  _|    _|    _|  _|    _|  _|    _|  _|        _|      _|_|  _|    _|  _|    _|  
#    _|_|_|    _|_|    _|    _|    _|  _|_|_|      _|_|_|  _|        _|  _|_|_|      _|_|    _|    _|  
#                                      _|                                                              
#                                      _|                                                              


# new decay chain: no D*
import X2DDpi: decay_matrix_element, AbstractxDD
struct πDD_noDˣ{T1} <: AbstractxDD
    ms::NamedTuple{(:m1,:m2,:m3),T1}
end
function decay_matrix_element(d::πDD_noDˣ,s,σ3,σ2)
	msq = d.ms^2
	v = (;s,s12=σ3,s13=σ2,msq)
# 	
	J12_I, J12_II = 1.0, 1.0
	J13_I, J13_II = 1.0, 1.0
# 	
	frakM = X2DDpi.A(v) * J12_I * J12_II +
			X2DDpi.B(v) * J13_I * J13_II +
			X2DDpi.C(v) * (J13_I * J12_II +  J12_I * J13_II)
    X2DDpi.f²*frakM/3/4
end

# create
π⁺D⁰D⁰_noDˣ = πDD_noDˣ((m1=mπ⁺,m2=mD⁰,m3=mD⁰))
# convolve
c3_noDˣ = project_convolv(π⁺D⁰D⁰_noDˣ, δm0_val; lims=lims_fig3)

# calculate
data_noDˣ = let c = c3_noDˣ
    intensity(m3) = func(c, m3)*polyakovfactor(m3)
    xv = range(lims_fig3..., length=100)
    yv = intensity.(xv)
    (; xv,yv)
end
# 
# plot
let scale = 0.0027
    plot(xlab=L"m(D^0\pi^+)\,\,[\mathrm{GeV}]", leg=:topleft)
    plot!(data_nominal.xv, data_nominal.yv*scale./ sum(data_nominal.yv),
        fill=0, α=0.3, c=2, lab=L"\mathrm{nominal}")
    plot!(data_nominal.xv, data_nominal.yv*scale./ sum(data_nominal.yv), lab="", l=(:black,1))
    # 
    plot!(data_noDˣ.xv, data_noDˣ.yv*scale./ sum(data_noDˣ.yv), 
        fill=0, α=0.3, c=3, lab=L"\mathrm{no}\,\,D^*")
    plot!(data_noDˣ.xv, data_noDˣ.yv*scale./ sum(data_noDˣ.yv), lab="", l=(:black,1))
    #
end
