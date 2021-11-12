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
using Cuba
using SpecialFunctions


theme(:wong2, size=(500,350), minorticks=true, grid=false, frame=:box,
    guidefontvalign=:top, guidefonthalign=:right,
    foreground_color_legend = nothing,
    legendfontsize=9, legend =:topright,
    xlim=(:auto,:auto), ylim=(:auto,:auto), lab="")
# 


###################################################################

@userplot Dalitz
@recipe function f(hp::Dalitz; what2apply=real)
    model, = hp.args
    iσx := 2
    iσy := 3
    density = σs->what2apply(X2DDpi.decay_matrix_element_squared(model,e2m(δm0_val)^2,σs.σ3,σs.σ2))
    (ThreeBodyMasses(m0=e2m(δm0_val), model.ms...), density)
end

###################################################################

#            _|                  _|              _|                _|  _|  
#    _|_|_|  _|_|_|    _|  _|_|        _|_|_|  _|_|_|_|    _|_|_|  _|  _|  
#  _|        _|    _|  _|_|      _|  _|_|        _|      _|    _|  _|  _|  
#  _|        _|    _|  _|        _|      _|_|    _|      _|    _|  _|  _|  
#    _|_|_|  _|    _|  _|        _|  _|_|_|        _|_|    _|_|_|  _|  _|  


function CrystallBall(x,α,n,xbar,σ) 
    μ = (x-xbar)/σ

    C = n/α/(n-1)*exp(-α^2/2)
    D = sqrt(2π)*erf(α/sqrt(2))
    N = 1/(2C+D)/σ
    # 
    abs(μ) < α && return N*exp(-(x-xbar)^2/(2σ^2))
    A = (n/α)^n * exp(-α^2/2)
    B = n/α - α
    return N*A*(B+abs(x-xbar)/σ)^(-n)
end


import AlgebraPDF: AbstractPDF, pars, lims, func, updatevalueorflag
struct convCrystallBall{T} <: AbstractPDF{1}
    source::T
    σ::Float64
    α::Float64
    n::Float64
end

pars(d::convCrystallBall, isfree::Bool) = pars(d.source, isfree)
lims(d::convCrystallBall) = lims(d.source)
updatevalueorflag(d::convCrystallBall, s::Symbol, isfree::Bool, v=getproperty(pars(d),s)) = 
    convCrystallBall(
        AlgebraPDF.ispar(d.source,s) ? updatevalueorflag(d.source,s,isfree,v) : d.source,
        d.σ, d.α, d.n)

function func(d::convCrystallBall, x::NumberOrTuple; p=pars(d))
    σ = func(d.σ, x; p)
    g(z) = CrystallBall(z,d.α,d.n,0,d.σ)
    f(z) = func(d.source, z; p)
    return quadgk(y->f(x-y) * g(y), -7*σ, +7*σ)[1]
end

# # tests
# t = Normalized(FunctionWithParameters((x;p)->x>0,∅), (-1,1))
# tc = convGauss(t,207.6e-3)
# tb = convCrystallBall(t,207.6e-3,1.33,4.58)

# plot(t)
# plot!(tc)
# plot!(tb)


#                                              _|_|                                
#  _|_|_|    _|  _|_|    _|_|    _|_|_|      _|      _|    _|  _|_|_|      _|_|_|  
#  _|    _|  _|_|      _|_|_|_|  _|    _|  _|_|_|_|  _|    _|  _|    _|  _|        
#  _|    _|  _|        _|        _|    _|    _|      _|    _|  _|    _|  _|        
#  _|_|_|    _|          _|_|_|  _|_|_|      _|        _|_|_|  _|    _|    _|_|_|  
#  _|                            _|                                                
#  _|                            _|                                                

using ThreeBodyDecay

function project_cos23(d, s, cos23)
    !(-1 < cos23 < 1) && return 0.0
    
    σ1_0, σ1_e = (d.ms[2]+d.ms[3])^2, (√s-d.ms[1])^2
        
    M²of1(σ1) = real(
        X2DDpi.decay_matrix_element_squared(d,s,
            X2DDpi.σ3of1(σ1, cos23, d.ms^2, s),
            X2DDpi.σ2of1(σ1, cos23, d.ms^2, s))
        ) *
        sqrt(λ(σ1,d.ms[2]^2,d.ms[3]^2)*λ(s,σ1,d.ms[1]^2)) / σ1
    return quadgk(M²of1, σ1_0, σ1_e)[1] / (2π*s)
end

function projectto1(d, s, σ1)

    σ1_0, σ1_e = (d.ms[2]+d.ms[3])^2, (√s-d.ms[1])^2
    !(σ1_0 < σ1 < σ1_e) && return 0.0
    # 
	σ3_0, σ3_e = X2DDpi.σ3of1_pm(σ1, d.ms^2, s)
    
    M²of3(σ3) = real(
        X2DDpi.decay_matrix_element_squared(d,s,
            σ3,sum(d.ms^2)+s-σ3-σ1))
    return quadgk(M²of3, σ3_0, σ3_e)[1] / (2π*s)
end


function projectto3(d, s, σ3)

    σ3_0, σ3_e = (d.ms[1]+d.ms[2])^2, (√s-d.ms[3])^2
    !(σ3_0 < σ3 < σ3_e) && return 0.0
    # 
	σ2_0, σ2_e = X2DDpi.σ2of3_pm(σ3, d.ms^2, s)
    
    M²of2(σ2) = real(X2DDpi.decay_matrix_element_squared(d,s,σ3,σ2))
    return quadgk(M²of2, σ2_0, σ2_e)[1] / (2π*s)
end
# 
function bw_e(s,e0,Γ)
    m = e2m(e0)
    return sqrt(m*(Γ*1e-3) / π) / (m^2-s-1im*m*(Γ*1e-3))
end
# 
function projectto3(d, σ3, lineshape, smin, smax)
    σ3_0 = (d.ms[1]+d.ms[2])^2
    σ3 < σ3_0 && return 0.0
    # 
    M²of2(s, σ2) = real(X2DDpi.decay_matrix_element_squared(d,s,σ3,σ2)) * abs2(lineshape(s))
    function M²of2(x)
        s = smin + x[1]*(smax-smin)
        # 
        σ3_e = (√s-d.ms[3])^2
        σ3 > σ3_e && return 0.0
        # 
        σ2_0, σ2_e = X2DDpi.σ2of3_pm(σ3, d.ms^2, s)
        σ2 = σ2_0 + x[2]*(σ2_e-σ2_0)
        M²of2(s, σ2)*(σ2_e-σ2_0)*(smax-smin) / (2π*s)
    end 
    return cuhre((x,f)->(f[1]=M²of2(x)), 2, 1)[1][1]
end

# @time projectto3(π⁺D⁰D⁰, e2m(δm0_val)^2, 2.007^2)
# @time projectto3(π⁺D⁰D⁰, 2.007^2)

# let d = π⁺D⁰D⁰, m = m_fig3
#     plot(e3->projectto3(d, m^2, e3^2),
#         2.004, 2.011, lab="")
# # 
#     plot!(e3->projectto3(d, e3^2)*1.19,
#         2.004, 2.011, lab="")
# end


polyakovfactor(m) = max(0, tanh((mDˣ⁺ - m ) / 350e-6))


#        _|              _|                
#    _|_|_|    _|_|_|  _|_|_|_|    _|_|_|  
#  _|    _|  _|    _|    _|      _|    _|  
#  _|    _|  _|    _|    _|      _|    _|  
#    _|_|_|    _|_|_|      _|_|    _|_|_|  

data_fig3 = CSV.read(joinpath("quickcheck", "Fig_3.csv"), DataFrame)
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
const Γ0 = 48e-3 # keV

π⁺D⁰D⁰ = πDD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), (m=mDˣ⁺,Γ=ΓDˣ⁺), (m=mDˣ⁺,Γ=ΓDˣ⁺))
# 
projectto3(d, σ3) = projectto3(d, σ3,
    s->bw_e(s,δm0_val,Γ0),
    e2m(δm0_val-2Γ0)^2, e2m(δm0_val+2Γ0)^2)

# # text
# let d = π⁺D⁰D⁰, m = m_fig3
#     plot(e3->projectto3(d, m^2, e3^2),
#         (d.ms[1]+d.ms[2]), (m-d.ms[3]), lab="")
# end
# let d = π⁺D⁰D⁰, m = m_fig3
#     plot(polyakovfactor, (d.ms[1]+d.ms[2]), (m-d.ms[3]),
#         lc=:red, lab="Polyakov factor", leg=:topleft)
# end

m_fig3 = e2m(δm0_val) # δm0_val
lims_fig3 = (2.004, 2.0105)


function project_convolv(d, e; lims)
    proj(e3) = projectto3(d, e2m(e)^2, e3^2)
    pdf = Normalized(FunctionWithParameters((x;p)->proj(x); p=∅), lims)
    c = convGauss(pdf, 363e-6)
    # c = convCrystallBall(pdf, 207.6e-6,1.33,5.58)
    return c
end

c3 = project_convolv(π⁺D⁰D⁰, δm0_val; lims=lims_fig3)
# 
data_nominal = let c = c3
    intensity(m3) = func(c, m3)*polyakovfactor(m3)
    xv = range(lims_fig3..., length=100)
    yv = intensity.(xv)
    (; xv,yv)
end

# plot
let scale = 0.051
    plot(xlab=L"m(D^0\pi^+) [\mathrm{GeV}]")
    N_nominal = 10e3/sum(data_nominal.yv) * scale
    plot!(data_nominal.xv, data_nominal.yv*N_nominal, lab="",
        fill=0, c=:red, α=0.3)
    plot!(data_nominal.xv, data_nominal.yv*N_nominal, lab="", l=(:black,1))
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
import X2DDpi: decay_matrix_element_squared, AbstractxDD
struct πDD_noDˣ{T1} <: AbstractxDD
    ms::NamedTuple{(:m1,:m2,:m3),T1}
end
function decay_matrix_element_squared(d::πDD_noDˣ,s,σ3,σ2)
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

λλ3(ms,s,σ) = λ(σ,ms[1]^2,ms[2]^2)*λ(s,σ,ms[3]^2)
ρ3(ms,s,σ) = (λλ = λλ3(ms,s,σ); λλ<0 ? 0.0 : sqrt(λλ)/(2*σ))
ρ3(σ) = ρ3(π⁺D⁰D⁰.ms,m_fig3^2,σ)
# 
data_phsp = let
    pdf = Normalized(FunctionWithParameters((x;p)->ρ3(x^2); p=∅), lims_fig3)
    c = convGauss(pdf, 363e-6)
    intensity(m3) = c(m3)*polyakovfactor(m3)
    xv = range(lims_fig3..., length=100)
    yv = intensity.(xv)
    (; xv,yv)
end

# plot
let scale = 0.055
    plot(xlab=L"m(D^0\pi^+)\,\,[\mathrm{GeV}]", leg=:topleft)
    N_nominal = 10e3/sum(data_nominal.yv) * scale
    plot!(data_nominal.xv, data_nominal.yv  .* N_nominal,
        fill=0, α=0.3, c=2, lab=L"\mathrm{nominal}")
    plot!(data_nominal.xv, data_nominal.yv .* N_nominal, lab="", l=(:black,1))
    
    N_noDˣ = 10e3/sum(data_noDˣ.yv) * scale
    plot!(data_noDˣ.xv, data_noDˣ.yv .* N_noDˣ, 
        fill=0, α=0.3, c=3, lab=L"\mathrm{broad}\,\,D^*")
    plot!(data_noDˣ.xv, data_noDˣ.yv .* N_noDˣ, lab="", l=(:black,1))
    # 
    N_phsp = 10e3/sum(data_phsp.yv) * scale
    plot!(data_phsp.xv, data_phsp.yv .* N_phsp, 
        fill=0, α=0.3, c=4, lab=L"\mathrm{phase\,\,space}")
    plot!(data_phsp.xv, data_phsp.yv .* N_phsp, lab="", l=(:black,1))
    #
    let ed = h.edges[1]
        xv = (ed[1:end-1]+ed[2:end]) ./ 2
        scatter!(xv, h.weights, yerr=sqrt.(h.weights), lab="Data", c=:black, leg=:topleft)
    end
end
savefig(joinpath("plots", "Dpi_noDstar.pdf"))

using DelimitedFiles
writedlm(joinpath("results", "nominal", "Dpi_noDstar.txt"),
    [data_nominal.xv data_nominal.yv data_noDˣ.xv data_noDˣ.yv data_phsp.xv data_phsp.yv])
# 


#  _|_|_|      _|    _|_|_|      _|    
#  _|    _|  _|  _|  _|    _|  _|  _|  
#  _|    _|  _|  _|  _|    _|  _|  _|  
#  _|    _|  _|  _|  _|    _|  _|  _|  
#  _|_|_|      _|    _|_|_|      _|    

# new decay chain: no DD resonance
import X2DDpi: decay_matrix_element_squared, AbstractxDD
struct πDD_DD{T} <: AbstractxDD
    ms::NamedTuple{(:m1,:m2,:m3),NTuple{3,Float64}}
    R::T
end
function decay_matrix_element_squared(d::πDD_DD,s,σ3,σ2)
	msq = d.ms^2
	v = (;s,s12=σ3,s13=σ2,msq)
    σ = X2DDpi.s23(v)
# 	
	J23_I, J23_II = J_I(σ,d.R), J_II(σ,d.R)
    qsq = λ(s,σ,msq[1]^2) / (4s)
# 	
	frakM = qsq * J23_I * J23_II
    return frakM
end


BW_D⁰D⁰(m,g) = BW_Swave(m,g,mD⁰,mD⁰)
BW_mΓ_D⁰D⁰(m,Γ) = BW_Swave(m,sqrt(m*Γ/sqrt(λ(m^2,mD⁰^2,mD⁰^2))*m^2),mD⁰,mD⁰)

resonance_DD = BW_mΓ_D⁰D⁰(2mD⁰+0.5e-3, 1.5e-3) #BW_D⁰D⁰(, 0.2)
let X = resonance_DD
    plot(e->real(λ(m_fig3^2,(1e-3*e+2mD⁰)^2,mπ⁺^2) *
        J_I( (1e-3*e+2mD⁰)^2,X)*
        J_II((1e-3*e+2mD⁰)^2,X)), 0, (m_fig3-mπ⁺-2mD⁰)*1e3,
        lab=L"\left|\mathrm{BW}_{D^0D^0}\right|^2",
        xlab=L"m(D^0D^0)-m_{D^0}-m_{D^0}\,\,[\mathrm{MeV}]")
end
#

π⁺D⁰D⁰_DD = πDD_DD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), resonance_DD)

dalitz(π⁺D⁰D⁰_DD, color=:lajolla)

# convolve
c3_DD = project_convolv(π⁺D⁰D⁰_DD, δm0_val; lims=lims_fig3)

# calculate
data_DD = let c = c3_DD
    intensity(m3) = func(c, m3)*polyakovfactor(m3)
    xv = range(lims_fig3..., length=100)
    yv = intensity.(xv)
    (; xv,yv)
end
# 

# plot
p1 = let scale = 0.055
    plot(xlab=L"m(D^0\pi^+)\,\,[\mathrm{GeV}]", leg=:topleft)
    N_nominal = 10e3/sum(data_nominal.yv) * scale
    plot!(data_nominal.xv, data_nominal.yv  .* N_nominal,
        fill=0, α=0.3, c=2, lab=L"\mathrm{nominal}")
    plot!(data_nominal.xv, data_nominal.yv .* N_nominal, lab="", l=(:black,1))
    # 
    N_DD = 10e3/sum(data_DD.yv) * scale
    plot!(data_DD.xv, data_DD.yv .* N_DD, 
        fill=0, α=0.3, c=:red, lab=L"DD\,\,\mathrm{resonance}")
    plot!(data_DD.xv, data_DD.yv .* N_DD, lab="", l=(:red,1))
    # 
    N_phsp = 10e3/sum(data_phsp.yv) * scale
    plot!(data_phsp.xv, data_phsp.yv .* N_phsp, 
        fill=0, α=0.3, c=4, lab=L"\mathrm{phase\,\,space}")
    plot!(data_phsp.xv, data_phsp.yv .* N_phsp, lab="", l=(:black,1))
    #
    let ed = h.edges[1]
        xv = (ed[1:end-1]+ed[2:end]) ./ 2
        scatter!(xv, h.weights, yerr=sqrt.(h.weights), lab="Data", c=:black, leg=:topleft)
    end
end
savefig(joinpath("plots", "Dpi_DDresonance.pdf"))

p2 = let scale = 1e3
    plot(xlab=L"m(D^0D^0)\,\,[\mathrm{GeV}]", leg=:topright)
    # 
    xv = range(2mD⁰, 3.740, length=300) # e2m(δm0_val)-mπ⁺
    # 
    d = π⁺D⁰D⁰
    proj(e1) = projectto1(d, e2m(δm0_val)^2, e1^2)
    yv = proj.(xv)
    N_nominal = 1/sum(yv) * scale
    plot!(xv, yv  .* N_nominal,
        fill=0, α=0.3, c=2, lab=L"\mathrm{nominal}")
    plot!(xv, yv .* N_nominal, lab="", l=(:black,1))
# 
    d = π⁺D⁰D⁰_DD
    proj(e1) = projectto1(d, e2m(δm0_val)^2, e1^2)
    yv = proj.(xv)
    N_DD = 1/sum(yv) * scale
    plot!(xv, yv .* N_DD, 
        fill=0, α=0.3, c=:red, lab=L"DD\,\,\mathrm{resonance}")
    plot!(xv, yv .* N_DD, lab="", l=(:red,1))
end
savefig(joinpath("plots", "DD_DDresonance.pdf"))


plot(p1,p2, size=(800,300), bottom_margin=5Plots.PlotMeasures.mm)
savefig(joinpath("plots", "Dpi-DD_DDresonance.pdf"))

# LLH scan
function llh_DD_scan(ΓDD=0.8e-3; ΔmDD = 0.5e-3) 
    resonance_DD = BW_mΓ_D⁰D⁰(2mD⁰+ΔmDD, ΓDD) #BW_D⁰D⁰(, 0.2)
    π⁺D⁰D⁰_DD = πDD_DD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), resonance_DD)

    c = project_convolv(π⁺D⁰D⁰_DD, δm0_val; lims=lims_fig3)

    intensity(m3) = func(c, m3)*polyakovfactor(m3)
    # quadgk(intensity, lims_fig3...)[1]
    const_integral = sum(intensity.(range(lims_fig3..., length=100)))
    # 
    xv = data_fig3.mass .* 1e-3
    yv = log.(intensity.(xv) ./ const_integral)
    return -sum(yv .* data_fig3.weight)
end

let ΔmDDv=range(0.2e-3,0.5e-3, length=5)
    Γv = 0.95:0.05:2.5
    # 
    plot(xlab=L"\Gamma_{D^0D^0}\,\,[\mathrm{MeV}]", ylab=L"\Delta \mathrm{NLL}")
    profiles = []
    for ΔmDD in ΔmDDv
        p = map(Γ->llh_DD_scan(Γ*1e-3; ΔmDD), Γv)
        push!(profiles, p)
    end
    pmin = minimum(vcat(profiles...))
    for (ΔmDD,p) in zip(ΔmDDv,profiles)
        plot!(Γv, p .- pmin, lab="ΔmDD=$(1e3*ΔmDD) MeV")
    end
    hline!([1.0], lc=:red)
    plot!()
end
savefig(joinpath("plots", "nllscan_DDresonance.pdf"))



p2 = let scale = 1e3
    plot(xlab=L"\cos\,\theta_{D^0D^0}", leg=:top)
    # 
    xv = range(-1+1e-5, 1-1e-5, length=300) # e2m(δm0_val)-mπ⁺
    # 
    proj(z) = project_cos23(d, e2m(δm0_val)^2, z)
    # 
    d = π⁺D⁰D⁰
    yv = proj.(xv)
    N_nominal = 1/sum(yv) * scale
    plot!(xv, yv  .* N_nominal,
        fill=0, α=0.3, c=2, lab=L"\mathrm{nominal}")
    plot!(xv, yv .* N_nominal, lab="", l=(:black,1))
    # 
    d = π⁺D⁰D⁰_DD
    proj(z) = project_cos23(d, e2m(δm0_val)^2, z)
    # 
    yv = proj.(xv)
    N_DD = 1/sum(yv) * scale
    plot!(xv, yv .* N_DD, 
        fill=0, α=0.3, c=:red, lab=L"DD\,\,\mathrm{resonance}")
    plot!(xv, yv .* N_DD, lab="", l=(:black,1))
end
savefig(joinpath("plots", "costhetaDD_DDresonance.pdf"))


plot(
    dalitz(π⁺D⁰D⁰, color=:lajolla, title=L"\mathrm{nominal}"),
    dalitz(π⁺D⁰D⁰_DD, color=:lajolla, title=L"DD\,\,\mathrm{resonance}"),
    size=(800,400), bottom_margin=5Plots.PlotMeasures.mm)
savefig(joinpath("plots", "dalitz_DDresonance.pdf"))

