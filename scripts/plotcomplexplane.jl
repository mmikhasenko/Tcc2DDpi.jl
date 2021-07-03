using Parameters
using Interpolations

import Plots.PlotMeasures.mm
using Plots
using LaTeXStrings

using X2DDpi

# read 

settings = transformdictrecursively!(
    readjson("settings.json"),
    ifstringgivemeasurement)
@unpack δm0 = settings["fitresults"]

# precalculated
@unpack mgrid, ppvalues = readjson(joinpath("results","nominal","pole_interpolation.json"))["pole_interpolation"]
@unpack pole_sv = transformdictrecursively!(
    readjson(joinpath("results","nominal","pole.json"))["pole_position"],
    ifstringgivemeasurement)

# inferred
ppsampled = d2nt.(ppvalues)
itr_m, itr_Γ =
	interpolate((mgrid,), getproperty.(ppsampled, :m_pole), Gridded(Linear())),
	interpolate((mgrid,), 2 .* getproperty.(ppsampled, :half_Γ_pole), Gridded(Linear()))
pole_sv = NamedTuple{(:m_pole, :Γ_pole)}([itr_m(δm0), itr_Γ(δm0)])

# used channels to get branch points
channels = [
    πDD((m1=mπ⁺,m2=mD⁰,m3=mD⁰), (m=mDˣ⁺,Γ=ΓDˣ⁺), (m=mDˣ⁺,Γ=ΓDˣ⁺)),
    πDD((m1=mπ⁰,m2=mD⁺,m3=mD⁰), (m=mDˣ⁺,Γ=ΓDˣ⁺), (m=mDˣ⁰,Γ=ΓDˣ⁰)),
    γDD((m1=mγ, m2=mD⁺,m3=mD⁰), (m=mDˣ⁺,Γ=ΓDˣ⁺), (m=mDˣ⁰,Γ=ΓDˣ⁰))]

#
let
    plot(grid=false)
	plot!(
        getproperty.(ppsampled, :m_pole),
		getproperty.(ppsampled, :half_Γ_pole), l=(1,:gray), lab="")
	massv = range(δm0.val-δm0.err, δm0.val+δm0.err, length=15)
	plot!(itr_m.(massv), itr_Γ.(massv)/2, l=(2,:black), lab="")
	scatter!([itr_m(δm0.val)], [itr_Γ(δm0.val)/2], m=(5,:black), lab="")
    #
	vline!([δm0.val], lab="", lc=2, frame=:origin)
    # 
	for bp in vcat(collect.(X2DDpi.branch_points.(channels))...)
		plot!([bp, complex(1.8,imag(bp))], lc=:red, lab="")
		scatter!([bp], m=(:o,3,:red), lab="")
	end
    #	
	plot!(xlab=L"\Re\,(\delta' m)\,\,[\mathrm{MeV}]",
		  ylab=L"\Im\,(\delta' m)\,\,[\mathrm{MeV}]")
end
savefig(joinpath("plots","nominal","complex_plane.pdf"))

# latex
pgfplotsx()
savefig(joinpath("plots","latex","complex_plane.tex"))