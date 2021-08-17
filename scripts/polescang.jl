using Pkg
cd(joinpath(@__DIR__,".."))
Pkg.activate(".")
# 
using X2DDpi
using Parameters
using Measurements
using CSV
using DataFrames
using ProgressLogging
# using Plots

settings = transformdictrecursively!(readjson("settings.json"), ifstringgivemeasurement)
#
@unpack δm0 = settings["fitresults"]

# retrieve the model
modelDict = readjson(joinpath("results","nominal","model.json"))
const ichannels = interpolated.(d2nt.(modelDict["ichannels"]))
const channels = getproperty.(ichannels, :channel)
# 

const ampX0 = Amplitude(Tuple(ichannels))
@time pole_position(ampX0, δm0.val)
# 


import X2DDpi:pole_position
pole_position(δm::Float64,g::Float64) = pole_position(
    Amplitude(Tuple(ichannels),
        e->(e2m(δm)^2-e2m(e)^2)/(g^2/2)),
    δm)
    #
gdm = CSV.File(joinpath("results", "scang.txt"))
nt_gdm = NamedTuple.(gdm)

extrapolate(x,x23,y23) = diff(y23)[1]/diff(x23)[1]*(x-x23[1])+y23[1]

push!(nt_gdm, (;nt_gdm[end]..., g=0.001,
    dm = extrapolate(0.001, [nt_gdm[end-2].g,  nt_gdm[end].g],
                            [nt_gdm[end-2].dm, nt_gdm[end].dm]) ))
#
plot(getproperty.(nt_gdm,:g), getproperty.(nt_gdm,:dm))
# 
pole_mddl = let 
    a = []
    @progress for row in nt_gdm
        push!(a, pole_position(row[2]/1000, row[1]))
    end
    collect(a)
end

pole_left = let 
    a = []
    @progress for row in nt_gdm
        push!(a, pole_position((row[2]-row[3])/1000, row[1]))
    end
    a
end
pole_rght = let 
    a = []
    @progress for row in nt_gdm
        push!(a, pole_position((row[2]+row[3])/1000, row[1]))
    end
    a
end

writejson(joinpath("results","nominal","scang.json"),
        Dict(
            :middle => pole_mddl,
            :middle_plus_error => pole_left,
            :middle_minus_error => pole_rght,
        )
    )
#
# let
#     plot()
#     plot!(getproperty.(pole_left,:m_pole),
#         getproperty.(pole_left,:half_Γ_pole))
#     #
#     plot!(getproperty.(pole_mddl,:m_pole),
#         getproperty.(pole_mddl,:half_Γ_pole))
#     #
#     plot!(getproperty.(pole_rght,:m_pole),
#         getproperty.(pole_rght,:half_Γ_pole))
# end