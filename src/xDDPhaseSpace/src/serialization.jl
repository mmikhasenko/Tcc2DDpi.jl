obj2nt(nt::NamedTuple) = nt

# deserilization
function constructchannel(ser::NamedTuple{(:R12, :R13, :ms, :type), T} where T)
    # Try to find the type in the parent module first, then in current module
    if startswith(ser.type, "γDD")
        type = isdefined(Main, :X2DDpi) ? Main.X2DDpi.γDD : γDD
    elseif startswith(ser.type, "DˣD")
        type = isdefined(Main, :X2DDpi) ? Main.X2DDpi.DˣD : DˣD
    else
        type = eval(Meta.parse(ser.type))
    end
    ms = NamedTuple{(:m1, :m2, :m3)}(ser.ms)
    R12 = constructlineshape(ser.R12)
    R13 = constructlineshape(ser.R13)
    type(ms, R12, R13)
end

constructlineshape(R::NamedTuple{(:m, :Γ), T} where T) = R
# 
function constructlineshape(R::NamedTuple)
    type = eval(Meta.parse(R.type))
    keysnotype = filter(x -> x != :type, keys(R))
    args = NamedTuple{keysnotype}(R)
    type(args...)
end