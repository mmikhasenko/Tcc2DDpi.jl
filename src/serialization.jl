# Read fields by name: JSON object order is not part of the cache format.
# Keep historical Tcc cache conventions here, outside the mass-based kernels.
function constructlineshape(x::NamedTuple)
    hasproperty(x, :type) || return (m = x.m, Γ = x.Γ)
    name = last(split(x.type, '.'))
    types = (; BW, BW_norm, BW_Swave, ZeroBW)
    hasproperty(types, Symbol(name)) || throw(ArgumentError("Unsupported lineshape: $(x.type)"))
    T = getproperty(types, Symbol(name))
    return T((getproperty(x, key) for key in fieldnames(T))...)
end

"""
    constructchannel(serialized::NamedTuple)

Restore a πDD, γDD, or DˣD channel independently of serialized field order.
Masses and widths are in GeV. Historical γDD caches without couplings use
this analysis's defaults, μ₊ and μ₀.
"""
function constructchannel(x::NamedTuple)
    name = last(split(first(split(x.type, '{')), '.'))
    ms = (; m1 = x.ms.m1, m2 = x.ms.m2, m3 = x.ms.m3)
    if name == "DˣD"
        return DˣD(ms, constructlineshape(x.R))
    end
    R12, R13 = constructlineshape(x.R12), constructlineshape(x.R13)
    name == "πDD" && return πDD(ms, R12, R13)
    name == "γDD" && return γDD(ms, R12, R13,
        get(x, :μ12, μ₊), get(x, :μ13, μ₀))
    throw(ArgumentError("Unsupported channel: $(x.type)"))
end

obj2nt(d::DˣD) = (; type = "DˣD", ms = d.ms, R = obj2nt(d.R))
