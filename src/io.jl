function writejson(path, obj)
    open(path, "w") do io
        JSON.print(io, obj, 4)
    end
end

function readjson(path)
    f = read(path, String)
    return JSON.parse(f)
end


function writetoml(path, obj)
    open(path, "w") do io
        TOML.print(io, obj)
    end
end

function readtoml(path)
    f = open(path)
    return TOML.parse(f)
end



function transformdictrecursively!(d::Dict, apply)
    for (k,v) in d
        if v isa Dict
            d[k] = transformdictrecursively!(v, apply)
            continue
        end
        if v isa NamedTuple
            d[k] = NamedTuple{keys(v)}(apply.(Tuple(v)))
            continue
        end
        d[k] = apply.(v)
    end
    return d
end

function ifstringgivemeasurement(x)
    if (x isa String) && (findfirst("±", x) !== nothing)
        return eval(Meta.parse(x))
    end
    return x
end

function ifmeasurementgivestring(x)
    if x isa Measurement
        return string(x)
    end
    return x
end

d2nt(d; process=identity) = NamedTuple{Tuple(Symbol.(keys(d)))}(process.(Base.values(d)))
