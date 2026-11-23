struct interpolated{T<:AbstractxDD,F1<:Real,F2<:Real,V}
    channel::T
    cutoff::F1
    cutoffratio::F2
    itr::V
end
function interpolated(channel::AbstractxDD, cutoff::Real; estep=0.01)
    eth = m2e(sum(masses(channel)))
    ev = eth:estep:(cutoff+2*estep)
    calv = ρ_thr.(Ref(channel),ev)
    itr = interpolate((ev,), calv, Gridded(Linear()))
    cutoffratio = real(ρ_thr(channel,cutoff))/ρ_tb(channel,cutoff)
    interpolated(channel,cutoff,cutoffratio,itr)
end
# 
function ρ_thr(d::interpolated, e::Real)
    e < m2e(sum(masses(d.channel))) && return 0.0
    e < d.cutoff ?
        d.itr(e) :
    ρ_tb(d.channel, e2m(e)) * d.cutoffratio
end
ρ_thr(d::interpolated, e::Complex) = ρ_thr(d.channel,e)


# serialization
obj2nt(ich::interpolated) =
    (channel = obj2nt(ich.channel),
    cutoff = ich.cutoff,
    cutoffratio = ich.cutoffratio,
    knots = ich.itr.knots[1],
    coefs = ich.itr.coefs)
#
# deserialization
interpolated(x::NamedTuple{(:coefs, :cutoff, :channel, :cutoffratio, :knots),T} where T) =
    interpolated(
        construct_cached_channel(x.channel),
        x.cutoff,
        x.cutoffratio,
        interpolate((x.knots,), x.coefs, Gridded(Linear()))
    )

# Older nominal caches predate explicit magnetic transition moments.
function construct_cached_channel(channel)
    if startswith(channel.type, "γDD") && !hasproperty(channel, :μ12)
        channel = NamedTuple{(:R12, :R13, :μ12, :μ13, :ms, :type)}(
            merge(channel, (; μ12 = 1.0, μ13 = -3.77)))
    end
    return constructchannel(channel)
end
