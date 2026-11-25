struct interpolated{T<:AbstractxDD,F1<:Real,F2<:Real,V}
    channel::T
    cutoff::F1
    cutoffratio::F2
    itr::V
end
# xDDPhaseSpace channels use absolute mass. The local single-channel model
# still uses the energy offset, just like the interpolation grid.
phase_space_argument(::AbstractxDD, e) = e2m(e)
phase_space_argument(::DˣD, e) = e

function interpolated(channel::AbstractxDD, cutoff::Real; estep=0.01)
    eth = m2e(sum(masses(channel)))
    ev = eth:estep:(cutoff+2*estep)
    calv = map(e -> ρ_thr(channel, phase_space_argument(channel, e)), ev)
    itr = interpolate((ev,), calv, Gridded(Linear()))
    arg = phase_space_argument(channel, cutoff)
    cutoffratio = real(ρ_thr(channel, arg)) / ρ_tb(channel, arg)
    interpolated(channel,cutoff,cutoffratio,itr)
end
# 
function ρ_thr(d::interpolated, e::Real)
    e < m2e(sum(masses(d.channel))) && return 0.0
    e < d.cutoff ?
        d.itr(e) :
    ρ_tb(d.channel, phase_space_argument(d.channel, e)) * d.cutoffratio
end
ρ_thr(d::interpolated, e::Complex) = ρ_thr(d.channel, phase_space_argument(d.channel, e))


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
