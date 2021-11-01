struct interpolated{T<:AbstractxDD,F1<:Real,F2<:Real,V}
    channel::T
    cutoff::F1
    cutoffratio::F2
    itr::V
end
function interpolated(channel::AbstractxDD, cutoff::Real; estep=0.01)
    eth = m2e(sum(channel.ms))
    ev = eth:estep:(cutoff+2*estep)
    calv = ρ_thr.(Ref(channel),ev)
    itr = interpolate((ev,), calv, Gridded(Linear()))
    cutoffratio = real(ρ_thr(channel,cutoff))/ρ_tb(channel,cutoff)
    interpolated(channel,cutoff,cutoffratio,itr)
end
# 
function ρ_thr(d::interpolated, e::Real)
    e < m2e(sum(d.channel.ms)) && return 0.0
    e < d.cutoff ?
        d.itr(e) :
    ρ_tb( d.channel, e) * d.cutoffratio
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
        constructchannel(x.channel),
        x.cutoff,
        x.cutoffratio,
        interpolate((x.knots,), x.coefs, Gridded(Linear()))
    )