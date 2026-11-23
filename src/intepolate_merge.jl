struct interpolated{T <: AbstractxDD, F1 <: Real, F2 <: Real, V}
    channel::T
    m_cutoff::F1
    cutoffratio::F2
    itr::V
end
function interpolated(channel::AbstractxDD, m_cutoff::Real; e_step = 0.01)
    mth = sum(masses(channel))
    m_step = e_step / 1e3
    mv = mth:m_step:(m_cutoff+2*m_step)
    calv = ρ_thr.(Ref(channel), mv)
    itr = interpolate((mv,), calv, Gridded(Linear()))
    cutoffratio = real(ρ_thr(channel, m_cutoff)) / ρ_tb(channel, m_cutoff)
    interpolated(channel, m_cutoff, cutoffratio, itr)
end
# 
function ρ_thr(d::interpolated, m::Real)
    m < sum(masses(d.channel)) && return 0.0
    m < d.m_cutoff ?
    d.itr(m) :
    ρ_tb(d.channel, m) * d.cutoffratio
end
ρ_thr(d::interpolated, m::Complex) = ρ_thr(d.channel, m)


# serialization
obj2nt(ich::interpolated) =
    (channel = obj2nt(ich.channel),
        m_cutoff = ich.m_cutoff,
        cutoffratio = ich.cutoffratio,
        knots = ich.itr.knots[1],
        coefs = ich.itr.coefs)
#
# deserialization
interpolated(x::NamedTuple{(:m_cutoff, :coefs, :channel, :cutoffratio, :knots), T} where T) =
    interpolated(
        constructchannel(x.channel),
        x.m_cutoff,
        x.cutoffratio,
        interpolate((x.knots,), x.coefs, Gridded(Linear())),
    )
