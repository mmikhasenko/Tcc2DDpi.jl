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
# Translate older energy grids at the cache boundary.
function interpolated(x::NamedTuple)
    channel = x.channel
    if startswith(channel.type, "γDD") && !hasproperty(channel, :μ12)
        channel = merge(channel, (; μ12 = 1.0, μ13 = -3.77))
    end
    channel = NamedTuple{Tuple(sort(collect(keys(channel))))}(channel)
    cutoff, knots = hasproperty(x, :m_cutoff) ?
        (x.m_cutoff, x.knots) : (e2m(x.cutoff), e2m.(x.knots))
    return interpolated(
        constructchannel(channel),
        cutoff,
        x.cutoffratio,
        interpolate((knots,), x.coefs, Gridded(Linear())),
    )
end
