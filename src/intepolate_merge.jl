struct interpolated{T <: AbstractxDD, F1 <: Real, F2 <: Real, V}
    channel::T
    m_cutoff::F1
    cutoffratio::F2
    itr::V
end
"""
    interpolated(channel::AbstractxDD, m_cutoff::Real; mstep=1e-5)

Tabulate phase space on an absolute-mass grid. Both `m_cutoff` and `mstep`
are in GeV, as in xDDPhaseSpace. Below the channel's three-body threshold
return zero; above the cutoff use the matched two-body phase space.
Complex masses are evaluated directly with the underlying channel.
"""
function interpolated(channel::AbstractxDD, m_cutoff::Real; mstep = 1e-5)
    mth = sum(masses(channel))
    isfinite(mstep) && mstep > 0 || throw(ArgumentError("mstep must be positive and finite (GeV)"))
    isfinite(m_cutoff) && m_cutoff > mth ||
        throw(ArgumentError("m_cutoff must exceed the channel threshold (GeV)"))
    proxy = ρ_tb(channel, m_cutoff)
    isfinite(proxy) && proxy > 0 ||
        throw(ArgumentError("m_cutoff must lie above the two-body threshold"))
    mv = mth:mstep:(m_cutoff+2*mstep)
    # Phase space vanishes exactly at threshold. Avoid roundoff in the
    # degenerate Dalitz mapping there, particularly for non-nominal masses.
    calv = [0.0; ρ_thr.(Ref(channel), mv[2:end])]
    itr = interpolate((mv,), calv, Gridded(Linear()))
    cutoffratio = real(ρ_thr(channel, m_cutoff)) / proxy
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
    channel = constructchannel(x.channel)
    cutoff, knots = hasproperty(x, :m_cutoff) ?
        (x.m_cutoff, x.knots) : (e2m(x.cutoff), e2m.(x.knots))
    return interpolated(
        channel,
        cutoff,
        x.cutoffratio,
        interpolate((knots,), x.coefs, Gridded(Linear())),
    )
end
