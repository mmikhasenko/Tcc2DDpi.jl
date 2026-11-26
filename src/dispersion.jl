"""
    dispersive(channel::interpolated, m)

Once-subtracted dispersion integral in s = m², with absolute mass m in GeV.
The subtraction is at s = 0. For real masses evaluate the upper-rim limit
using m + 1e-9im GeV (the former 1e-6im MeV regulator).
The channel threshold and spectral density determine the result; no Tcc
reference mass enters this calculation.
"""
dispersive(d::interpolated, m::Real) = dispersive(d, m + 1e-9im)
function dispersive(d::interpolated, m)
    s = m^2
    mth = sum(masses(d.channel))
    function integrand(s′)
        m′ = sqrt(s′)
        ρ_thr(d, m′) / s′ / (s′ - s)
    end
    s / π * quadgk(integrand, mth^2, Inf)[1]
end
