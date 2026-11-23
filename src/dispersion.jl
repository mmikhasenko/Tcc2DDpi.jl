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
