struct DˣD{T1, T2} <: AbstractxDD
    ms::NamedTuple{(:m1, :m2, :m3), T1}
    R::T2
end

function decay_matrix_element_squared(d::DˣD, s, σ3, σ2)
    A = xDDPhaseSpace.λ(σ3, d.ms[1]^2, d.ms[2]^2) / (4 * σ3)
    return Jᴵ(σ3, d.R) * Jᴵᴵ(σ3, d.R) * xDDPhaseSpace.f² / 3 / 4 * A
end

function integrand_mapped_thr(d::DˣD, s, x)
    method = mapdalitzmethod(d)
    (σ3, σ2), jac = mapdalitz(method, (x, 0.0), d.ms, s)
    return decay_matrix_element_squared(d, s, σ3, σ2) / (2π * s) * jac
end

function ρ_thr(d::DˣD, m::Complex)
    integrand(x) = integrand_mapped_thr(d, m^2, x)
    v = quadgk(integrand, 0, 1)[1]
    v / (8π)^2
end

function ρ_thr(d::DˣD, m::Real)
    integrand(x) = real(integrand_mapped_thr(d, m^2, x))
    v = quadgk(integrand, 0, 1)[1]
    v[1] / (8π)^2
end

function ρ_tb(d::DˣD, m::Real)
    mR, mB = d.R.m, d.ms.m3
    sqrts = m
    sqrts < mR + mB ? 0.0 :
    sqrt(λ(m^2, mR^2, mB^2)) / m^2
end


