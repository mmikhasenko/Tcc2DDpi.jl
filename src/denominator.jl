struct Amplitude{T, F}
    ik::T
    corrections::F
end

Amplitude(ik) = Amplitude(ik, zero)

"""
    denominator_I(a::Amplitude, e, δm)

First-sheet inverse amplitude. `e` and the subtraction point `δm` are MeV
relative to the D⁰D*⁺ reference threshold. The correction receives `e` in
MeV; phase-space and dispersion calls receive absolute masses in GeV.
"""
function denominator_I(a::Amplitude{<:interpolated}, e, δm)
    m = e2m(e)
    m0 = e2m(δm)
    -(dispersive(a.ik, m) - real(dispersive(a.ik, m0))) + a.corrections(e)
end

"""
    denominator_II(a::Amplitude, e, δm)

Continue the inverse amplitude to sheet II below the real axis. Both
energy arguments retain the same MeV convention as `denominator_I`.
"""
function denominator_II(a::Amplitude{<:interpolated}, e, δm)
    imag(e) ≥ 0 ?
    denominator_I(a, e, δm) :
    denominator_I(a, e, δm) - 2im * ρ_thr(a.ik, e2m(e))
end

# sum of several channels
denominator_I(a::Amplitude{T}, e, δm) where {T <: Tuple{Vararg{interpolated}}} =
    sum(-(dispersive(ik, e2m(e)) - real(dispersive(ik, e2m(δm)))) for ik in a.ik) +
    a.corrections(e)
denominator_II(a::Amplitude{T}, e, δm) where {T <: Tuple{Vararg{interpolated}}} =
    imag(e) ≥ 0 ?
    denominator_I(a, e, δm) :
    denominator_I(a, e, δm) - 2im * sum(ρ_thr(ik, e2m(e)) for ik in a.ik)
#
function pole_position(a::Amplitude{T}, δm; init = [δm, 0.1 * δm]) where {T <: Tuple{Vararg{interpolated}}}
    abs2D(x) = abs2(1000 * denominator_II(a, complex(x...), δm))
    fitres = optimize(abs2D, init)
    return NamedTuple{(:m_pole, :half_Γ_pole, :invD)}(
        [fitres.minimizer..., fitres.minimum])
end


struct NonRelBW
end
function denominator_I(d::NonRelBW, e, δm)
    ik(e) = 1im * sqrt(λ(e2m(e)^2, mDˣ⁺^2, mD⁰^2)) / (2 * e2m(e))
    ike = ik(e)
    ik0 = ik(δm + 1e-6im)
    return -ike + real(ik0)
end
