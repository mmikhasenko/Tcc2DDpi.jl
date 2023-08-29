struct Amplitude{T,F}
	ik::T
	corrections::F
end

Amplitude(ik) = Amplitude(ik, zero)

denominator_I(a::Amplitude{<:interpolated},e,δm) =
	-(dispersive(a.ik,e)-real(dispersive(a.ik,δm))) + a.corrections(e)

denominator_II(a::Amplitude{<:interpolated},e,δm) = imag(e) ≥ 0 ?
	denominator_I(a,e,δm) :
	denominator_I(a,e,δm)-2im*ρ_thr(a.ik,e)

# sum of several channels
denominator_I(a::Amplitude{T},e,δm) where {T <: Tuple{Vararg{interpolated}}} =
	sum(-(dispersive(ik,e)-real(dispersive(ik,δm))) for ik in a.ik) +
		a.corrections(e)
denominator_II(a::Amplitude{T},e,δm) where {T <: Tuple{Vararg{interpolated}}} =
	imag(e) ≥ 0 ?
		denominator_I(a,e,δm) :
		denominator_I(a,e,δm)-2im*sum(ρ_thr(ik,e) for ik in a.ik)
#

function pole_position(a::Amplitude{T},δm; init = [δm,0.1*δm]) where {T <: Tuple{Vararg{interpolated}}}
	abs2D(x) = abs2(1000*denominator_II(a, complex(x...),δm))
	fitres = optimize(abs2D, init)
	return NamedTuple{(:m_pole, :half_Γ_pole, :invD)}(
		[fitres.minimizer..., fitres.minimum])
end


struct NonRelBW
end
function denominator_I(d::NonRelBW, e, δm)
	ik(e) = 1im*sqrt(λ(e2m(e)^2,mDˣ⁺^2,mD⁰^2))/(2*e2m(e))
	ike = ik(e)
	ik0 = ik(δm+1e-6im)
	return -ike + real(ik0)
end
