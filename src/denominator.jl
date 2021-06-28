denominator_I(d::interpolated,e,δm) = - (dispersive(d,e)-real(dispersive(d,δm)))

denominator_II(d::interpolated,e,δm) = imag(e) ≥ 0 ?
	denominator_I(d,e,δm) :
	denominator_I(d,e,δm)+2im*ρ_thr(d,e)

# sum of several channels
denominator_I(ds::Tuple,e,δm) = sum(denominator_I(d,e,δm) for d in ds)
denominator_II(ds::Tuple,e,δm) = sum(denominator_II(d,e,δm) for d in ds)

intensity(e; δm) =
	1/abs2(denominator_II(Tuple(ichannels),e,δm))*ρ_thr(ichannels[1],e)

function pole_position(d,δm; init = [δm,-0.1*δm])
	abs2D(x) = abs2(denominator_II(d, complex(x...),δm))
	fitres = optimize(abs2D, init)
	return NamedTuple{(:m_pole, :half_Γ_pole, :invD)}(
		[fitres.minimizer..., fitres.minimum])
end