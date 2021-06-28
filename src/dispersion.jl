dispersive(d::interpolated,e::Real) = dispersive(d,e+1e-6im)
function dispersive(d::interpolated,e)
	s = e2m(e)^2
	eth = m2e(sum(d.channel.ms))
	function integrand(s′)
		e′ = m2e(sqrt(s′))
		ρ_thr(d,e′)/s′/(s′-s)
	end
	s/π * quadgk(integrand, e2m(eth)^2, Inf)[1]
end