struct DˣD{T1,T2} <: AbstractxDD
    ms::NamedTuple{(:m1,:m2,:m3),T1}
    R::T2
end

function decay_matrix_element_squared(d::DˣD,s,σ3,σ2)
	A = λ(σ3,d.ms[1]^2,d.ms[2]^2)/(4*σ3)
	return J_I(σ3,d.R) * J_II(σ3,d.R) * f²/3/4 * A
end

function integrand_mapped_thr(d::DˣD,s,x)
	method = HookSqrtDalitzMapping()
	# method = LinearDalitzMapping()
	(σ3,σ2), jac = mapdalitz(method, (x,0.0), d.ms, s)
	return decay_matrix_element_squared(d,s,σ3,σ2) / (2π*s) * jac
end

function ρ_thr(d::DˣD, e::Complex)
	integrand(x) = integrand_mapped_thr(d,e2m(e)^2,x)
	v = quadgk(integrand, 0, 1)[1]
	complex(v...) / (8π)^2
end

function ρ_thr(d::DˣD, e::Real)
	integrand(x) = real(integrand_mapped_thr(d,e2m(e)^2,x))
	v = quadgk(integrand, 0, 1)[1]
	v[1] / (8π)^2
end


function ρ_tb(d::DˣD, e::Real)
	M,m = d.R.m, d.ms.m3
	sqrts = e2m(e)
	sqrts < M+m ? 0.0 :
    	sqrt(λ(e2m(e)^2,M^2,m^2))/e2m(e)^2
end	
