struct DˣD{T1,T2} <: AbstractxDD
    ms::NamedTuple{(:m1,:m2,:m3),T1}
    R::T2
end

function decay_matrix_element_squared(d::DˣD,s,σ3,σ2)
	A = λ(σ3,d.ms[1]^2,d.ms[2]^2)/(4*σ3)
	return J_I(σ3,d.R) * J_II(σ3,d.R) * f²/3/4 * A
end

function integrand_mapped_thr(d::DˣD,s,x)
	# 	
	σ3_0, σ3_e = (d.ms[1]+d.ms[2])^2, (√s-d.ms[3])^2
    # 
    # cut straight
    σ3 = σ3_0 + x[1]*(σ3_e-σ3_0) # straight path
	jac = (σ3_e-σ3_0)
    # 
    # cut down
    # σ3_i = real((√s-d.ms[3])^2)+sign(imag(s))*1e-6im
    # σ3 = x[1] < 0.5 ? 
    #         σ3_0 + ( x[1]     / 0.5)*(σ3_i-σ3_0) : # straight path
	#         σ3_i + ((x[1]-0.5)/ 0.5)*(σ3_e-σ3_i) # straight path
	# #
	# jac = x[1] < 0.5 ?
    #         (σ3_e-σ3_0) / 0.5 :
    #         (σ3_e-σ3_i) / 0.5
    #
    σ2 = 0.0
	othervar = 
		sqrt((sqrt(s)-d.ms[3])^2-σ3) *
		sqrt((sqrt(s)+d.ms[3])^2-σ3) *
		sqrt(σ3-(d.ms[1]-d.ms[2])^2) *
		sqrt(σ3-(d.ms[1]+d.ms[2])^2) / σ3
	decay_matrix_element_squared(d,s,σ3,σ2) / (2π*s) * jac * othervar * 2
end

function ρ_thr(d::DˣD, e::Complex)
	integrand(x) = integrand_mapped_thr(d,e2m(e)^2,[x])
	v = quadgk(integrand, 0, 1)[1]
	complex(v...) / (8π)^2
end

function ρ_thr(d::DˣD, e::Real)
	integrand(x) = real(integrand_mapped_thr(d,e2m(e)^2,[x]))
	v = quadgk(integrand, 0, 1)[1]
	v[1] / (8π)^2
end

function ρ_tb(d::DˣD, e::Real)
	M,m = d.R.m, d.ms.m3
	sqrts = e2m(e)
	sqrts < M+m ? 0.0 :
    	sqrt(λ(e2m(e)^2,M^2,m^2))/e2m(e)^2
end	
