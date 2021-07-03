σ2of3_pm(σ,msq,s) = (msq[1]+msq[3]+
	(σ+msq[1]-msq[2])*(s-σ-msq[3])/(2*σ)) .+ [-1,1] .*
		sqrt(λ(σ,msq[1],msq[2])*λ(s,σ,msq[3]))/(2*σ)

function integrand_mapped_thr(d::AbstractxDD,s,x)
	# 	
	σ3_0, σ3_e = (d.ms[1]+d.ms[2])^2, (√s-d.ms[3])^2
	σ3 = σ3_0 + x[1]*(σ3_e-σ3_0) # straight path
	# 		
	σ2_0, σ2_e = σ2of3_pm(σ3, d.ms^2, s)
	σ2 = σ2_0 + x[2]*(σ2_e-σ2_0) # straight path
	#
	jac = (σ3_e-σ3_0)*(σ2_e-σ2_0)
	decay_matrix_element(d,s,σ3,σ2) / (2π*s) * jac
end

function ρ_thr(d::AbstractxDD, e)
	integrand(x,f) = f[1:2] .= reim(integrand_mapped_thr(d,e2m(e)^2,x))
	v = cuhre(integrand, 2, 2)[1]
	complex(v...) / (8π)^2
end

