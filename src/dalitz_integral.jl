σ2of3_pm(σ,msq,s) = (msq[1]+msq[3]+
	(σ+msq[1]-msq[2])*(s-σ-msq[3])/(2*σ)) .+ [-1,1] .*
		sqrt(λ(σ,msq[1],msq[2])*λ(s,σ,msq[3]))/(2*σ)
#
σ3of1_pm(σ,msq,s) = (msq[2]+msq[1]+
	(σ+msq[2]-msq[3])*(s-σ-msq[1])/(2*σ)) .+ [-1,1] .*
		sqrt(λ(σ,msq[2],msq[3])*λ(s,σ,msq[1]))/(2*σ)
# 

σ3of1(σ1, cos23, msq, s) = msq[1]+msq[2]+
    (s-σ1-msq[1])*(σ1+msq[2]-msq[3])/(2σ1) + 
    cos23*sqrt(λ(s,σ1,msq[1])*λ(σ1,msq[2],msq[3])) / σ1
# 
σ2of1(σ1, cos23, msq, s) = msq[1]+msq[3]+
    (s-σ1-msq[1])*(σ1-msq[2]+msq[3])/(2σ1) - 
    cos23*sqrt(λ(s,σ1,msq[1])*λ(σ1,msq[2],msq[3])) / σ1
#
 
function integrand_mapped_thr(d::AbstractxDD,s,x)
	# 	
	σ3_0, σ3_e = (d.ms[1]+d.ms[2])^2, (√s-d.ms[3])^2
	σ3 = σ3_0 + x[1]*(σ3_e-σ3_0) # straight path
	# 		
	σ2_0, σ2_e = σ2of3_pm(σ3, d.ms^2, s)
	σ2 = σ2_0 + x[2]*(σ2_e-σ2_0) # straight path
	#
	jac = (σ3_e-σ3_0)*(σ2_e-σ2_0)
	decay_matrix_element_squared(d,s,σ3,σ2) / (2π*s) * jac
end

function ρ_thr(d::AbstractxDD, e::Complex)
	integrand(x,f) = f[1:2] .= reim(integrand_mapped_thr(d,e2m(e)^2,x))
	v = cuhre(integrand, 2, 2)[1]
	complex(v...) / (8π)^2
end

function ρ_thr(d::AbstractxDD, e::Real)
	integrand(x,f) = f[1] = real(integrand_mapped_thr(d,e2m(e)^2,x))
	v = cuhre(integrand, 2, 1)[1]
	v[1] / (8π)^2
end
