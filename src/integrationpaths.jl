function σ2of3_pm(σ3,msq,s)
	regularpart = 
		(msq[1]+msq[3]+
		(σ3+msq[1]-msq[2])*(s-σ3-msq[3])/(2*σ3))
	# 
	sqrt_σ3 = sqrt(σ3)
	ms = sqrt.(msq)
	sqrt_s = sqrt(s)
	# 
	irregularpart = # sqrt(λλ)/σ
		sqrt((sqrt_s-ms[3])-sqrt_σ3) *  # upper end
		sqrt((sqrt_s-ms[3])+sqrt_σ3) *
		sqrt((sqrt_s+ms[3])-sqrt_σ3) *
		sqrt((sqrt_s+ms[3])+sqrt_σ3) *
		sqrt(sqrt_σ3-(ms[1]-ms[2])) *
		sqrt(sqrt_σ3+(ms[1]-ms[2])) *
		sqrt(sqrt_σ3-(ms[1]+ms[2])) * # lower end
		sqrt(sqrt_σ3+(ms[1]+ms[2])) / σ3
	# 
	return regularpart .+ [-1,1] .* (irregularpart / 2)
end
# 
function σ3of1_pm(σ1,msq,s)
	regularpart = 
		(msq[1]+msq[2]+
		(σ1+msq[2]-msq[3])*(s-σ1-msq[1])/(2*σ1))
	# 
	sqrt_σ1 = sqrt(σ1)
	ms = sqrt.(msq)
	sqrt_s = sqrt(s)
	# 
	irregularpart = # sqrt(λλ)/σ
		sqrt((sqrt_s-ms[1])-sqrt_σ1) *  # upper end
		sqrt((sqrt_s-ms[1])+sqrt_σ1) *
		sqrt((sqrt_s+ms[1])-sqrt_σ1) *
		sqrt((sqrt_s+ms[1])+sqrt_σ1) *
		sqrt(sqrt_σ1-(ms[2]-ms[3])) *
		sqrt(sqrt_σ1+(ms[2]-ms[3])) *
		sqrt(sqrt_σ1-(ms[2]+ms[3])) * # lower end
		sqrt(sqrt_σ1+(ms[2]+ms[3])) / σ1
	# 
	return regularpart .+ [-1,1] .* (irregularpart / 2)
end
# 
function σ3of2_pm(σ2,msq,s)
	regularpart = 
		(msq[1]+msq[2]+
		(σ2+msq[1]-msq[3])*(s-σ2-msq[2])/(2*σ2))
	# 
	sqrt_σ2 = sqrt(σ2)
	ms = sqrt.(msq)
	sqrt_s = sqrt(s)
	# 
	irregularpart = # sqrt(λλ)/σ
		sqrt((sqrt_s-ms[2])-sqrt_σ2) *  # upper end
		sqrt((sqrt_s-ms[2])+sqrt_σ2) *
		sqrt((sqrt_s+ms[2])-sqrt_σ2) *
		sqrt((sqrt_s+ms[2])+sqrt_σ2) *
		sqrt(sqrt_σ2-(ms[3]-ms[1])) *
		sqrt(sqrt_σ2+(ms[3]-ms[1])) *
		sqrt(sqrt_σ2-(ms[3]+ms[1])) * # lower end
		sqrt(sqrt_σ2+(ms[3]+ms[1])) / σ2
	# 
	return regularpart .+ [-1,1] .* (irregularpart / 2)
end
# 
σ3of1(σ1, cos23, msq, s) = msq[1]+msq[2]+
    (s-σ1-msq[1])*(σ1+msq[2]-msq[3])/(2σ1) + 
    cos23*sqrt(λ(s,σ1,msq[1])*λ(σ1,msq[2],msq[3])) / (2σ1)
# 
σ2of1(σ1, cos23, msq, s) = msq[1]+msq[3]+
    (s-σ1-msq[1])*(σ1-msq[2]+msq[3])/(2σ1) - 
    cos23*sqrt(λ(s,σ1,msq[1])*λ(σ1,msq[2],msq[3])) / (2σ1)
# 

abstract type AbstractDalitzMapping end 
struct LinearDalitzMapping <: AbstractDalitzMapping end
struct HookSqrtDalitzMapping{K} <: AbstractDalitzMapping end

function mapdalitz(method::LinearDalitzMapping, x, ms, s)
	# 	
	σ3_0, σ3_e = (ms[1]+ms[2])^2, (√s-ms[3])^2
	σ3 = σ3_0 + x[1]*(σ3_e-σ3_0) # straight path
	#
	jacobian = (σ3_e-σ3_0)
	# 
	σ2_0, σ2_e = σ2of3_pm(σ3, ms^2, s)
	σ2 = σ2_0 + x[2]*(σ2_e-σ2_0) # straight path
	#
	jacobian *= (σ2_e-σ2_0)
	# 
	return (σ3,σ2), jacobian
end

function mapdalitz(method::HookSqrtDalitzMapping{3}, x, ms, s)
	#
	sqrt_σ3_0, sqrt_σ3_e = (ms[1]+ms[2]), (√s-ms[3])
    # 
    sqrt_σ3_i = real(sqrt_σ3_e) #+sign(imag(s))*nextfloat(0.0)
    sqrt_σ3 = x[1] < 0.5 ? 
            sqrt_σ3_0 + ( x[1]     / 0.5)*(sqrt_σ3_i-sqrt_σ3_0) : # straight path
	        sqrt_σ3_i + ((x[1]-0.5)/ 0.5)*(sqrt_σ3_e-sqrt_σ3_i) # straight path
	#
	jacobian = x[1] < 0.5 ?
            (sqrt_σ3_i-sqrt_σ3_0) / 0.5 :
            (sqrt_σ3_e-sqrt_σ3_i) / 0.5
	# 
	jacobian *= 2sqrt_σ3 # since the integration is in √σ
    #
	σ3 = sqrt_σ3^2
	# 
	σ2_0, σ2_e = σ2of3_pm(σ3, ms^2, s)
	σ2 = σ2_0 + x[2]*(σ2_e-σ2_0) # straight path
	#
	jacobian *= (σ2_e-σ2_0)
	# 
	return (σ3,σ2), jacobian
end

function mapdalitz(method::HookSqrtDalitzMapping{2}, x, ms, s)
	#
	sqrt_σ2_0, sqrt_σ2_e = (ms[1]+ms[3]), (√s-ms[2])
    # 
    sqrt_σ2_i = real(sqrt_σ2_e)+sign(imag(s))*nextfloat(0.0)
    sqrt_σ2 = x[1] < 0.5 ? 
            sqrt_σ2_0 + ( x[1]     / 0.5)*(sqrt_σ2_i-sqrt_σ2_0) : # straight path
	        sqrt_σ2_i + ((x[1]-0.5)/ 0.5)*(sqrt_σ2_e-sqrt_σ2_i) # straight path
	#
	jacobian = x[1] < 0.5 ?
            (sqrt_σ2_i-sqrt_σ2_0) / 0.5 :
            (sqrt_σ2_e-sqrt_σ2_i) / 0.5
	# 
	jacobian *= 2sqrt_σ2 # since the integration is in √σ
    #
	σ2 = sqrt_σ2^2
	# 
	σ3_0, σ3_e = σ3of2_pm(σ2, ms^2, s)
	σ3 = σ3_0 + x[2]*(σ3_e-σ3_0) # straight path
	#
	jacobian *= (σ3_e-σ3_0)
	# 
	return (σ3,σ2), jacobian
end