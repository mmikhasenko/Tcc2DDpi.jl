abstract type AbstractxDD end
#
function ρ_tb(d::AbstractxDD, e::Real)
	M,m = d.R13.m, d.ms.m2
	sqrts = e2m(e)
	sqrts < M+m ? 0.0 :
    	sqrt(λ(e2m(e)^2,M^2,m^2))/e2m(e)^2
end	

# πDD
struct πDD{T1,T2,T3} <: AbstractxDD
    ms::NamedTuple{(:m1,:m2,:m3),T1}
    R12::NamedTuple{T2}
    R13::NamedTuple{T3}
end

# γDD
struct γDD{T1,T2,T3} <: AbstractxDD
    ms::NamedTuple{(:m1,:m2,:m3),T1}
    R12::NamedTuple{T2}
    R13::NamedTuple{T3}
end


function decay_matrix_element(d::πDD,s,σ3,σ2)
	msq = d.ms^2
	v = (;s,s12=σ3,s13=σ2,msq)
# 	
	J12_I, J12_II = J_I(σ3,d.R12), J_II(σ3,d.R12)
	J13_I, J13_II = J_I(σ2,d.R13), J_II(σ2,d.R13)
# 	
	frakM = A(v) * J12_I * J12_II +
			B(v) * J13_I * J13_II +
			C(v) * (J13_I * J12_II +  J12_I * J13_II)
	f²*frakM/3/4
end


function decay_matrix_element(d::γDD,s,σ3,σ2)
	msq = d.ms^2
	v = (;s,s12=σ3,s13=σ2,msq)
# 	
	J12_I, J12_II = J_I(σ3,d.R12), J_II(σ3,d.R12)
	J13_I, J13_II = J_I(σ2,d.R13), J_II(σ2,d.R13)
# 	
	_p1_p2 = p1_p2(v)
	_p1_p3 = p1_p3(v)
	_G = G(v)
# 	
	frakM = μ₊^2 * (_p1_p2^2+_G) * J12_I * J12_II +
			μ₀^2 * (_p1_p3^2+_G) * J13_I * J13_II -
			μ₊*μ₀ * (_p1_p2*_p1_p3 - _G) *
		(J13_I * J12_II +  J12_I * J13_II)
# 	
	h²*frakM/3
end

branch_points(d::AbstractxDD) = (
	m2e(d.ms[3] + sqrt(pole_position(d.R12))),
	m2e(d.ms[2] + sqrt(pole_position(d.R13))))


# serilization
obj2nt(ch::AbstractxDD) =
    (type=string(typeof(ch)),
        ms=ch.ms, R12=ch.R12, R13=ch.R13)

# deserilization
function constructchannel(ser::NamedTuple{(:R12, :R13, :ms, :type),T} where T)
    type = eval(Meta.parse(ser.type))
    ms = NamedTuple{(:m1,:m2,:m3)}(ser.ms)
    type(ms, ser.R12, ser.R13)
end
