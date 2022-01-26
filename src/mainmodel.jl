
# πDD
struct πDD{T1,T2,T3} <: AbstractxDD
    ms::NamedTuple{(:m1,:m2,:m3),T1}
    R12::T2
    R13::T3
end

# γDD
struct γDD{T1,T2,T3} <: AbstractxDD
    ms::NamedTuple{(:m1,:m2,:m3),T1}
    R12::T2
    R13::T3
end
# 





function covertapply(𝔐²,d::AbstractxDD,s,σ3,σ2)
	msq = masses(d)^2
	v = (;s,s12=σ3,s13=σ2,msq)
# 	
	J₁₂ᴵ, J₁₂ᴵᴵ = Jᴵ(σ3,d.R12), Jᴵᴵ(σ3,d.R12)
	J₁₃ᴵ, J₁₃ᴵᴵ = Jᴵ(σ2,d.R13), Jᴵᴵ(σ2,d.R13)
# 	
	return 𝔐²(v, (J₁₂ᴵ, J₁₂ᴵᴵ), (J₁₃ᴵ, J₁₃ᴵᴵ))
end


# πDD
"""
	πDD_𝔐²_nonana3(v, J₁₂, J₁₃ᴵ)

The function has the both poles in σ₃,
	and only the top pole in σ₂
"""
function πDD_𝔐²_nonana3(v, (J₁₂ᴵ, J₁₂ᴵᴵ), (J₁₃ᴵ, J₁₃ᴵᴵ))
	𝔐² = A(v) * J₁₂ᴵ * J₁₂ᴵᴵ +
             C(v) * J₁₃ᴵ * J₁₂ᴵᴵ
	return f² * 𝔐²/3/4
end

"""
	πDD_𝔐²_nonana2(v, J₁₂, J₁₃ᴵ)

The function has the both poles in σ₂,
	and only the top pole in σ₃
"""
function πDD_𝔐²_nonana2(v, (J₁₂ᴵ, J₁₂ᴵᴵ), (J₁₃ᴵ, J₁₃ᴵᴵ))
	𝔐² = B(v) * J₁₃ᴵ * J₁₃ᴵᴵ +
             C(v) * J₁₂ᴵ * J₁₃ᴵᴵ
	return f² * 𝔐²/3/4
end

decay_matrix_element_squared(d::πDD,s,σ3,σ2) = covertapply(
		(v,J₁₂,J₁₃)->πDD_𝔐²_nonana3(v,J₁₂,J₁₃)+
		             πDD_𝔐²_nonana2(v,J₁₂,J₁₃),
	    d, s,σ3,σ2)

# γDD

function γDD_𝔐²_nonana3(v, (J₁₂ᴵ, J₁₂ᴵᴵ), (J₁₃ᴵ, J₁₃ᴵᴵ))
	_p1_p2 = p1_p2(v)
	_p1_p3 = p1_p3(v)
	_G = G(v)
# 	
	𝔐² = μ₊^2 * (_p1_p2^2+_G) * J₁₂ᴵ * J₁₂ᴵᴵ -
			μ₊*μ₀ * (_p1_p2*_p1_p3 - _G) * J₁₃ᴵ * J₁₂ᴵᴵ
# 	
	return h² * 𝔐² / 3
end

function γDD_𝔐²_nonana2(v, (J₁₂ᴵ, J₁₂ᴵᴵ), (J₁₃ᴵ, J₁₃ᴵᴵ))
	_p1_p2 = p1_p2(v)
	_p1_p3 = p1_p3(v)
	_G = G(v)
# 
	𝔐² = μ₀^2 * (_p1_p3^2+_G) * J₁₃ᴵ * J₁₃ᴵᴵ -
			μ₊*μ₀ * (_p1_p2*_p1_p3 - _G) * J₁₂ᴵ * J₁₃ᴵᴵ
# 
	return h² * 𝔐² / 3
end

decay_matrix_element_squared(d::γDD,s,σ3,σ2) = covertapply(
		(v,J₁₂,J₁₃)->γDD_𝔐²_nonana3(v,J₁₂,J₁₃)+
		             γDD_𝔐²_nonana2(v,J₁₂,J₁₃),
	    d, s,σ3,σ2)

branch_points(d::Union{πDD,γDD}) = (
	m2e(d.ms[3] + sqrt(pole_position(d.R12))),
	m2e(d.ms[2] + sqrt(pole_position(d.R13))))
# 

function ρ_tb(d::Union{πDD,γDD}, e::Real)
	M,m = d.R13.m, d.ms.m2
	sqrts = e2m(e)
	sqrts < M+m ? 0.0 :
    	sqrt(λ(e2m(e)^2,M^2,m^2))/e2m(e)^2
end	




#                                _|  _|  _|                        _|      _|                      
#    _|_|_|    _|_|    _|  _|_|      _|      _|_|_|_|    _|_|_|  _|_|_|_|        _|_|    _|_|_|    
#  _|_|      _|_|_|_|  _|_|      _|  _|  _|      _|    _|    _|    _|      _|  _|    _|  _|    _|  
#      _|_|  _|        _|        _|  _|  _|    _|      _|    _|    _|      _|  _|    _|  _|    _|  
#  _|_|_|      _|_|_|  _|        _|  _|  _|  _|_|_|_|    _|_|_|      _|_|  _|    _|_|    _|    _|  


#
obj2nt(ch::Union{πDD,γDD}) =
    (type=string(typeof(ch)),
        ms=ch.ms,
		R12=obj2nt(ch.R12), R13=obj2nt(ch.R13))
#
obj2nt(nt::NamedTuple) = nt

# deserilization
function constructchannel(ser::NamedTuple{(:R12, :R13, :ms, :type),T} where T)
    type = eval(Meta.parse(ser.type))
    ms = NamedTuple{(:m1,:m2,:m3)}(ser.ms)
	R12 = constructlineshape(ser.R12)
	R13 = constructlineshape(ser.R13)
    type(ms, R12, R13)
end

constructlineshape(R::NamedTuple{(:m,:Γ),T} where T) = R
# 
function constructlineshape(R::NamedTuple)
    type = eval(Meta.parse(R.type))
	keysnotype =  filter(x->x!=:type, keys(R))
	args = NamedTuple{keysnotype}(R)
	type(args...)
end
