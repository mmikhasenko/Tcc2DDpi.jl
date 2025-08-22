abstract type AbstractxDD end
# 
# methods with general implementation:
#  - decay_matrix_element_squared
#  - integrand_mapped_thr
#  - mapdalitzmethod
#  - ρ_thr
#  - masses
# 
#  no general implementation:
#  - ρ_tb
#  - obj2nt
#  - branch_points


mapdalitzmethod(d::AbstractxDD) = HookSqrtDalitzMapping{3}()
masses(d::AbstractxDD) = d.ms

#  _|              _|                                                _|      _|                      
#      _|_|_|    _|_|_|_|    _|_|      _|_|_|  _|  _|_|    _|_|_|  _|_|_|_|        _|_|    _|_|_|    
#  _|  _|    _|    _|      _|_|_|_|  _|    _|  _|_|      _|    _|    _|      _|  _|    _|  _|    _|  
#  _|  _|    _|    _|      _|        _|    _|  _|        _|    _|    _|      _|  _|    _|  _|    _|  
#  _|  _|    _|      _|_|    _|_|_|    _|_|_|  _|          _|_|_|      _|_|  _|    _|_|    _|    _|  
#                                          _|                                                        
#                                      _|_|                                                          

function integrand_mapped_thr(d::AbstractxDD,s,x)
	method = mapdalitzmethod(d)
	(σ3,σ2), jac = mapdalitz(method, x, masses(d), s)
	return decay_matrix_element_squared(d,s,σ3,σ2) / (2π*s) * jac
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
