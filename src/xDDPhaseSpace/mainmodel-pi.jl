struct πDD{T1, T2, T3} <: AbstractxDD
    ms::NamedTuple{(:m1, :m2, :m3), T1}
    R12::T2
    R13::T3
end

function covertapply(𝔐², d::AbstractxDD, s, σ3, σ2)
    msq = masses(d)^2
    v = (; s, s12 = σ3, s13 = σ2, msq)
    # 	
    J₁₂ᴵ, J₁₂ᴵᴵ = Jᴵ(σ3, d.R12), Jᴵᴵ(σ3, d.R12)
    J₁₃ᴵ, J₁₃ᴵᴵ = Jᴵ(σ2, d.R13), Jᴵᴵ(σ2, d.R13)
    # 	
    return 𝔐²(v, (J₁₂ᴵ, J₁₂ᴵᴵ), (J₁₃ᴵ, J₁₃ᴵᴵ))
end

const f² = 282.42


"""
πDD_𝔐²_nonana3(v, J₁₂, J₁₃ᴵ)

The function has the both poles in σ₃,
and only the top pole in σ₂
"""
function πDD_𝔐²_nonana3(v, (J₁₂ᴵ, J₁₂ᴵᴵ), (J₁₃ᴵ, J₁₃ᴵᴵ))
    𝔐² = A(v) * J₁₂ᴵ * J₁₂ᴵᴵ +
         C(v) * J₁₃ᴵ * J₁₂ᴵᴵ
    return f² * 𝔐² / 3 / 4
end

"""
πDD_𝔐²_nonana2(v, J₁₂, J₁₃ᴵ)

The function has the both poles in σ₂,
and only the top pole in σ₃
"""
function πDD_𝔐²_nonana2(v, (J₁₂ᴵ, J₁₂ᴵᴵ), (J₁₃ᴵ, J₁₃ᴵᴵ))
    𝔐² = B(v) * J₁₃ᴵ * J₁₃ᴵᴵ +
         C(v) * J₁₂ᴵ * J₁₃ᴵᴵ
    return f² * 𝔐² / 3 / 4
end

decay_matrix_element_squared(d::πDD, s, σ3, σ2) = covertapply(
    (v, J₁₂, J₁₃) -> πDD_𝔐²_nonana3(v, J₁₂, J₁₃) +
                     πDD_𝔐²_nonana2(v, J₁₂, J₁₃),
    d, s, σ3, σ2)

branch_points(d::πDD) = (
    m2e(d.ms[3] + sqrt(pole_position(d.R12))),
    m2e(d.ms[2] + sqrt(pole_position(d.R13))))

function ρ_tb(d::πDD, e::Real)
    M, m = d.R13.m, d.ms.m2
    sqrts = e2m(e)
    sqrts < M + m ? 0.0 :
    sqrt(λ(e2m(e)^2, M^2, m^2)) / e2m(e)^2
end

#                                _|  _|  _|                        _|      _|                      
#    _|_|_|    _|_|    _|  _|_|      _|      _|_|_|_|    _|_|_|  _|_|_|_|        _|_|    _|_|_|    
#  _|_|      _|_|_|_|  _|_|      _|  _|  _|      _|    _|    _|    _|      _|  _|    _|  _|    _|  
#      _|_|  _|        _|        _|  _|  _|    _|      _|    _|    _|      _|  _|    _|  _|    _|  
#  _|_|_|      _|_|_|  _|        _|  _|  _|  _|_|_|_|    _|_|_|      _|_|  _|    _|_|    _|    _|  

#
obj2nt(ch::πDD) =
    (type = string(typeof(ch)),
        ms = ch.ms,
        R12 = obj2nt(ch.R12), R13 = obj2nt(ch.R13))
#
