###################################################################

#            _|                  _|              _|                _|  _|  
#    _|_|_|  _|_|_|    _|  _|_|        _|_|_|  _|_|_|_|    _|_|_|  _|  _|  
#  _|        _|    _|  _|_|      _|  _|_|        _|      _|    _|  _|  _|  
#  _|        _|    _|  _|        _|      _|_|    _|      _|    _|  _|  _|  
#    _|_|_|  _|    _|  _|        _|  _|_|_|        _|_|    _|_|_|  _|  _|  


function CrystallBall(x,α,n,xbar,σ) 
    μ = (x-xbar)/σ

    C = n/α/(n-1)*exp(-α^2/2)
    D = sqrt(2π)*erf(α/sqrt(2))
    N = 1/(2C+D)/σ
    # 
    abs(μ) < α && return N*exp(-(x-xbar)^2/(2σ^2))
    A = (n/α)^n * exp(-α^2/2)
    B = n/α - α
    return N*A*(B+abs(x-xbar)/σ)^(-n)
end


import AlgebraPDF: AbstractPDF, pars, lims, func, updatevalueorflag
struct convCrystallBall{T} <: AbstractPDF{1}
    source::T
    σ::Float64
    α::Float64
    n::Float64
end

pars(d::convCrystallBall, isfree::Bool) = pars(d.source, isfree)
lims(d::convCrystallBall) = lims(d.source)
updatevalueorflag(d::convCrystallBall, s::Symbol, isfree::Bool, v=getproperty(pars(d),s)) = 
    convCrystallBall(
        AlgebraPDF.ispar(d.source,s) ? updatevalueorflag(d.source,s,isfree,v) : d.source,
        d.σ, d.α, d.n)

function func(d::convCrystallBall, x::NumberOrTuple; p=pars(d))
    σ = func(d.σ, x; p)
    g(z) = CrystallBall(z,d.α,d.n,0,d.σ)
    f(z) = func(d.source, z; p)
    return quadgk(y->f(x-y) * g(y), -7*σ, +7*σ)[1]
end

# # tests
# t = Normalized(FunctionWithParameters((x;p)->x>0,∅), (-1,1))
# tc = convGauss(t,207.6e-3)
# tb = convCrystallBall(t,207.6e-3,1.33,4.58)

# plot(t)
# plot!(tc)
# plot!(tb)
