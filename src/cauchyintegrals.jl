

circleintegral(f, x₀::Real, r::Real) =
    quadgk(ϕ->f(x₀+r*cis(ϕ))*r*cis(ϕ), -7π/8, 9π/8)[1]/(2π) # cancel 1im
# 
circleintegral(f, x₀::Real, r::Real, N::Int) =
    sum(f(x₀+r*cis(ϕ))*r*cis(ϕ) for ϕ in range(-7π/8, 9π/8, length=N+1)[2:end])[1]/N # 2π i cancels
#
cauchyintegral(  F, x₀::Real, r) = circleintegral(x′-> F(x′+x₀)/x′,   0.0, r)
cauchyintegral′( F, x₀::Real, r) = circleintegral(x′-> F(x′+x₀)/x′^2, 0.0, r)
cauchyintegral′′(F, x₀::Real, r) = circleintegral(x′->2F(x′+x₀)/x′^3, 0.0, r)
# 


abstract type CircularIntegralMethod end
struct CircularIntegral <: CircularIntegralMethod
    r::Real
end
struct CircularSum <: CircularIntegralMethod
    r::Real
    N::Int
end

circleintegral(F, x₀::Real, method::CircularIntegral) = circleintegral(F, x₀, method.r)
circleintegral(F, x₀::Real, method::CircularSum) = circleintegral(F, x₀,method.r, method.N)







abstract type EffectiveRangeExpansionMethod end
struct ComplexBranchPointExpansion <: EffectiveRangeExpansionMethod
    cim::CircularIntegralMethod
end

effectiverangeexpansion(f, k, r::Float64) =
    effectiverangeexpansion(f, k, ComplexBranchPointExpansion(CircularIntegral(r)))

"""
    f(s) / N = a⁻¹ + r k(s)^2 / 2 - i k(s)

Is calculated using cauchy integral theorem.
"""
function effectiverangeexpansion(f, k, method::EffectiveRangeExpansionMethod)
    N = circleintegral(f,0.0,method.cim) /
        circleintegral(k,0.0,method.cim) / (-1im)
    # 
    f̂₀  = cauchyintegral(  x->f(x)/N + 1im*k(x), 0.0, method.cim)
    f̂′  = cauchyintegral′( x->f(x)/N + 1im*k(x), 0.0, method.cim)
    # 
    # k²₀  = cauchyintegral(  x->k(x)^2, 0.0, r)
    k²′  = cauchyintegral′( x->k(x)^2, 0.0, method.cim)
    # 
    a⁻¹ = f̂₀
    r = 2 * f̂′/k²′
    #
    return (; a⁻¹, r, N)
end