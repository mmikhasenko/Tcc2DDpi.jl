
circleintegral(f,r) = quadgk(ϕ->f(r*cis(ϕ))*r*cis(ϕ), -π, π)[1]/(2π) # cancel 1im

cauchy(F,x₀,r) = circleintegral(x′->F(x′+x₀)/x′, r)
cauchy′(F,x₀,r) = circleintegral(x′->F(x′+x₀)/x′^2, r)

"""
    f(s) / N = a⁻¹ + r k(s)^2 / 2 - i k(s)

Is calculated using cauchy integral theorem.
"""
function effectiverangeexpansion(f, k, r)
    N = circleintegral(f,r) / circleintegral(k,r) / (-1im)
    # 
    f̂₀  = cauchy(  x->f(x)/N + 1im*k(x), 0.0, r)
    f̂′  = cauchy′( x->f(x)/N + 1im*k(x), 0.0, r)
    # 
    k²₀  = cauchy(  x->k(x)^2, 0.0, r)
    k²′  = cauchy′( x->k(x)^2, 0.0, r)
    # 
    a⁻¹ = f̂₀
    r = 2 * f̂′/k²′
    #
    return (; a⁻¹, r)
end