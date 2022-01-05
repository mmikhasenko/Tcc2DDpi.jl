using QuadGK
using Test

circleintegral(f,r) = quadgk(ϕ->f(r*cis(ϕ))*r*cis(ϕ), -π, π)[1]/(2π) # cancel 1im

cauchy(F,x₀,r) = circleintegral(x′->F(x′+x₀)/x′, r)
cauchy′(F,x₀,r) = circleintegral(x′->F(x′+x₀)/x′^2, r)
cauchy′′(F,x₀,r) = circleintegral(x′->2F(x′+x₀)/x′^3, r)

