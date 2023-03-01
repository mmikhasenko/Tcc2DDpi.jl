
ere(k; a⁻¹, r, N) = N*(a⁻¹ + r * k^2 / 2 - 1im* k)
ere(k, p) = ere(k, a⁻¹ = p.a⁻¹, r = p.r, N = p.N)

struct EffectiveRangeFit{T<:AbstractArray, S} <: EffectiveRangeExpansionMethod
    evaluationrange::T
    start::S
end

function effectiverangeexpansion(f, k, method::EffectiveRangeFit)

    @unpack evaluationrange, start = method
    dD(x; a⁻¹, r, N) = f(x)-ere(k(x); a⁻¹, r, N)
    n = sum(abs2, f.(evaluationrange))

    function mismatch(x)
        a⁻¹, r, N = x[1:3] .+ 1im .* x[4:6]
        χ² = sum(abs2, dD.(evaluationrange; a⁻¹, r, N))
        return χ² / n
    end

    p = NamedTuple{(:a⁻¹, :r, :N)}(start)
    x0 = [real([p...])..., imag.([p...])...]
    optim_result = Optim.optimize(mismatch, x0, Optim.BFGS())
    
    xmin = optim_result.minimizer
    a⁻¹, r, N = xmin[1:3] .+ 1im .* xmin[4:6]

    return (; a⁻¹, r, N, minimum=optim_result.minimum, χ0 = mismatch(x0))
end

