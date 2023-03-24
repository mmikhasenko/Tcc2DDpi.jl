
struct EffectiveRangeFit{T<:AbstractArray,S} <: EffectiveRangeExpansionMethod
    evaluationrange::T
    start::S
end

function effectiverangeexpansion(f, k, method::EffectiveRangeFit)

    @unpack evaluationrange, start = method
    ev = evaluationrange
    yv = f.(ev)
    # 
    n = sum(abs2, yv)

    function mismatch(x)
        a⁻¹, r, N = x[1:3] .+ 1im .* x[4:6]
        yv′ = ere.(k.(ev); a⁻¹, r, N)
        χ² = sum(abs2, yv - yv′)
        return χ² / n
    end

    p = NamedTuple{(:a⁻¹, :r, :N)}(start)
    x0 = [real.([p...])..., imag.([p...])...]
    optim_result = Optim.optimize(mismatch, x0, Optim.BFGS())

    xmin = optim_result.minimizer
    a⁻¹, r, N = xmin[1:3] .+ 1im .* xmin[4:6]

    return (; a⁻¹, r, N, minimum=optim_result.minimum, χ0=mismatch(x0))
end

