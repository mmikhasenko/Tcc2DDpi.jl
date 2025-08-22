# WithThrE type for channel-specific thresholds
struct WithThrE{T <: Real} <: Real
    val::T
    threshold::Float64
end

# Make WithThrE behave like a number for conversion and comparison only
Base.convert(::Type{T}, e::WithThrE) where T <: Number = convert(T, e.val)
Base.promote_rule(::Type{<:Number}, ::Type{<:WithThrE}) = promote_rule(Number, Number)

# Constructor for complex values
WithThrE(val::Complex, threshold::Float64) = complex(WithThrE(real(val), threshold), WithThrE(imag(val), threshold))

# conversion functions
e2m(e::WithThrE) = e.threshold + e.val * 1e-3
e2m(e::Complex{<:WithThrE}) = e.re.threshold + e.re.val * 1e-3 + e.im.val * 1e-3 * 1im

