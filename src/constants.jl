
# masses
const mπ⁰ = 0.1349768
const mπ⁺ = 0.13957039
const mD⁰ = 1.86483
const mD⁺ = 1.86965
# 
const mDˣ⁺ = 1.86483 + 145.4258e-3 # m(D) + Δm(D*,D) from PDG
const mDˣ⁰ = 2.00685
# 
const ΓDˣ⁺ = 83.4e-6
const ΓDˣ⁰ = 55.2e-6
const mγ = 0.0
# 

const thr_D⁰Dˣ⁺ = mD⁰ + mDˣ⁺

# Backward compatibility
e2m(e) = (mD⁰ + mDˣ⁺) + e * 1e-3
m2e(m) = (m - mD⁰ - mDˣ⁺) * 1e3

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



# NR quantities
const μDˣ⁺D⁰ = mDˣ⁺ * mD⁰ / (mDˣ⁺ + mD⁰)
const Eᵦˣ⁺ = m2e(sqrt(mDˣ⁺^2 - 1im * mDˣ⁺ * ΓDˣ⁺) + mD⁰)

# coupling constants
const h² = 20.13e-3
const f² = 282.42
const μ₀ = -3.77
const μ₊ = 1.0

# Δllh = Δ² / 2σ²
const ΔNLL_90CL = 1.352 # 90% → Δ=1.644σ; 
const ΔNLL_95CL = 1.921 # 95% → Δ=1.959σ
const ΔNLL_68CL = 0.5   # 68% → Δ=1σ

# units conversion
const fm_times_mev = 197.3269804
# 
function tophysicsunits(p::NamedTuple)
    @unpack a⁻¹, r = p
    a_fm = 1e-3 * fm_times_mev / real(p.a⁻¹)
    r_fm = 1e-3 * fm_times_mev * p.r
    (; a_fm, r_fm)
end


"""
    k3b(e)

Full relativistic expression for Dˣ⁺D⁰ break-up momentum
"""
function k3b(e)
    m = e2m(e)
    M = sqrt(mDˣ⁺^2 - 1im * mDˣ⁺ * ΓDˣ⁺) # taken at the Dˣ⁺ pole
    p = cis(π / 4) * sqrt((m - (M + mD⁰)) * cis(-π / 2)) *
        sqrt(m + (M + mD⁰)) * sqrt(m - (M - mD⁰)) * sqrt(m + (M - mD⁰)) / (2 * m)  # branch cut down
    return p
end

"""
    kNR(e)

Non-relativistic expression for Dˣ⁺D⁰ break-up momentum
"""
function kNR(e)
    μ = μDˣ⁺D⁰
    Eᵦ = Eᵦˣ⁺
    # 
    p = cis(π / 4) * sqrt(2μ * (e - Eᵦ) / 1e3 * cis(-π / 2)) # branch cut down
    return p
end
