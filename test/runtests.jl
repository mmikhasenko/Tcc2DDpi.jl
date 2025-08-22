using X2DDpi
using Test
using QuadGK
using Parameters


@testset "Break up momentum of Dˣ⁺D⁰" begin

    @test kNR(Eᵦˣ⁺) + 1 ≈ 1.0 + 0.0im
    @test k3b(Eᵦˣ⁺) + 1 ≈ 1.0 + 0.0im

    @test let
        a, b = kNR.((Eᵦˣ⁺ - 1e3im * ΓDˣ⁺ / 4) .+ [-1, 1] .* 1e-6)
        abs(1 - b / a) > 1.0
    end
    @test let
        a, b = kNR.((Eᵦˣ⁺ + 1e3im * ΓDˣ⁺ / 4) .+ [-1, 1] .* 1e-6)
        abs(1 - b / a) < 1e-4
    end
    @test let
        a, b = k3b.((Eᵦˣ⁺ - 1e3im * ΓDˣ⁺ / 4) .+ [-1, 1] .* 1e-6)
        abs(1 - b / a) > 1.0
    end
    @test let
        a, b = k3b.((Eᵦˣ⁺ + 1e3im * ΓDˣ⁺ / 4) .+ [-1, 1] .* 1e-6)
        abs(1 - b / a) < 1e-4
    end
end


import X2DDpi.xDDPhaseSpace: HookSqrtDalitzMapping, LinearDalitzMapping
import X2DDpi.xDDPhaseSpace: BW

@testset "Dalitz Mappings" begin
    e = -0.1 - 0.1im
    ch = DˣD((m1 = mπ⁺, m2 = mD⁰, m3 = mD⁰), BW(m = mDˣ⁺, Γ = ΓDˣ⁺))
    #
    method1 = HookSqrtDalitzMapping{3}()
    set1 = ChannelWithIntegrationMethod(ch, method1)
    value1 = ρ_thr(set1, e)
    # 
    method2 = LinearDalitzMapping()
    set2 = ChannelWithIntegrationMethod(ch, method1)
    value2 = ρ_thr(set2, e)
    # 
    @test value1 ≈ value2
end


# Effective range Expansion
@testset "Circular integrals" begin
    @test circleintegral(x -> x^2, 0, 0.1) + 1 ≈ 1.0 + 0.0im
    @test circleintegral(x -> 1 / x, 0, 0.1) ≈ 1.0 + 0.0im
    @test circleintegral(x -> 1 / x, 0, 0.1) ≈ circleintegral(x -> 1 / x, 0, 0.01)
end
# 
@testset "Cauchy integrals" begin
    @test cauchyintegral(x -> x, 1, 0.1) ≈ 1.0 + 0.0im
    @test cauchyintegral′(x -> x, 0, 0.1) ≈ 1.0 + 0.0im
    @test cauchyintegral′′(x -> x^2, 0, 0.1) ≈ 2.0 + 0.0im
end


@testset "Cauchy sums" begin
    v0 = X2DDpi.cauchyintegral(x -> x^10, 2.0, 0.1)
    v1 = X2DDpi.cauchyintegral(x -> x^10, 2.0, CircularIntegral(0.1))
    v2 = X2DDpi.cauchyintegral(x -> x^10, 2.0, CircularSum(0.1, 250))
    v3 = X2DDpi.cauchyintegral(x -> x^10, 2.0, CircularSum(0.1, 6))
    # 
    @test v0 == v1
    @test v1 ≈ v2
    @test !(v3 ≈ v2)
end

circleintegral(x -> x^2, 2.0, 0.1, 1000)

@testset "Effective range Expansion" begin
    a₀_fm = -7.0 # fm
    r₀_fm = 3.0 # fm
    # 
    a⁻¹ = 1 / (a₀_fm / (1e-3 * fm_times_mev))
    r = r₀_fm / (1e-3 * fm_times_mev)
    #
    k = k3b
    D(e) = a⁻¹ + r / 2 * k(e)^2 - 1im * k(e)
    #
    # expansion is at 
    Eᵦ = m2e(sqrt(mDˣ⁺^2 - 1im * mDˣ⁺ * ΓDˣ⁺) + mD⁰)
    # 
    effrangepars =
        effectiverangeexpansion(
            Δe -> D(Δe + Eᵦ),
            Δe -> k(Eᵦ + Δe),
            abs(imag(Eᵦ)) / 20)
    # 
    @unpack a_fm, r_fm = tophysicsunits(effrangepars)
    # 
    @test a_fm ≈ a₀_fm
    @test real(r_fm) ≈ r₀_fm
end

@testset "Channel With Integration Method" begin

    e_safe = Eᵦˣ⁺ - 0.1
    e_pole = Eᵦˣ⁺ + 0.1

    ch1 = πDD((m1 = mπ⁺, m2 = mD⁰, m3 = mD⁰), BW(m = mDˣ⁺, Γ = ΓDˣ⁺), BW(m = mDˣ⁺, Γ = ΓDˣ⁺))

    set3 = ChannelWithIntegrationMethod(ch1, HookSqrtDalitzMapping{3}())
    set2 = ChannelWithIntegrationMethod(ch1, HookSqrtDalitzMapping{2}())

    @test set3 isa ChannelWithIntegrationMethod{T, HookSqrtDalitzMapping{3}} where T <: πDD
    @test set2 isa ChannelWithIntegrationMethod{T, HookSqrtDalitzMapping{2}} where T <: πDD

    full_safe = ρ_thr(ch1, e_safe)
    ρ3_safe = ρ_thr(set3, e_safe)
    ρ2_safe = ρ_thr(set2, e_safe)
    sum23_safe = ρ2_safe + ρ3_safe
    @test abs((full_safe - sum23_safe) / full_safe) < 1e-4

    full_pole = ρ_thr(ch1, e_pole)
    ρ3_pole = ρ_thr(set3, e_pole)
    ρ2_pole = ρ_thr(set2, e_pole)
    sum23_pole = ρ2_pole + ρ3_pole
    @test abs((full_pole - sum23_pole) / full_pole) > 1e-4
end
