# WithThrE is a conversion wrapper supplied by xDDPhaseSpace, not an
# arithmetic quantity. Convert to absolute mass before calling mass APIs.
using Test
using X2DDpi

@testset "Explicit reference thresholds" begin
    for threshold in (3.0, 7.0), e in (1.5, 1.5 + 0.1im)
        tagged = WithThrE(e, threshold)
        @test e2m(tagged) ≈ threshold + e / 1000
        @test m2e(e2m(tagged)) ≈ (threshold - thr_D⁰Dˣ⁺) * 1000 + e
    end
    @test convert(Float64, WithThrE(1.0, 3.0)) == 1.0
    @test e2m(1.5) ≈ thr_D⁰Dˣ⁺ + 0.0015
    @test m2e(e2m(1.5 + 0.1im)) ≈ 1.5 + 0.1im
end
