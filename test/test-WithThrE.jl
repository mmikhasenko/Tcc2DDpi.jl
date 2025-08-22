# Test file for WithThrE implementation
# Run with: julia --project -e "include(\"test/test-WithThrE.jl\")"

using Test
using X2DDpi
using Parameters

# Test 1: WithThrE struct creation
@testset "WithThrE struct" begin
    magic_real = WithThrE(1.5, 7.0)
    @test magic_real.val == 1.5
    @test magic_real.threshold == 7.0
    @test typeof(magic_real) == WithThrE{Float64}

    # Test constructor with complex value
    magic_complex = WithThrE(1.5 + 0.1im, 7.0)
    @test typeof(magic_complex) == Complex{WithThrE{Float64}}
    @test magic_complex.re.val == 1.5
    @test magic_complex.im.val == 0.1
    @test magic_complex.re.threshold == 7.0
    @test magic_complex.im.threshold == 7.0
end

# Test 2: e2m function with WithThrE
@testset "e2m with WithThrE" begin
    magic_real = WithThrE(1.5, 7.0)
    result = e2m(magic_real)
    expected = 7.0 + 1.5 * 1e-3
    @test result ≈ expected
end

# Test 3: e2m function with Complex{WithThrE}
@testset "e2m with Complex{WithThrE}" begin
    magic_complex = complex(WithThrE(1.5, 7.0), WithThrE(0.1, 7.0))
    result = e2m(magic_complex)
    expected_real = 7.0 + 1.5 * 1e-3
    expected_imag = 0.1 * 1e-3
    expected = expected_real + expected_imag * 1im
    @test result ≈ expected
end

# Test 4: Backward compatibility with regular numbers
@testset "Backward compatibility" begin
    # Test with regular Real
    result_real = e2m(1.5)
    expected_real = (mD⁰ + mDˣ⁺) + 1.5 * 1e-3
    @test result_real ≈ expected_real

    # Test with regular Complex
    result_complex = e2m(1.5 + 0.1im)
    expected_complex = (mD⁰ + mDˣ⁺) + 1.5 * 1e-3 + 0.1 * 1e-3 * 1im
    @test result_complex ≈ expected_complex
end

# Test 5: Real/Complex dispatch compatibility
@testset "Dispatch compatibility" begin
    # Define test functions that expect Real/Complex
    function test_real_dispatch(e::Real)
        return "Real: $(typeof(e))"
    end

    function test_complex_dispatch(e::Complex)
        return "Complex: $(typeof(e))"
    end

    # Test with WithThrE
    magic_real = WithThrE(1.5, 7.0)
    result_real = test_real_dispatch(magic_real)
    @test result_real == "Real: WithThrE{Float64}"

    # Test with Complex{WithThrE}
    magic_complex = complex(WithThrE(1.5, 7.0), WithThrE(0.1, 7.0))
    result_complex = test_complex_dispatch(magic_complex)
    @test result_complex == "Complex: Complex{WithThrE{Float64}}"
end

# Test 6: k3b and kNR functions with WithThrE
@testset "k3b and kNR with WithThrE" begin
    magic_real = WithThrE(1.5, 7.0)

    # Test k3b
    result_k3b = k3b(magic_real)
    @test typeof(result_k3b) <: Complex

    # Test kNR
    result_kNR = kNR(magic_real)
    @test typeof(result_kNR) <: Complex
end

# Test 7: Different thresholds
@testset "Different thresholds" begin
    threshold1 = 3.0
    threshold2 = 7.0

    magic1 = WithThrE(1.5, threshold1)
    magic2 = WithThrE(1.5, threshold2)

    result1 = e2m(magic1)
    result2 = e2m(magic2)

    expected1 = threshold1 + 1.5 * 1e-3
    expected2 = threshold2 + 1.5 * 1e-3

    @test result1 ≈ expected1
    @test result2 ≈ expected2
    @test result1 != result2  # Different thresholds should give different results
end

# Test 8: Conversion operations (what we actually need)
@testset "Conversion operations" begin
    m1 = WithThrE(1.0, 3.0)

    # Test conversion to Float64
    converted = convert(Float64, m1)
    @test converted isa Float64
    @test converted == 1.0

    # Test comparison (which requires conversion)
    @test m1 < 2.0
    @test m1 == 1.0
end


