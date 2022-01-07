using X2DDpi
using Test
using QuadGK

@testset "scalar production" begin
	s,s12,s13,msq = 1,2,3,(2,3,5)
	v = (;s,s12,s13,msq)
	@test X2DDpi.p12_p3(v)+v.s12 == X2DDpi.p_p12(v)
	@test X2DDpi.p13_p2(v)+v.s13 == X2DDpi.p_p13(v)
	@test X2DDpi.p12_p2(v) + X2DDpi.p2_p3(v) == X2DDpi.p_p2(v)
	@test X2DDpi.p13_p2(v) + msq[2] == X2DDpi.p_p2(v)
end

@testset "ABCDEG" begin
	s,s12,s13,msq = 1,2,3,(2,3,3)
	v = (;s,s12,s13,msq)
	@test X2DDpi.A(v) != 0
	@test X2DDpi.B(v) != 0
	@test X2DDpi.C(v) != 0
	@test X2DDpi.C((;s,s12,s13,msq)) ≈ X2DDpi.C((;s12=s13,s13=s12,s,msq))
	@test X2DDpi.D(v) != 0
	@test X2DDpi.E(v) != 0
	@test X2DDpi.G(v) != 0
end



struct TestCh{T} <: X2DDpi.AbstractxDD
	ms::T
end
import X2DDpi:decay_matrix_element_squared
decay_matrix_element_squared(d::TestCh,s,σ3,σ2) = 1.0

@testset "implementation of the phase-space integral" begin

	e_test = 0.0
	ms_test = (m1=mπ⁺,m2=mD⁰,m3=mD⁰)

	t = TestCh(ms_test)
	ρ_2dim = ρ_thr(t,e_test)

	s = e2m(e_test)^2
	ms = ms_test
	ρ_1dim = quadgk(σ1->sqrt(X2DDpi.λ(s,ms.m1^2,σ1)*X2DDpi.λ(σ1,ms.m2^2,ms.m3^2))/σ1,
		(ms.m2+ms.m3)^2, (sqrt(s)-ms.m1)^2)[1]  / (8π)^2 / (2π*s)

	@test abs(ρ_2dim-ρ_1dim)/ρ_1dim < 1e-4
end



@testset "Break up momentum of Dˣ⁺D⁰" begin

	Eᵦ = X2DDpi.Eᵦˣ⁺
	# 
	@test kNR(Eᵦ) + 1 ≈ 1.0 + 0.0im
	@test k3b(Eᵦ) + 1 ≈ 1.0 + 0.0im

	@test let
		a,b = kNR.((Eᵦ - 1e3im*ΓDˣ⁺/4) .+ [-1,1].*1e-6)
		abs(1-b/a) > 1.0
	end
	@test let
		a,b = kNR.((Eᵦ + 1e3im*ΓDˣ⁺/4) .+ [-1,1].*1e-6)
		abs(1-b/a) < 1e-4
	end
	@test let
		a,b = k3b.((Eᵦ - 1e3im*ΓDˣ⁺/4) .+ [-1,1].*1e-6)
		abs(1-b/a) > 1.0
	end
	@test let
		a,b = k3b.((Eᵦ + 1e3im*ΓDˣ⁺/4) .+ [-1,1].*1e-6)
		abs(1-b/a) < 1e-4
	end
end







@testset "Cauchy integrals" begin
    @test cauchy(x->x, 1.0, 0.1) ≈ 1.0 + 0.0im
    @test cauchy′(x->x, 0, 0.1) ≈ 1.0 + 0.0im
    @test cauchy′′(x->x^2, 0, 0.1) ≈ 2.0 + 0.0im
end
# 
@testset "Circular integrals" begin
    @test circleintegral(x->x^2, 0.1) + 1 ≈ 1.0 + 0.0im
    @test circleintegral(x->1/x, 0.1) ≈ 1.0 + 0.0im
    @test circleintegral(x->1/x, 0.1) ≈ circleintegral(x->1/x, 0.01)
end
