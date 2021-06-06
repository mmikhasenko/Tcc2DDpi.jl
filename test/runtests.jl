using X2DDpi
using Test

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