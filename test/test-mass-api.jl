# A smooth spectral density with an analytic dispersion integral. Changing its
# threshold catches accidental dependence on the nominal Tcc reference mass.
struct MassTestChannel <: X2DDpi.AbstractxDD
    ms::NamedTuple{(:m1, :m2, :m3), Tuple{Float64, Float64, Float64}}
end
X2DDpi.ρ_thr(d::MassTestChannel, m::Real) = m <= sum(d.ms) ? 0.0 : 1 - (sum(d.ms) / m)^2
X2DDpi.ρ_thr(d::MassTestChannel, m::Complex) = 1 - (sum(d.ms) / m)^2
X2DDpi.ρ_tb(d::MassTestChannel, m::Real) = X2DDpi.ρ_thr(d, m)

@testset "Mass-parametric interpolation and dispersion" begin
    for threshold in (2.0, 4.7)
        ch = MassTestChannel((m1 = 0.2, m2 = 0.8, m3 = threshold - 1))
        cutoff = threshold + 0.1
        ich = interpolated(ch, cutoff; mstep = 0.0005)
        @test first(ich.itr.knots[1]) ≈ threshold
        @test ich.itr.knots[1][2] - ich.itr.knots[1][1] ≈ 0.0005
        @test ρ_thr(ich, threshold - 0.01) == 0
        @test ρ_thr(ich, threshold) == 0
        for m in (threshold + 0.02525, cutoff, cutoff + 0.2)
            @test isapprox(ρ_thr(ich, m), ρ_thr(ch, m); rtol = 2e-5)
        end
        @test ρ_thr(ich, cutoff - 0.01im) == ρ_thr(ch, cutoff - 0.01im)
        exact(m) = (1 + (threshold^2 / m^2 - 1) * log(1 - m^2 / threshold^2)) / π
        for m in (threshold - 0.1 + 0.0im, cutoff + 0.03im, cutoff - 0.03im)
            @test isapprox(dispersive(ich, m), exact(m); rtol = 2e-5)
        end
        @test isapprox(dispersive(ich, cutoff), exact(cutoff + 1e-9im); rtol = 2e-5)
        @test isapprox(imag(dispersive(ich, cutoff)), ρ_thr(ich, cutoff); rtol = 2e-5)
        @test_throws ArgumentError interpolated(ch, cutoff; mstep = 0)
        @test_throws ArgumentError interpolated(ch, cutoff; mstep = -0.1)
        @test_throws ArgumentError interpolated(ch, threshold - 0.1)
    end
end

reversefields(x) = x
reversefields(x::NamedTuple) = (; (k => reversefields(getproperty(x, k)) for k in reverse(keys(x)))...)

@testset "Cache round trips preserve masses, lineshapes, and couplings" begin
    ms = (m1 = 0.15, m2 = 1.8, m3 = 1.9)
    shapes = (BW(m = 2.1, Γ = 0.001), BW_norm(m = 2.2, Γ = 0.002, n = 0.7),
        BW_Swave(m = 2.1, g = 0.1, ma = 1.8, mb = 0.15), ZeroBW(m = 2.2, Γ = 0.002))
    for R in shapes
        for ch in (πDD(ms, R, shapes[1]), γDD(ms, R, shapes[1], 2.3, -0.8), DˣD(ms, R))
            nt = X2DDpi.obj2nt(ch)
            for ser in (nt, reversefields(nt), d2nt(X2DDpi.JSON.parse(X2DDpi.JSON.json(nt))))
                restored = constructchannel(ser)
                @test X2DDpi.obj2nt(restored) == nt
            end
        end
    end
    model = readjson(joinpath(@__DIR__, "..", "results", "nominal", "model.json"))
    for entry in model["ichannels"]
        ich = interpolated(d2nt(entry))
        restored = interpolated(reversefields(X2DDpi.obj2nt(ich)))
        @test restored.itr.coefs == ich.itr.coefs
        @test restored.itr.knots == ich.itr.knots
        @test X2DDpi.obj2nt(restored.channel) == X2DDpi.obj2nt(ich.channel)
        for m in (e2m(-0.3592), e2m(1.5), ich.m_cutoff + 0.001)
            @test ρ_thr(restored, m) ≈ ρ_thr(ich, m)
        end
        # Direct phase space checks the stored values as well as the coordinates.
        @test isapprox(ρ_thr(ich, e2m(-0.3592)), ρ_thr(ich.channel, e2m(-0.3592)); rtol = 2e-3)
    end
    γ = γDD(ms, shapes[1], shapes[1])
    old = (; type = "γDD", ms, R12 = (Γ = 0.001, m = 2.1), R13 = (m = 2.1, Γ = 0.001))
    restored = constructchannel(old)
    @test (restored.μ12, restored.μ13) == (γ.μ12, γ.μ13)
end

@testset "MeV amplitude boundary and sheet continuation" begin
    ch = MassTestChannel((m1 = 0.2, m2 = 1.6, m3 = 2.0))
    ich = interpolated(ch, 3.9; mstep = 0.001)
    correction(e) = 0.01e
    single, multiple = Amplitude(ich, correction), Amplitude((ich,), correction)
    δm = -0.3592
    for e in (-0.7, 0.2 + 0.1im, -0.3 - 0.05im)
        expected = -(dispersive(ich, e2m(e)) - real(dispersive(ich, e2m(δm)))) + correction(e)
        @test denominator_I(single, e, δm) ≈ expected
        @test denominator_I(multiple, e, δm) ≈ expected
        expected_II = imag(e) >= 0 ? expected : expected - 2im * ρ_thr(ich, e2m(e))
        @test denominator_II(single, e, δm) ≈ expected_II
        @test denominator_II(multiple, e, δm) ≈ expected_II
    end
    @test isapprox(real(denominator_I(Amplitude(ich), δm, δm)), 0; atol = 1e-14)
end

@testset "Photon integration paths retain channel couplings" begin
    ch = γDD((m1 = mγ, m2 = mD⁺, m3 = mD⁰),
        BW(m = mDˣ⁺, Γ = ΓDˣ⁺), BW(m = mDˣ⁰, Γ = ΓDˣ⁰), 2.0, -0.5)
    m = e2m(-0.1 - 0.1im)
    s = m^2
    σ3, σ2 = 4.0 + 0.001im, 4.01 - 0.001im
    ch2 = ChannelWithIntegrationMethod(ch, HookSqrtDalitzMapping{2}())
    ch3 = ChannelWithIntegrationMethod(ch, HookSqrtDalitzMapping{3}())
    f = X2DDpi.decay_matrix_element_squared
    @test f(ch2, s, σ3, σ2) + f(ch3, s, σ3, σ2) ≈ f(ch, s, σ3, σ2)
end

@testset "All historical cache families load" begin
    function checkcaches(x)
        if x isa Dict && haskey(x, "knots") && haskey(x, "coefs")
            restored = interpolated(d2nt(x))
            @test restored.m_cutoff ≈ e2m(x["cutoff"])
            @test restored.itr.coefs == x["coefs"]
            @test all(isfinite, ρ_thr.(Ref(restored), e2m.((-0.3592, 0.0, 5.0))))
        elseif x isa Dict
            foreach(checkcaches, values(x))
        elseif x isa Vector
            foreach(checkcaches, x)
        end
    end
    for file in ("nominal/X_model.json", "nominal/setofmodelg1g2.json", "Dstar_p3/model.json")
        checkcaches(readjson(joinpath(@__DIR__, "..", "results", file)))
    end
end

@testset "Nominal dispersion agrees with the energy-grid formulation" begin
    model = readjson(joinpath(@__DIR__, "..", "results", "nominal", "model.json"))
    for entry in model["ichannels"]
        ich = interpolated(d2nt(entry))
        ev = m2e.(ich.itr.knots[1])
        olditr = X2DDpi.interpolate((ev,), ich.itr.coefs, X2DDpi.Gridded(X2DDpi.Linear()))
        oldrho(e) = e < first(ev) ? 0.0 : e < m2e(ich.m_cutoff) ? olditr(e) :
            X2DDpi.ρ_tb(ich.channel, e2m(e)) * ich.cutoffratio
        function olddisp(e)
            z = e isa Real ? e + 1e-6im : e
            s = e2m(z)^2
            quadgk(t -> oldrho(m2e(sqrt(t))) / t / (t - s), sum(masses(ich.channel))^2, Inf)[1] * s / π
        end
        for e in (-0.3592, 0.2 + 0.05im, -0.3 - 0.05im)
            @test isapprox(dispersive(ich, e2m(e)), olddisp(e); rtol = 2e-6)
        end
        @test isapprox(imag(dispersive(ich, e2m(-0.3592))), ρ_thr(ich, e2m(-0.3592)); rtol = 2e-4)
    end
end

@testset "Single-channel construction with non-nominal masses" begin
    ch = DˣD((m1 = 0.15, m2 = 1.8, m3 = 1.9), BW(m = 2.0, Γ = 0.001))
    ich = interpolated(ch, 4.0; mstep = 0.01)
    @test ρ_thr(ch, sum(masses(ch))) == 0
    @test ρ_thr(ch, sum(masses(ch)) - 0.01) == 0
    @test isfinite(ρ_thr(ch, 3.95))
    @test ρ_thr(ich, 3.95) ≈ ρ_thr(ch, 3.95)
    @test ρ_tb(ch, 3.8) == 0
    @test ρ_tb(ch, 4.0) ≈ sqrt(λ(4.0^2, 2.0^2, 1.9^2)) / 4.0^2
end
