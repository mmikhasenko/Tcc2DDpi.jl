
@userplot DalitzPlot
@recipe function f(hp::DalitzPlot; what2apply=real)
    model, δm = hp.args
    iσx := 2
    iσy := 3
    density = σs->what2apply(
        decay_matrix_element_squared(model,e2m(δm)^2,σs.σ3,σs.σ2))
    (ThreeBodyMasses(m0=e2m(δm0_val), model.ms...), density)
end
