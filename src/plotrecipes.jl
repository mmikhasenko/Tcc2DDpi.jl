
@userplot DalitzPlot
@recipe function f(hp::DalitzPlot; what2apply=real)
    model, δm = hp.args
    iσx := 2
    iσy := 3
    density = σs->what2apply(
        decay_matrix_element_squared(model,e2m(δm)^2,σs.σ3,σs.σ2))
    (ThreeBodyMasses(m0=e2m(δm0_val), model.ms...), density)
end

@userplot DpDzSpectrum
@recipe function f(hp::DpDzSpectrum; what2apply=real)
    xv, y2v, y3v, = hp.args
    
    # 
    n = sum(y2v+y3v) / (length(xv)-1)
    @series begin
        label --> L"\gamma D^+D^0"
        fillrange := 0
        seriescolor := :blue
        seriesalpha := 0.4
        (xv, y3v ./ n)
    end
    @series begin
        label := ""
        linecolor --> :black
        linewidth := 1.5
        (xv, y3v ./ n)
    end
    # 
    @series begin
        label --> L"\pi^0 D^+D^0"
        fillrange := y3v ./ n
        seriescolor := :red
        seriesalpha := 0.4
        (xv, (y2v+y3v) ./ n)
    end
    @series begin
        label := ""
        linecolor --> :black
        linewidth := 1.5
        (xv, (y2v+y3v) ./ n)
    end
    label := ""        
    ()
end
