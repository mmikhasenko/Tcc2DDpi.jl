using Plots

function simpleBW(x,m=1.01,Γ=0.01)
    xth = 1.0
    ρ = sqrt(x-xth)
    m*Γ/(m^2-x^2-1im*m*Γ*ρ)
end

let
    plot()
    for Γ in exp.(range(log(0.01), log(0.5), length=10))
        plot!(x->abs2(simpleBW(x,1.01,Γ)), 1, 1.03)
    end
    plot!()
end