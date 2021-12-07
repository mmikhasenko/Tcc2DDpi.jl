using Pkg
cd(joinpath(@__DIR__,".."))
Pkg.activate(".")
# 
using X2DDpi
using Parameters
using Measurements

using Plots
using LaTeXStrings
theme(:wong2, frame=:box, grid=false, minorticks=true, 
    guidefontvalign=:top, guidefonthalign=:right,
    # xlim=(:auto,:auto), ylim=(:auto,:auto),
    xlab=L"\delta m_{D^0D^0\pi^+}\,\,[\mathrm{MeV}]", lw=1.2, lab="")

ch = πDD((m1=mπ⁺,m2=mπ⁺,m3=mπ⁺), BW(m=0.77,Γ=0.05), BW(m=0.77,Γ=0.05))
# 

const d = ch;
# const s = e2m(-0.1)^2
# 	
function σ3_σ2(x,s,d)
    σ3_0, σ3_e = (d.ms[1]+d.ms[2])^2, (√s-d.ms[3])^2
    σ3 = σ3_0 + x[1]*(σ3_e-σ3_0) # straight path
    # 
    σ2_0, σ2_e = X2DDpi.σ2of3_pm(σ3, d.ms^2, s)
    σ2 = σ2_0 + x[2]*(σ2_e-σ2_0) # straight path
    σ3, σ2
end

# m2e(mD⁰+mD⁰+mπ⁺)
(2mπ⁺)^2

let #e = -1.0-0.01im
    s = (3mπ⁺+0.01-0.06im)^2 #e2m(e)^2
    # 
    σ3_0, _ = σ3_σ2([0.0, 0.1],s,ch)
    σ3_3, _ = σ3_σ2([1.0, 0.1],s,ch)
    #
    plot()
    scatter!([σ3_0, σ3_3], c=1, ms=6)
    # 
    σ2_xv = [σ3_σ2([x, 0.0],s,ch)[2] for x in 0.01:0.01:0.99]
    σ2_yv = [σ3_σ2([x, 1.0],s,ch)[2] for x in 0.01:0.01:0.99]
    # 
    for (x,y) in zip(σ2_xv,σ2_yv)
        plot!([x, y], lab="", c=2)
        scatter!([x, y], lab="", c=2)
    end
    # plot()
    display(VSCodeServer.InlineDisplay(), "image/svg+xml", plot!())
end
plot()


let #e = -1.0-0.01im
    s = (3mπ⁺+0.01-0.043im)^2 #e2m(e)^2
    # 
    σ3_0, _ = σ3_σ2([0.0, 0.1],s,ch)
    σ3_3, _ = σ3_σ2([1.0, 0.1],s,ch)
    # 
    f(σ) = X2DDpi.σ2of3_pm(σ,d.ms^2,s)[1]
    # 
    xv = range(sort(real.([σ3_0, σ3_3]))..., length=300)
    yv = range(sort(imag.([σ3_0, σ3_3]))..., length=300)
    calv = imag.(f.(xv' .+ 1im .* yv))
    heatmap(xv,yv,calv)
    plot!([σ3_0, σ3_3], c=1)
    scatter!([σ3_0, σ3_3], c=1, ms=6)
    # plot()
    display(VSCodeServer.InlineDisplay(), "image/svg+xml", plot!())
    # calv
end

#
# jac = (σ3_e-σ3_0)*(σ2_e-σ2_0)
decay_matrix_element_squared(d,s,σ3,σ2) / (2π*s) * jac


σ2of3_pm(σ,msq,s) = (msq[1]+msq[3]+
	(σ+msq[1]-msq[2])*(s-σ-msq[3])/(2*σ)) .+ [-1,1] .*
		sqrt(λ(σ,msq[1],msq[2])*λ(s,σ,msq[3]))/(2*σ)
