using Plots
import Plots.PlotMeasures:mm
using LaTeXStrings
theme(:wong2, size=(500,350), minorticks=true, grid=false, frame=:box,
    guidefontvalign=:top, guidefonthalign=:right,
    foreground_color_legend = nothing,
    legendfontsize=9, legend =:topright,
    xlim=(:auto,:auto), ylim=(:auto,:auto),
    xlab=L"\cos\theta")
# 
using SymPy
using SymPy.PyCall
PyCall.pyimport_conda("sympy.physics.wigner",       "sympy")
PyCall.pyimport_conda("sympy.physics.optics",       "sympy")
PyCall.pyimport_conda("sympy.physics.quantum.spin", "sympy")
#
import_from(sympy.physics.quantum.spin)
import_from(sympy.physics.quantum.spin, (:WignerD,), typ=:Any)
import_from(sympy.physics.wigner)

clgn(j1,λ1,j2,λ2,j,λ) = clebsch_gordan(Sym(j1),Sym(j2),j,λ1,λ2,λ)
wignerd(j,λ1,λ2,θ) = WignerD(j,λ1,λ2,0,θ,0).doit()

θ, z = @vars θ z real=true

A(λ; p) = sqrt((2p.L+1)*(2p.j+1)/(2p.j0+1))*clgn(p.L,0,p.j,λ,p.j0,λ)*wignerd(p.j,λ,0,θ)
I(; p) = simplify(sum(abs2(A(λ; p)) for λ in -p.j:p.j))

Iz(; p) = I(; p).subs(cos(θ),z).subs(sin(θ),sqrt(1-z^2))

@show I(; p=(j=1,L=1,j0=1))
@show Iz(; p=(j=1,L=1,j0=1))

Iz(; p=(L=0,j=1,j0=1))

function L_lab(;p)
    l = ('S','P','D')[p.j+1]
    L = ('S','P','D')[p.L+1]
    parity = isodd(p.L+p.j+1) ? '-' : '+'
    latexstring("(D^0\\pi^+)_{\\mathrm{$l}}D^0\\,\\,$(L)\\textrm{-}\\mathrm{wave}: J^P=$(p.j0)^$(parity)")
end

p1 = let
    plot(ylim=(0,5))
    p=(j=1,L=0,j0=1); plot!(Iz(; p)+1e-6z, -1,1, lab=L_lab(;p))
end
# 
p2 = let
    plot(leg=:topright)
    p=(j=1,L=1,j0=0); plot!(Iz(; p), -1,1, lab=L_lab(;p))
    p=(j=1,L=1,j0=1); plot!(Iz(; p), -1,1, lab=L_lab(;p))
    p=(j=1,L=1,j0=2); plot!(Iz(; p), -1,1, lab=L_lab(;p))
end

p3 = let
    plot(leg=:topright)
    p=(j=1,L=2,j0=1); plot!(Iz(; p), -1,1, lab=L_lab(;p))
    p=(j=1,L=2,j0=2); plot!(Iz(; p), -1,1, lab=L_lab(;p))
    p=(j=1,L=2,j0=3); plot!(Iz(; p), -1,1, lab=L_lab(;p))
end


plot(p1,p2,p3, layout=grid(3,1), size=(400,900), left_margin=4mm)
savefig(joinpath("plots","nominal","angular_distributions.pdf"))
