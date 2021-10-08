using X2DDpi

using Plotly, LaTeXStrings

# get data
fin = readjson(joinpath("results","nominal","logA_2d_grid.json"))
xi = Vector{Float64}(fin["x"])
yi = Vector{Float64}(fin["y"])
zi = Matrix{Float64}(hcat(fin["z"]...))

yp, xp = let ij = findmax(zi)[2]
    yi[ij[1]], xi[ij[2]] 
end

zexts = extrema(zi)
dz = (zexts[2]-zexts[1])
zlims = (zexts[1]-0.2*dz, zexts[2]-0.05*dz)

zcut = zexts[1]+0.05*dz
yc, xc = -81.3e-3/2, 0.0
# let ij = findmin(zi)[2]
#     yi[ij[1]], xi[ij[2]] 
# end


layout = Layout(
    title=L"T_{cc}^+\to D^{0}D^{0}\pi^+\,\,\mathrm{amplitude\,\,in\,\,the\,\,complex\,\,plane.}",
    # autosize=true,
    automargin=true,
    # paper_bgcolor="LightSteelBlue",
    scene_camera_eye=attr(x=0.8, y=-2.0, z=0.4),
    width=650, height=550,
    scene_camera_center = attr(x=-0.1, y=0, z=-0.25),
    scene = attr(xaxis_title="Re √s [MeV]",
                 yaxis_title="Im √s [MeV]",
                 zaxis_title="log |A|²",
                 xaxis_backgroundcolor="rgb(245,245,255)",
                 xaxis_showbackground=true,
                 yaxis_backgroundcolor="rgb(245,245,255)",
                 yaxis_showbackground=true,
                 #  
                 zaxis_range=zlims,
                 #  
                 yaxis_zerolinewidth=4,
                 yaxis_zerolinecolor="red",
            annotations=[
                attr(
                    x=xc,
                    y=yc,
                    z=zcut,
                    text=L"D^{*+}D^0\,\,\mathrm{branch\,\,point}",
                    textangle=0,
                    ax=0,
                    ay=-55,
                    font=attr(
                        color="black",
                        size=12
                    ),
                    arrowcolor="#FDE725",
                    arrowsize=2,
                    arrowwidth=1,
                    arrowhead=2),
                attr(
                    x=xp,
                    y=yp,
                    z=zlims[2],
                    text=L"T_{cc}^{+}\,\,\mathrm{pole}",
                    textangle=0,
                    ax=60,
                    ay=0,
                    font=attr(
                        color="black",
                        size=14
                    ),
                    arrowcolor="#95D840",
                    arrowsize=2,
                    arrowwidth=1,
                    arrowhead=2)],
        )
)
f = plot(surface(
    x=xi,
    y=yi,
    z=zi',
    showscale=false,
    colorscale="Viridis",
    zmin=zexts[1]-0.15dz,
    zmax=zlims[2],
    # reversescale=true,
    contours=attr(
        y=attr(show=true, start=0.0, size=0.1, color="red", project_y=true),
        y_end=0.01,
        z=attr(show=true, start=zlims[1], size=0.3, usecolormap=true, project_z=true),
        z_end=zlims[2]-0.2dz,
    ),
), layout)
# savefig(f, joinpath("plots","nominal","complexplane_model_without_interference.html"))
# savefig(f, joinpath("plots","nominal","complexplane_model_without_interference.png"))
# savefig(f, joinpath("plots","nominal","complexplane_model_without_interference.pdf"))

# Plotly.signin()
remote_plot = post(f)

