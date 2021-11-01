using Plots

gr()

dy = 0.9
r = 0.5
loop_point = 0.334π
δ = 0.06

lDˣ = 5
lD = 2
lπ = 2

ylim = (-1,dy+r+0.1)

p1 = let 
    plot(xaxis=false, xticks=false,
         yaxis=false, yticks=false, size=(200,150))
    plot!(cis.(π-loop_point:0.01:2π+loop_point), l=(lD,:black), ylim = ylim)#
    plot!(cis.(0:0.01:loop_point), l=(lDˣ,:black))#
    plot!(cis.(π-loop_point:0.01:π), l=(lDˣ,:black))#
    scatter!([-1+0im, 1+0im], lab="", m=(:black, 6))
    scatter!([cis(loop_point), cis(π-loop_point)], lab="", m=(:black, 5))
    plot!(dy*1im .+ r .* cis.(0:0.01:π), l=(lD,:black))
    plot!(dy*1im .+ r .* cis.(π:0.01:2π), l=(lπ,:black, :dash))
    # 
    plot!([-1.1,-1] .+   0im, l=(2,:black))
    plot!([-1.1,-1] .+ δ*1im, l=(2,:black))
    plot!([-1.1,-1] .- δ*1im, l=(2,:black))
    plot!([+1.1,+1] .+   0im, l=(2,:black))
    plot!([+1.1,+1] .+ δ*1im, l=(2,:black))
    plot!([+1.1,+1] .- δ*1im, l=(2,:black))
    # plot!([-1.2im, (dy+r+0.2)*1im], l=(1,:red))
    # annotate!([
    #     (0,-1+δ,text(L"D", :bottom)),
    #     (0,dy+r-δ,text(L"D", :top)),
    #     (0,dy-r-δ,text(L"\pi\,\,\mathrm{or}\,\,\gamma", :top)),
    #     (reim(cis(loop_point/2)+δ)...,text(L"D^*", :bottom, :left)),
    #     (reim(cis(π-loop_point/2)+δ)...,text(L"D^*", :bottom, :right)),
    # ])
    plot!(xlab="", ylab="")
end
savefig(joinpath("plots","directloop_nolab.pdf"))

p2 = let dy = 0.9, r = 0.5, loop_point = 0.334π, δ=0.06, lDˣ=5, lD=2, lπ=2
    plot(xaxis=false, xticks=false,
         yaxis=false, yticks=false, size=(200,150))
    plot!(cis.(0:0.01:2π), l=(lD,:black), xlim=(-1.1,1.1))#
    plot!(cis.(-loop_point:0.01:0), l=(lDˣ,:black))#
    plot!(cis.(π-loop_point:0.01:π), l=(lDˣ,:black))#
    scatter!([-1+0im, 1+0im], lab="", m=(:black, 6))
    scatter!([cis(-loop_point), cis(π-loop_point)], lab="", m=(:black, 5))
    plot!([cis(π-loop_point), cis(-loop_point)], l=(lπ,:black,:dash))#
    # 
    plot!([-1.1,-1] .+   0im, l=(2,:black))
    plot!([-1.1,-1] .+ δ*1im, l=(2,:black))
    plot!([-1.1,-1] .- δ*1im, l=(2,:black))
    plot!([+1.1,+1] .+   0im, l=(2,:black))
    plot!([+1.1,+1] .+ δ*1im, l=(2,:black))
    plot!([+1.1,+1] .- δ*1im, l=(2,:black))
    #
    # plot!([-1.2im, 1.2im], l=(1,:red))
    # annotate!([
    #     (reim(cis(loop_point-π)+δ)...,text(L"D", :bottom, :left)),
    #     (reim(cis(loop_point)-δ)...,text(L"D", :top, :right)),
    #     (2δ, 0,text(L"\pi\,\,\mathrm{or}\,\,\gamma", :left)),
    #     (reim(cis(-loop_point/2)-δ)...,text(L"D^*", :top, :left)),
    #     (reim(cis(π-loop_point/2)+δ)...,text(L"D^*", :bottom, :right)),
    # ])
    plot!(xlab="", ylab="", ylim = ylim)
end
savefig(joinpath("plots","crossedloop_nolab.pdf"))

plot(p1,p2,layout=grid(1,2), size=(600,200), link=:y)


m2e(mD⁺+mD⁺+mπ⁺)

m2e(mD⁰+mD⁰+mπ⁺)