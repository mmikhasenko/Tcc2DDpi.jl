### A Pluto.jl notebook ###
# v0.16.1

using Markdown
using InteractiveUtils

# ╔═╡ 61b01d0b-1f8f-4e9d-805a-b533cd39410c
begin
    import Pkg
    Pkg.activate(@__DIR__)
end

# ╔═╡ bb940830-336f-11ec-1c51-716442c1f7cf
using Symbolics

# ╔═╡ 8bb045d3-c242-4c43-bdc4-bf3af643c342
using LinearAlgebra

# ╔═╡ 3410c558-a6af-4e84-a9f4-e6e4efda8462
md"""
# K-matrix with non-diagonal G
using `Symbolics`
"""

# ╔═╡ 40e5fc56-7ce1-4b64-89f6-bb4a114a35e2
begin
	@variables G_11 G_12 G_22
	@variables g_1 g_2
end

# ╔═╡ d36aac19-bf1b-4db6-9563-14b200db2fb0
begin
	G = [G_11 0
		 0 G_22]
	#
	K = [g_1^2 g_1*g_2
		 g_1*g_2 g_2^2]
end

# ╔═╡ 305a3f35-40fc-4634-9304-a8ca9532bc86
begin
	D = [1 0; 0 1] - G*K
	invD = simplify.(inv(D))
	T = K*invD
end;

# ╔═╡ b7a6b157-ca85-4f1e-8bcf-4c5569dc0981
simplify_fractions(simplify(inv(D)[1,1]))

# ╔═╡ 95b92ae9-14e4-4f4d-96fb-be5671aa85bc
invD

# ╔═╡ 310cf756-d638-4a4d-8b46-2695982dc213
simplify_fractions(T[1,2])

# ╔═╡ c9f8da5c-1d32-43ce-8ea0-299849412ea5
T[1,1]

# ╔═╡ Cell order:
# ╠═61b01d0b-1f8f-4e9d-805a-b533cd39410c
# ╟─3410c558-a6af-4e84-a9f4-e6e4efda8462
# ╠═bb940830-336f-11ec-1c51-716442c1f7cf
# ╠═40e5fc56-7ce1-4b64-89f6-bb4a114a35e2
# ╠═d36aac19-bf1b-4db6-9563-14b200db2fb0
# ╠═8bb045d3-c242-4c43-bdc4-bf3af643c342
# ╠═b7a6b157-ca85-4f1e-8bcf-4c5569dc0981
# ╠═305a3f35-40fc-4634-9304-a8ca9532bc86
# ╠═95b92ae9-14e4-4f4d-96fb-be5671aa85bc
# ╠═310cf756-d638-4a4d-8b46-2695982dc213
# ╠═c9f8da5c-1d32-43ce-8ea0-299849412ea5
