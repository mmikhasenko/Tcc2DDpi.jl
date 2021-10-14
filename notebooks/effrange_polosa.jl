### A Pluto.jl notebook ###
# v0.16.1

using Markdown
using InteractiveUtils

# ╔═╡ 8a8ef3e6-1277-407c-b906-7fb173230a7d
using Measurements

# ╔═╡ f1c958c0-ce58-44da-bcd2-678ed39ab858
using Unitful

# ╔═╡ 85b31d10-2c1b-11ec-3920-2b888371d326
mX = -7.18u"MeV"

# ╔═╡ 033a5372-bdd0-4065-8147-7a9a375ed3d4
gLHCb = 0.108

# ╔═╡ a32cb39c-5ad0-457e-8e78-6dc1afe795b8
begin
	const mπ⁰ = 0.1349768u"GeV"
	const mπ⁺ = 0.13957039u"GeV"
	const mD⁰ = 1.86483u"GeV"
	const mD⁺ = 1.86965u"GeV"
	# 
	const mDˣ⁺ = (1.86483+145.4258e-3)*1u"GeV" # m(D) + Δm(D*,D) from PDG
	const mDˣ⁰ = 2.00685u"GeV"
	# 
	const ΓDˣ⁺ = 83.4e-6u"GeV"
	const ΓDˣ⁰ = 55.2e-6u"GeV"
	const mγ = 0.0u"GeV"
	#
	const fm_times_mev = 197.3269804u"fm*MeV"
end

# ╔═╡ 1c80959c-ed0a-44e9-8da1-b16e1d0e3bc0
μ = mD⁰*mDˣ⁰ / (mD⁰+mDˣ⁰)

# ╔═╡ 18c1704e-c2e4-49f2-9c44-e5c7aced7bc3
μ⁺ = mD⁺*mDˣ⁺ / (mD⁺+mDˣ⁺)

# ╔═╡ 9d1aa84b-e5b1-473f-91e2-52d361a5360a
δ = 8.2u"MeV"

# ╔═╡ 6dcc6e60-ad11-48ea-a7fc-3efc6e961898
1u"MeV"/u"GeV"

# ╔═╡ e73dbf75-53d8-4ec4-b23e-3a15abb3a28b
md"""
## Effective range
"""

# ╔═╡ 0e275bfc-9d93-4d02-957b-de9fc77c1964
r0_1, r0_2 =
	-2/(μ*gLHCb),
	-sqrt(μ⁺/(2*μ^2*δ))

# ╔═╡ 782412e1-b150-459a-b973-3e7e27c20408
r0 = uconvert(u"fm", (r0_1 + r0_2) * fm_times_mev)

# ╔═╡ 1eb7633e-4c75-4d18-96d5-bf2fadf8e5cd
md"""
## Scattering length
"""

# ╔═╡ cc0af199-3b96-4f8b-a31a-706bdf716365
inva0_1, inva0_2 = -2mX/gLHCb, -sqrt(2μ⁺*δ)

# ╔═╡ 5f9294d6-68a6-4cde-b59d-4a9790d63668
inva0 = uconvert(u"MeV", inva0_1 + inva0_2)

# ╔═╡ a35aabc0-1e04-4b8e-9883-610209029d1f
md"""
## Compositeness
$1-Z = \sqrt{\frac{1}{1+2|r/a|}}$
"""

# ╔═╡ e6a09750-cce8-4728-8b8f-4fccfa2d258d
1-sqrt(1/(1+abs(2r0*inva0/fm_times_mev)))

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
Measurements = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[compat]
Measurements = "~2.6.0"
Unitful = "~1.9.0"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.7.0-beta4"
manifest_format = "2.0"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Calculus]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f641eb0a4f00c343bbc32346e1217b86f3ce9dad"
uuid = "49dc2e85-a5d0-5ad3-a950-438e2897f1b9"
version = "0.5.1"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f74e9d5388b8620b4cee35d4c5a618dd4dc547f4"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.3.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.Measurements]]
deps = ["Calculus", "LinearAlgebra", "Printf", "RecipesBase", "Requires"]
git-tree-sha1 = "31c8c0569b914111c94dd31149265ed47c238c5b"
uuid = "eff96d63-e80a-5855-80a2-b1b0885c5ab7"
version = "2.6.0"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Random]]
deps = ["Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
git-tree-sha1 = "44a75aa7a527910ee3d1751d1f0e4148698add9e"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.1.2"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "4036a3bd08ac7e968e27c203d45f5fff15020621"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.1.3"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.Unitful]]
deps = ["ConstructionBase", "Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "a981a8ef8714cba2fd9780b22fd7a469e7aaf56d"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.9.0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
"""

# ╔═╡ Cell order:
# ╠═8a8ef3e6-1277-407c-b906-7fb173230a7d
# ╠═f1c958c0-ce58-44da-bcd2-678ed39ab858
# ╠═85b31d10-2c1b-11ec-3920-2b888371d326
# ╠═033a5372-bdd0-4065-8147-7a9a375ed3d4
# ╠═a32cb39c-5ad0-457e-8e78-6dc1afe795b8
# ╠═1c80959c-ed0a-44e9-8da1-b16e1d0e3bc0
# ╠═18c1704e-c2e4-49f2-9c44-e5c7aced7bc3
# ╠═9d1aa84b-e5b1-473f-91e2-52d361a5360a
# ╠═6dcc6e60-ad11-48ea-a7fc-3efc6e961898
# ╟─e73dbf75-53d8-4ec4-b23e-3a15abb3a28b
# ╠═0e275bfc-9d93-4d02-957b-de9fc77c1964
# ╠═782412e1-b150-459a-b973-3e7e27c20408
# ╟─1eb7633e-4c75-4d18-96d5-bf2fadf8e5cd
# ╠═cc0af199-3b96-4f8b-a31a-706bdf716365
# ╠═5f9294d6-68a6-4cde-b59d-4a9790d63668
# ╟─a35aabc0-1e04-4b8e-9883-610209029d1f
# ╠═e6a09750-cce8-4728-8b8f-4fccfa2d258d
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
