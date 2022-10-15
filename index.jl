### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# ╔═╡ 176dc896-4012-44d1-812f-a8041889959d
# ╠═╡ show_logs = false
begin
	cd(@__DIR__)
	import Pkg
	Pkg.activate(".")
	Pkg.instantiate()
	# 
	using PlutoSliderServer
	using Dates
	using Markdown
end

# ╔═╡ 2d94ee0d-df8e-4643-859f-dccc770362e4
md"""
# Studies of the $T_{cc}^+$ tetraquark
"""

# ╔═╡ 5b1ed9c1-32c1-44c3-b1c4-848f2e3b04e9
let
	notebooklist = PlutoSliderServer.find_notebook_files_recursive(@__DIR__)
	sort!(notebooklist)
	filter!(x->x!="index.jl", notebooklist)
	# 
	content = "## Table of Notebooks\n"
	for (i,name) in enumerate(notebooklist)
		_name = PlutoSliderServer.without_pluto_file_extension(name)
		_name_html = _name * ".html"
		content *= " $(i). [**`$(_name)`**]($(_name_html))\n"
	end
	Markdown.parse(content)
end

# ╔═╡ Cell order:
# ╟─2d94ee0d-df8e-4643-859f-dccc770362e4
# ╟─5b1ed9c1-32c1-44c3-b1c4-848f2e3b04e9
# ╟─176dc896-4012-44d1-812f-a8041889959d
