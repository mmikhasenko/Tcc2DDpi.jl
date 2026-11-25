### A Pluto.jl notebook ###
# v0.16.1

using Markdown
using InteractiveUtils

# ╔═╡ 61b01d0b-1f8f-4e9d-805a-b533cd39410c
begin
    import Pkg
    Pkg.activate(@__DIR__)
end

# ╔═╡ d433d680-3375-11ec-26dd-c359f027b37e
using SymPy

# ╔═╡ 4bcd6f88-d35a-41c2-aa15-915f0e06b721
md"""
# K-matrix with non-diagonal G
using `SymPy`
"""

# ╔═╡ 0499f433-9cd0-475f-9974-e0c041ebde1d
begin
	@syms G_11::real G_12::real G_22::real
	@syms g_1::real g_2::real
end

# ╔═╡ 9112eb24-62e4-40dd-a6a5-e942f5ccda91
begin
	G = [G_11 0
		 0 G_22]
	#
	K = [g_1^2 g_1*g_2
		 g_1*g_2 g_2^2]
end

# ╔═╡ 348f8ea1-e675-46ff-8a42-7e987e468746
begin
	𝕀₂ = [1 0; 0 1]
	D = 𝕀₂ - G*K
	invD = inv(D)
	T = K*invD
end;

# ╔═╡ 15af2ecd-9823-4fa2-992c-edbc9e2666fa
simplify.(T)

# ╔═╡ c22aaa1c-d390-4493-ad47-f1063953742e
simplify.([g_1 g_2] * invD)

# ╔═╡ 48cc47f5-cf11-47a8-a2a4-010cbdf71f6f
md"""
### Expression in the paper
$\begin{align}
    \mathcal{A}_{\mathrm{U}} \begin{pmatrix} \phantom{-} g\\- g  \end{pmatrix}
= [1- K G]^{-1}P\,.
\end{align}$
"""

# ╔═╡ 2312e552-1cbf-4b74-8014-08d30e337833
begin
	P = transpose([g_1 g_2])
	A_UP = inv(𝕀₂ - K*G)*P
	simplify.(A_UP)
end

# ╔═╡ Cell order:
# ╠═61b01d0b-1f8f-4e9d-805a-b533cd39410c
# ╟─4bcd6f88-d35a-41c2-aa15-915f0e06b721
# ╠═d433d680-3375-11ec-26dd-c359f027b37e
# ╠═0499f433-9cd0-475f-9974-e0c041ebde1d
# ╠═9112eb24-62e4-40dd-a6a5-e942f5ccda91
# ╠═348f8ea1-e675-46ff-8a42-7e987e468746
# ╠═15af2ecd-9823-4fa2-992c-edbc9e2666fa
# ╠═c22aaa1c-d390-4493-ad47-f1063953742e
# ╟─48cc47f5-cf11-47a8-a2a4-010cbdf71f6f
# ╠═2312e552-1cbf-4b74-8014-08d30e337833
